library(rTensor)
library(MASS)
library(glmnet)
rm(list=ls())
set.seed(1)

#first get all data
aff.path="D:\\desktop\\simulation new\\data set\\FGNET\\FGNET\\small29n"
setwd(aff.path)
files <- list.files(path = aff.path, pattern = "*.csv" )
n=length(files)
P=29
p=20
n1=floor(n*9/10)
XX=matrix(rep(0),n1,P*P)
Y=matrix(rep(0),n1,1)
YY=matrix(rep(0),n1,1)
i=1
for (i in 1:n1)
{
  XX[i,]=matrix(t(as.matrix(read.csv(files[i],header = F))),1,P^2)
  Y[i]=as.numeric(substr(files[i],5,6))
}

#########
data0=data.frame(XX,Y)

beta=betaols=glm(Y ~ XX-1, data=data0,family='gaussian')[1]$coefficients
beta[which(is.na(beta)==TRUE)]=mean(beta[which(is.na(beta)==F)])
betaols=beta
beta_ols=as.tensor(matrix(betaols,P,P))
beta_ols_dec=tucker(beta_ols,ranks=c(p,p) ,tol=1e-6,max_iter = 50)
beta_ols_eval=as.numeric(beta_ols_dec$Z@data)

ne=floor(p^2-10)
m=1
ne2=0
multi=1
lambda2=8
sprs=F
lambda=50/multi
iter=10000
betamat=NULL
betamat1=NULL
eigvalmat=NULL
eigvalmat1=NULL
ll=matrix(rep(0),iter,1)
ll_bar=matrix(rep(0),iter,1)
stop=0.01
beta2=betaols

i=1#i+1
for(i in 1:(iter)){
  noise=NULL
  beta_1=as.tensor(matrix(beta,P,P))
  beta_1_dec=tucker(beta_1,ranks=c(p,p) ,tol=1e-6,max_iter = 20)
  eigval1=beta_1_dec$Z@data
  eigvec=beta_1_dec$U
  
  eigvalmat=rbind(eigvalmat,eigval1)
  if(i>600){
    sign=apply(sign( eigvalmat[(i-600):i,]),2,median)
    eigval2=matrix(apply( abs(eigvalmat[(i-600):i,]),2,median),p^2,m)
    eigval=eigval2*sign
  }else{
    eigval=(eigval1)
  }
  eigvalmat1=rbind(eigvalmat1,as.numeric(eigval))
  smalllizt <- list('mat1' = matrix(eigvec[[2]],P,p),
                    'mat2' = matrix(eigvec[[1]],P,p))
  beta_k=kronecker_list(smalllizt)
  eigvec_expand=beta_k
  ## set some numerical thresholds
  thred1=1e-6
  thred2=1e+6
  eigval=eigval*(abs(eigval)>thred1)+thred1*(abs(eigval)<=thred1) 
  ## generate noises
  for(j in 1:ne){
    linrad=rnorm(p^2,0,sqrt(lambda/abs(as.numeric(eigval))^4))
    linrad=linrad*(abs(linrad)<thred2)+thred2*(abs(linrad)>=thred2)*sign(linrad)
    lincomb=matrix((linrad),p^2,1)
    
    noise_j=matrix(matrix(eigvec_expand%*%lincomb,P,P),1,P^2)
    noise=rbind(noise,noise_j)
  }
  
  noisex=NULL
  noisey=NULL
  for(k in 1:multi){
    noisex=rbind(noisex,noise)
    noisey=rbind(noisey,matrix(0,ne,m))
  }
  xe=rbind(XX,noisex)
  ye=rbind(Y,noisey)
  
  order=sample(1:nrow(xe),nrow(xe))
  xe=xe[order,]
  ye=ye[order]
  
  datae=data.frame(xe,ye)
  rg=lm(ye ~ xe-1, data=datae)
  beta1=rg[1]$coefficients
  
  ## record the beta estimation
  betamat=rbind(betamat,as.numeric(beta1))
  betaold=beta
  ll[i]=logLik(rg)
  ## moving average
  if(i>600){
    beta=matrix(apply(betamat[(i-600):i,],2,median),P^2,m)
    ll_bar[i]=median(ll[(i-600):i,])
  }else{
    if(i==iter){
      beta=matrix(apply(betamat[(i-600):i,],2,median),P^2,m)
      ll_bar[i]=median(ll[(i-600):i,])
    }
    else{
      beta=beta1
      ll_bar[i]=ll[i,]
    }
  }
  betamat1=rbind(betamat1,as.numeric(beta))
  
  if(!i%%400){print(i)}
  converge=0
  if (i>4000)
  {
    if (abs(ll_bar[i]-ll_bar[i-1])<stop)
    #if (sum(abs(betamat1[i,]-betamat1[i-1,]))<0.001)
    {
      converge=1
      break
    }
  }
}

#Doing prediction
setwd(aff.path)
XX_p=matrix(rep(0),n-n1,P*P)
Y_p=matrix(rep(0),n-n1,1)
for (i in 1:(n-n1))
{
  XX_p[i,]=matrix(t(as.matrix(read.csv(files[n1+i],header = F))),1,P^2)
  Y_p[i]=as.numeric(substr(files[n1+i],5,6))
}
diff=Y_p-XX_p%*%matrix(beta,P*P,1)
diff_ols=Y_p-XX_p%*%matrix(betaols,P*P,1)
XX%*%matrix(betaols,P*P,1)
sum(abs(diff))/(n-n1)
sum(abs(diff_ols))/(n-n1)

# save(beta,file="beta_our.Rdata")
# save(diff,file="diff_our.Rdata")
# save(betaols,file="beta_ols.Rdata")
# save(diff_ols,file="diff_ols.Rdata")
