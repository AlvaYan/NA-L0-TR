library(rTensor)
library(MASS)
#library(extraDistr)
library(glmnet)
rm(list=ls())
set.seed(1)

#---------------------first get all data
aff.path="D:\\desktop\\simulation new\\data set\\FGNET\\FGNET\\small29n"
setwd(aff.path)
files <- list.files(path = aff.path, pattern = "*.csv" )
n=length(files)
P=29
n1=floor(n*9/10)# n1-n2:validation set
n2=floor(n*8/10)# 0-n2:training set
#X=array(rep(0),dim=c(n,P,P))
XX=matrix(rep(0),n2-1,P*P)
Y=matrix(rep(0),n2-1,1)
i=1
for (i in 1:(n2-1))
{
  #X[i,,]=t(as.matrix(read.csv(files[i],header = F)))
  XX[i,]=matrix(t(as.matrix(read.csv(files[i],header = F))),1,P^2)
  Y[i]=as.numeric(substr(files[i],5,6))
}

XX_v=matrix(rep(0),n1-n2,P*P)
Y_v=matrix(rep(0),n1-n2,1)
for (i in n2:(n1-1))
{
  #X[i,,]=t(as.matrix(read.csv(files[i],header = F)))
  XX_v[i-n2+1,]=matrix(t(as.matrix(read.csv(files[i],header = F))),1,P^2)
  Y_v[i-n2+1]=as.numeric(substr(files[i],5,6))
}

#-----------------------
XX_t=rbind(XX,XX_v)
Y_t=rbind(Y,Y_v)
data0=data.frame(XX_t,Y_t)
beta=betaols=glm(Y_t ~ XX_t-1, data=data0,family='gaussian')[1]$coefficients
beta[which(is.na(beta)==TRUE)]=mean(beta[which(is.na(beta)==F)])
betaols=beta
#beta_ols=as.tensor(matrix(betaols,P,P))
#beta_ols_dec=tucker(beta_ols,ranks=c(p,p) ,tol=1e-6,max_iter = 50)
#beta_ols_eval=as.numeric(beta_ols_dec$Z@data)

#-------define model
our=function(betaols,p, ne){
  #ne=floor(p^2-1)
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
  beta=betaols
  i=1#i+1
  for(i in 1:(iter)){
    noise=NULL
    
    ## beta is stored in 2d matrix, therefore to use it as a tensor, we need to fold it first
    beta_1=as.tensor(matrix(beta,P,P))
    ## apply tucker decomposition to get the core tensor and unitary matrices
    #set.seed(1)
    beta_1_dec=tucker(beta_1,ranks=c(p,p) ,tol=1e-6,max_iter = 20)
    ## treat the core tnsor as eigen values and unitary matrices as eigenvectors
    eigval1=beta_1_dec$Z@data
    eigvec=beta_1_dec$U
    
    eigvalmat=rbind(eigvalmat,eigval1)
    
    if(i>600){
      sign=apply(sign( eigvalmat[(i-600):i,]),2,median)
      eigval2=matrix(apply(abs(eigvalmat[(i-600):i,]),2,median),p^2,m)
      #sign2=sign(eigval2)
      eigval=eigval2*sign
    }else{
      #sign=sign(beta_0_eval)
      #sign2=sign(eigval1)
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
    
    ## do linear regression on vectorized tensor data
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
    #ll[i]=sum(dnorm(ye,mean=xe%*%matrix(beta,p^3,1),sd=sigma,log=TRUE))
    #ll[i]=sum((matrix(ye,n+ne,1)-matrix(xe,n+ne,p^3)%*%matrix(beta,p^3,1))^2)
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
  return(list(beta,converge,ll[i]))
}
#-------------------validation
parameters=list(c(29,840),c(29,831),c(29,821),c(25,623),c(25,615),c(25,605),c(20,399),c(20,390),c(20,380))
beta_all=NULL
err_all=NULL
i=1
i=i+1

i=8
for (i in 1:length(parameters)){
  p=unlist(parameters[i])[1]
  ne=unlist(parameters[i])[2]
  res=our(betaols,p,ne)
  cur_beta=matrix(unlist(res[1]),P*P,1)
  beta_all=rbind(beta_all,unlist(res[1]))
  
  diff=Y_v-XX_v%*%matrix(cur_beta,P*P,1)
  err=sum(abs(diff))/(n1-n2)
  err_all=cbind(err_all,sum(abs(diff))/(n1-n2))
}
#-------------------Doing prediction
setwd(aff.path)
beta=which.min(err_all)
beta=matrix(beta,P*P,1)
XX_p=matrix(rep(0),n-n1+1,P*P)
Y_p=matrix(rep(0),n-n1+1,1)
for (i in n1:n)
{
  #X[i,,]=t(as.matrix(read.csv(files[i],header = F)))
  XX_p[i-n1+1,]=matrix(t(as.matrix(read.csv(files[i],header = F))),1,P^2)
  Y_p[i-n1+1]=as.numeric(substr(files[i],5,6))
}
diff=Y_p-XX_p%*%matrix(beta,P*P,1)
diff_ols=Y_p-XX_p%*%matrix(betaols,P*P,1)
XX%*%matrix(betaols,P*P,1)
sum(abs(diff))/(n-n1)
sum(abs(diff_ols))/(n-n1)
print(parameters)
#sum(abs(Y-XX%*%matrix(betaols,P*P,1)))/n1
setwd("D:\\desktop\\simulation new\\data set\\FGNET\\FGNET\\redo2")
save(beta,file="beta_our.Rdata")
save(diff,file="diff_our.Rdata")
save(betaols,file="beta_ols.Rdata")
save(diff_ols,file="diff_ols.Rdata")
save(cur_beta,file="beta_vali8.Rdata")
save(err,file="err_vali8.Rdata")
# 5 is the best, 25
#############################################
######## trials
#aff.path="/afs/crc.nd.edu/user/t/tyan/R"
#setwd(aff.path)
#files <- list.files(path = aff.path, pattern = "*.csv" )
#i=1
#P=320
#w=matrix(t(as.matrix(read.csv(files[i],header = F))),1,P^2)
#save(w[1,],file="trial.Rdata")
#w=matrix(abs(beta),P,P)
#w[which(w<0.35)]=0
#image(w)

#6.031962, p=29, ne=floor(p^2-10), iter=4000, small29n, lamda=50
#5.955997, p=20, ne=floor(p^2-10), iter=4000, small29n, lamda=50
#6.478786, p=15, ne=floor(p^2-10), iter=4000, small29n, lamda=50
#4.658921, p=20, ne=floor(p^2-1), iter=4000, small29n, lamda=50
#4.66273, p=20, ne=floor(p^2-1), iter=7000, small29n, lamda=50
#4.9920, p=20, ne=floor(p^2-1), iter=4000, small29n, lamda=500
#6.648552, p=20, ne=floor(p^2-1), iter=4000, 29p, lamda=50
#6.392828, p=20, ne=floor(p^2-10), iter=4000, 29p, lamda=50
#5.782411, p=20, ne=floor(p^2-10), iter=4000, small29n, lamda=50
#10.83 , p=20, ne=floor(p^2-1), iter=4000, small29n, lamda=50
#3.9727 , p=20, ne=floor(p^2-10), iter=3000, small29n, lamda=50
#4.574104 , p=20, ne=floor(p^2-1), iter=3000, small29n, lamda=50