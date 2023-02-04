library(rTensor)
library(MASS)
library(glmnet)
rm(list=ls())
set.seed(1)

#first get all data
aff.path="D:\\desktop\\simulation new\\data set\\pointing04\\small"
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
  name=strsplit(files[i], split=".",fixed=T)[[1]][1]
  sp=unlist(strsplit(name,split="(?=[+-])",perl=T))
  Y[i]=as.numeric(paste(sp[4],sp[5],collapse="",sep=""))
}

U=array(rnorm(P*p*2,0,1),dim=c(2,P,p))
G=as.tensor(array(rnorm(p^2,0,1),dim=c(p,p)))
stop=0.01
ite=4000
jump=0
lam=10^(-10)

ll=matrix(rep(0),ite,1)
betamat=NULL
i=1
for (i in 1:ite)
{
  for (j in 1:2)
  {
    xx=matrix(rep(0),n1,P*p)#is the data that can be used as GLM reg to update Uj
    if(j==1)
    {
      #using tucker decom to turn optimization of Uj into a GLM problem
      kp=matrix(U[2,,],P,p)
      Gn=k_unfold(G,1)  
      Gn=(Gn@data)
      for (k in 1:n1)
      {
        Xn=k_unfold(as.tensor(array(XX[k,],dim=c(P,P))),1)  
        Xn=(Xn@data)
        xx[k,]=matrix(Xn%*%kp%*%t(Gn),1,P*p)
      }
      datanew=data.frame(xx,Y)
      U[1,,]=lm(Y ~ xx-1, data=datanew)[1]$coefficients
    }else if (j==2){
      kp=matrix(U[1,,],P,p)
      Gn=k_unfold(G,2)  
      Gn=(Gn@data)
      for (k in 1:n1)
      {
        Xn=k_unfold(as.tensor(array(XX[k,],dim=c(P,P))),2)  
        Xn=(Xn@data)
        xx[k,]=matrix(Xn%*%kp%*%t(Gn),1,P*p)
      }
      datanew=data.frame(xx,Y)
      U[2,,]=lm(Y ~ xx-1, data=datanew)[1]$coefficients
    }
  }
  #using tucker decom to optimize core tensor
  xx=matrix(rep(0),n1,p*p)#is the data that can be used as GLM reg to update Uj
  L=list('mat1'=matrix(U[2,,],P,p),'mat2'=matrix(U[1,,],P,p))
  kp=kronecker_list(L)
  k=1
  for (k in 1:n1)
  {
    Xn=vec(as.tensor(array(XX[k,],dim=c(P,P))))
    xx[k,]=t(t(kp)%*%Xn)
  }
  datanew=data.frame(xx,Y)
  reg=lm(Y ~ xx-1, data=datanew)
  G=as.tensor(array(reg$coefficients,dim=c(p,p)))
  
  print(i)
  beta=kp%*%matrix(G@data,p^2,1)
  betamat=rbind(betamat,as.numeric(beta))
  ll[i,]=logLik(reg)
  convergence=0
  if (i>1)
  {
    if (sum(abs(betamat[i,]-betamat[i-1,]))<stop)
    #if ((ll[i,1]-ll[i-1,1])<stop)
    {
      beta=kp%*%matrix(G@data,p^2,1)
      convergence=1
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
  name=strsplit(files[n1+i], split=".",fixed=T)[[1]][1]
  sp=unlist(strsplit(name,split="(?=[+-])",perl=T))
  Y_p[i]=as.numeric(paste(sp[4],sp[5],collapse="",sep=""))
}
diff=Y_p-XX_p%*%matrix(beta,P*P,1)
sum(abs(diff))/(n-n1)

save(beta,file="beta_tuc.Rdata")
save(diff,file="diff_tuc.Rdata")