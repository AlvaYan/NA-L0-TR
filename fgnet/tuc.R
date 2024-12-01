library(rTensor)
library(MASS)
#library(extraDistr)
library(glmnet)
rm(list=ls())
set.seed(1)

#first get all data
aff.path="D:\\desktop\\simulation new\\data set\\FGNET\\FGNET\\small29n"
setwd(aff.path)
files <- list.files(path = aff.path, pattern = "*.csv" )
n=length(files)
P=29
p=25
n1=floor(n*9/10)#length(files)
XX=matrix(rep(0),n1,P*P)
Y=matrix(rep(0),n1,1)
YY=matrix(rep(0),n1,1)
i=1
for (i in 1:n1)
{
  #X[i,,]=t(as.matrix(read.csv(files[i],header = F)))
  XX[i,]=matrix(t(as.matrix(read.csv(files[i],header = F))),1,P^2)
  Y[i]=as.numeric(substr(files[i],5,6))
}

U=array(rnorm(P*p*2,0,1),dim=c(2,P,p))
G=as.tensor(array(rnorm(p^2,0,1),dim=c(p,p)))
stop=0.01
ite=5000
jump=0
lam=10^(-10)
betamat=NULL
ll=matrix(rep(0),ite,1)#llh function
i=1
i=i+1
for (i in 1:ite)
{
  #update each U(egvector), according to Tucker Tensor Regression and Neuroimaging Analysis
  #j=2
  for (j in 1:2)
  {
    xx=matrix(rep(0),n1,P*p)#is the data that can be used as GLM reg to update Uj
    if(j==1)
    {
      #using tucker decom to turn optimization of Uj into a GLM problem
      #L=list('mat1'=matrix(U[3,,],p,p),'mat2'=matrix(U[2,,],p,p))
      #kp=kronecker_list(L)
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
      #Gn=t(Gn)
    }else if (j==2){
      
      #using tucker decom to turn optimization of Uj into a GLM problem
      #L=list('mat1'=matrix(U[3,,],p,p),'mat2'=matrix(U[1,,],p,p))
      #kp=kronecker_list(L)
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
    #Xn=(Xn@data)
    xx[k,]=t(t(kp)%*%Xn)
  }
  #datanew=data.frame(xx,y)
  #next minimize regularized loss function
  #first select best lambda
  datanew=data.frame(xx,Y)
  reg=lm(Y ~ xx-1, data=datanew)
  G=as.tensor(array(reg$coefficients,dim=c(p,p)))
  
  print(i)
  beta=kp%*%matrix(G@data,p^2,1)
  betamat=rbind(betamat,as.numeric(beta))
  #ll[i,1]=sum(dnorm(Y,mean=mean,sd=sigma,log=TRUE))
  ll[i,]=logLik(reg)
  convergence=0
  if (i>1)
  {
    #if ((ll[i,1]-ll[i-1,1])<stop)
    if (sum(abs(betamat[i,]-betamat[i-1,]))<stop)
    {
      beta=kp%*%matrix(G@data,p^2,1)
      convergence=1
      break
    }
  }
}
#return(list(beta,convergence,ll[i,1]))
#Doing prediction
setwd(aff.path)
XX_p=matrix(rep(0),n-n1,P*P)
Y_p=matrix(rep(0),n-n1,1)
for (i in 1:(n-n1))
{
  #X[i,,]=t(as.matrix(read.csv(files[i],header = F)))
  XX_p[i,]=matrix(t(as.matrix(read.csv(files[n1+i],header = F))),1,P^2)
  Y_p[i]=as.numeric(substr(files[n1+i],5,6))
}
diff=Y_p-XX_p%*%matrix(beta,P*P,1)
sum(abs(diff))/(n-n1)
setwd("D:\\desktop\\simulation new\\data set\\FGNET\\FGNET\\redo1")
save(beta,file="beta_tuc.Rdata")
save(diff,file="diff_tuc.Rdata")
