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
#p=29
R=25

#free=p*P+p*P+p*p-p^2-p^2
#free-R*(P+P)-R^2

n1=floor(n*9/10)#length(files)
#X=array(rep(0),dim=c(n,P,P))
XX=matrix(rep(0),n1,P*P)
Y=matrix(rep(0),n1,1)
#YY=matrix(rep(0),n1,1)
i=1
for (i in 1:n1)
{
  #X[i,,]=t(as.matrix(read.csv(files[i],header = F)))
  XX[i,]=matrix(t(as.matrix(read.csv(files[i],header = F))),1,P^2)
  Y[i]=as.numeric(substr(files[i],5,6))
}
###############
#n=200;p=4 ## n by p matrix , x???? n by p^3 matrix
#P=10
#R=3
#beta0=matrix(rep(rnorm(P^2/2,0,1),2),P^2,1)
#XX=matrix(rnorm(P^2*n,0,1),n,P^2)
#Y=XX%*%beta0+matrix(rnorm(n,0,0.5),n,1)#Y is n*1. Original tensor without decompositon
########################
U=array(rnorm(P*R*3,0,1),dim=c(3,P,R))
stop=0.01
ite=4000
ll=matrix(rep(0),ite,1)#llh function
betamat=NULL
i=1
i=i+1
jump=0
lam=10^(-10)
for (i in 1:ite)
{
  #update each U(egvector), according to Tucker Tensor Regression and Neuroimaging Analysis
  #j=2
  for (j in 1:2)
  {
    xx=matrix(rep(0),n1,P*R)#is the data that can be used as GLM reg to update Uj
    if(j==1)
    {
      #using tucker decom to turn optimization of Uj into a GLM problem
      #L=list('mat1'=matrix(U[3,,],p,R),'mat2'=matrix(U[2,,],p,R))
      #kp=khatri_rao_list(L)
      kp=matrix(U[2,,],P,R)
      for (k in 1:n1)
      {
        Xn=k_unfold(as.tensor(array(XX[k,],dim=c(P,P))),1)  
        Xn=(Xn@data)
        xx[k,]=matrix(Xn%*%kp,1,R*P)
      }
      datanew=data.frame(xx,Y)
      U[1,,]=lm(Y ~ xx-1, data=datanew)[1]$coefficients
      #Gn=t(Gn)
    }else if (j==2){
      
      #using tucker decom to turn optimization of Uj into a GLM problem
      #L=list('mat1'=matrix(U[3,,],p,R),'mat2'=matrix(U[1,,],p,R))
      #kp=khatri_rao_list(L)
      kp=matrix(U[1,,],P,R)
      for (k in 1:n1)
      {
        Xn=k_unfold(as.tensor(array(XX[k,],dim=c(P,P))),2)  
        Xn=(Xn@data)
        xx[k,]=matrix(Xn%*%kp,1,R*P)
      }
      datanew=data.frame(xx,Y)
      reg=lm(Y ~ xx-1, data=datanew)
      U[2,,]=reg[1]$coefficients
    }
  }
  #using tucker decom to optimize core tensor
  
  print(i)
  #loss[i,]=sum((x%*%beta0-xx%*%matrix(G@data,p*p*p,1))^2)+sum(abs(G@data))
  L=list('mat1'=matrix(U[2,,],P,R),'mat2'=matrix(U[1,,],P,R))
  kp=khatri_rao_list(L)
  beta=kp%*%matrix(rep(1),R,1)
  betamat=rbind(betamat,as.numeric(beta))
  ll[i,1]=logLik(reg)
  convergence=0
  if (i>1)
  {
    if (sum(abs(betamat[i,]-betamat[i-1,]))<stop)
    #if ((ll[i,1]-ll[i-1,1])<stop)
    {
      beta=kp%*%matrix(rep(1),R,1)
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
setwd("D:\\desktop\\simulation new\\data set\\FGNET\\FGNET\\redo2")
save(beta,file="beta_cp.Rdata")
save(diff,file="diff_cp.Rdata")
#17.4337, R=20