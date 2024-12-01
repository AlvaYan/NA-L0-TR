library(rTensor)
library(MASS)
#library(extraDistr)
library(glmnet)
library(RNifti)
library(dplyr)
library(purrr)
library(abind)
rm(list=ls())

#-------------------------------------
setwd('C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\data set\\kth\\Rprocessed')
load(file = "XX.RData")
load(file = "Y.RData")
load(file = "Xn.RData")
load(file = "Xn1.RData")
load(file = "Xn2.RData")
load(file = "Xn3.RData")
load(file = "Xn3.RData")
load(file = "n.RData")
load(file = "n_train.RData")

P1=4#61
P2=4#73
P3=4#61
aff.path='C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\data set\\kth\\preprocessed'

#-------------------------------------------
#be careful of constant term
p1=4
p2=4
p3=4

B=array(rnorm(P1*P2*P3,0,1),dim=c(P1*P2*P3,1))
Bd=tucker(as.tensor(array(B,c(P1,P2,P3))),ranks=c(p1,p2,p3) ,tol=1e-6,max_iter = 50)
U1=Bd$U[[1]]
U2=Bd$U[[2]]
U3=Bd$U[[3]]
G=Bd$Z@data

stop=0.01
ite=1000
jump=0
lam=10^(-10)
ll=matrix(rep(0),ite,1)#loss function
i=1
set.seed(1)
for (i in 1:ite)
{
  #update each U(egvector), according to Tucker Tensor Regression and Neuroimaging Analysis
  
  #-------------
  #using tucker decom to turn optimization of Uj into a GLM problem
  xx=matrix(rep(0),n_train,P1*p1)#is the data that can be used as GLM reg to update Uj
  L=list('mat1'=matrix(U3,P3,p3),'mat2'=matrix(U2,P2,p2))
  kp=kronecker_list(L)
  #kp=matrix(U[2,,],P,p)
  Gn=k_unfold(as.tensor(G),1)  
  Gn=(Gn@data)
  k=1
  for (k in 1:n_train)
  {
    xx[k,]=matrix(Xn1[k,,]%*%kp%*%t(Gn),1,P1*p1)
  }
  datanew=data.frame(xx,Y[1:n_train])
  mod=glm(Y.1.n_train. ~ xx-1, data=datanew,family='binomial')#, offset = Z[1:n_train,]%*%R)
  Unew=array(mod[1]$coefficients,dim=c(P1,p1))
  
  # mod=glmnet(xx,Y[1:n_train],alpha=1,lambda=0,
  #            standardize=FALSE,family='binomial',intercept=F)
  #            #offset = Z[1:n_train,]%*%(R))
  # Unew=array(mod$beta,dim=c(P1,p1))
  
  Uold=U1
  Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
  #Unew[which(Unew==0)]=Uold[which(Unew==0)]
  U1=Unew
  
  #logLik(mod)
  #Gn=t(Gn)
  
  #-------------
  #using tucker decom to turn optimization of Uj into a GLM problem
  xx=matrix(rep(0),n_train,P2*p2)#is the data that can be used as GLM reg to update Uj
  L=list('mat1'=matrix(U3,P3,p3),'mat2'=matrix(U1,P1,p1))
  kp=kronecker_list(L)
  #kp=matrix(U[2,,],P,p)
  Gn=k_unfold(as.tensor(G),2)  
  Gn=(Gn@data)
  for (k in 1:n_train)
  {
    xx[k,]=matrix(Xn2[k,,]%*%kp%*%t(Gn),1,P2*p2)
  }
  datanew=data.frame(xx,Y[1:n_train])
  mod=glm(Y.1.n_train. ~ xx-1, data=datanew,family='binomial')#, offset = Z[1:n_train,]%*%(R))
  Unew=array(mod[1]$coefficients,dim=c(P2,p2))
  
  # mod=glmnet(xx,Y[1:n_train],alpha=1,lambda=0,
  #            standardize=FALSE,family='binomial',intercept=F)
  #            #offset = Z[1:n_train,]%*%(R))
  # Unew=array(mod$beta,dim=c(P2,p2))
  
  Uold=U2
  Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
  #Unew[which(Unew==0)]=Uold[which(Unew==0)]
  U2=Unew
  
  #-------------
  #using tucker decom to turn optimization of Uj into a GLM problem
  xx=matrix(rep(0),n_train,P3*p3)#is the data that can be used as GLM reg to update Uj
  L=list('mat1'=matrix(U2,P2,p2),'mat2'=matrix(U1,P1,p1))
  kp=kronecker_list(L)
  Gn=k_unfold(as.tensor(G),3)
  Gn=(Gn@data)
  for (k in 1:n_train)
  {
    xx[k,]=matrix(Xn3[k,,]%*%kp%*%t(Gn),1,P3*p3)
  }
  datanew=data.frame(xx,Y[1:n_train])
  mod=glm(Y.1.n_train. ~ xx-1, data=datanew,family='binomial')#, offset = Z[1:n_train,]%*%(R))
  Unew=array(mod[1]$coefficients,dim=c(P3,p3))
  
  # mod=glmnet(xx,Y[1:n_train],alpha=1,lambda=0,
  #            standardize=FALSE,family='binomial',intercept=F)
  #            #offset = Z[1:n_train,]%*%(R))
  # Unew=array(mod$beta,dim=c(P3,p3))
  
  Uold=U3
  Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
  #Unew[which(Unew==0)]=Uold[which(Unew==0)]
  U3=Unew
  #logLik(mod)
  
  #---------------------------
  #using tucker decom to optimize core tensor
  xx=matrix(rep(0),n_train,p1*p2*p3)#is the data that can be used as GLM reg to update Uj
  L=list('mat1'=matrix(U3,P3,p3),'mat2'=matrix(U2,P2,p2),'mat3'=matrix(U1,P1,p1))
  kp=kronecker_list(L)
  k=1
  for (k in 1:n_train)
  {
    xx[k,]=t(t(kp)%*%matrix(Xn[k,],P1*P2*P3,1))
  }
  datanew=data.frame(xx,Y[1:n_train])
  reg=glm(Y.1.n_train. ~ xx-1, data=datanew,family='binomial')#, offset = Z[1:n_train,]%*%(R))
  Gnew=array(reg$coefficients,dim=c(p1,p2,p3))
  
  # reg=glmnet(xx,Y[1:n_train],alpha=1,lambda=0,
  #            standardize=FALSE,family='binomial',intercept=F)
  #            #offset = Z[1:n_train,]%*%(R))
  # Gnew=array(reg$beta,dim=c(p1,p2,p3))
  Gold=G
  Gnew[which(is.na(Gnew)==TRUE)]=Gold[which(is.na(Gnew)==TRUE)]
  #Gnew[which(Gnew==0)]=Gold[which(Gnew==0)]
  G=Gnew

  #-----------------
  Bold=B
  B=kp%*%matrix(G,p1*p2*p3,1)
  #---------------------
  print(i)
  #ll[i,]=logit_logLik(Y[1:n_train],predict(reg, xx, type="response"))
  convergence=0
  if (i>1)
  {
    #if (abs(ll[i,1]-ll[i-1,1])<stop)
    if (sum(abs(Bold-B))/length(B)<stop)
    {
      #beta=kp%*%matrix(G@data,p1*p2*p3,1)
      convergence=1
      break
    }
  }
}
#-------------------
#Doing prediction
lg=Xn[(n_train+1):n,]%*%B#+Z[(n_train+1):n,]%*%R
p=1/(1+exp(-lg))
sum((p>0.5)==Y[(n_train+1):n])/(n-n_train)

setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\data set\\kth")
save(B,file="beta_tuc.Rdata")
#save(p,file="p_tuc.Rdata")
