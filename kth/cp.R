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
R=6
set.seed(1)
# U1=array(rnorm(P1*p1,0,1),dim=c(P1,p1))
# U2=array(rnorm(P2*p2,0,1),dim=c(P2,p2))
# U3=array(rnorm(P3*p3,0,1),dim=c(P3,p3))
# G=array(rnorm(p1*p2*p3,0,1),dim=c(p1,p2,p3))
B=array(rnorm(P1*P2*P3,0,1),dim=c(P1*P2*P3,1))
Bd=cp(as.tensor(array(B,c(P1,P2,P3))),num_components = R,  max_iter = 50, tol = 1e-05)
U1=Bd$U[[1]]
U2=Bd$U[[2]]
U3=Bd$U[[3]]

stop=0.01
ite=1000
jump=0
lam=10^(-10)
ll=matrix(rep(0),ite,1)#loss function
i=1
i=i+1
for (i in 1:ite)
{
    xx=matrix(rep(0),n_train,P1*R)#is the data that can be used as GLM reg to update Uj
    L=list('mat1'=matrix(U3,P3,R),'mat2'=matrix(U2,P2,R))
    kp=khatri_rao_list(L)
    for (k in 1:n_train)
    {
      #Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),1)  
      #Xn=(Xn@data)
      xx[k,]=matrix(Xn1[k,,]%*%kp,1,R*P1)
    }
    
    mod=glm(Y[1:n_train] ~ xx-1, family='binomial')
    Unew=array(mod[1]$coefficients,dim=c(P1,R))
    
    # mod=glmnet(xx,Y[1:n_train],alpha=1,lambda=0,
    #            standardize=FALSE,family='binomial',intercept=F,
    #            offset = Z[1:n_train,]%*%(RR))
    # Unew=array(mod$beta,dim=c(P1,R))
    Uold=U1
    Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
    #Unew[which(Unew==0)]=Uold[which(Unew==0)]
    U1=Unew
    
    #---------------
    #using tucker decom to turn optimization of Uj into a GLM problem
    xx=matrix(rep(0),n_train,P2*R)#is the data that can be used as GLM reg to update Uj
    L=list('mat1'=matrix(U3,P3,R),'mat2'=matrix(U1,P1,R))
    kp=khatri_rao_list(L)
    for (k in 1:n_train)
    {
      #Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),2)  
      #Xn=(Xn@data)
      xx[k,]=matrix(Xn2[k,,]%*%kp,1,R*P2)
    }
    mod=glm(Y[1:n_train] ~ xx-1, family='binomial')
    Unew=array(mod[1]$coefficients,dim=c(P2,R))
    
    # mod=glmnet(xx,Y[1:n_train],alpha=1,lambda=0,
    #            standardize=FALSE,family='binomial',intercept=F,
    #            offset = Z[1:n_train,]%*%(RR))
    # Unew=array(mod$beta,dim=c(P2,R))
    Uold=U2
    Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
    #Unew[which(Unew==0)]=Uold[which(Unew==0)]
    U2=Unew
    #---------------
    #using tucker decom to turn optimization of Uj into a GLM problem
    xx=matrix(rep(0),n_train,P3*R)#is the data that can be used as GLM reg to update Uj
    L=list('mat1'=matrix(U2,P2,R),'mat2'=matrix(U1,P1,R))
    kp=khatri_rao_list(L)
    for (k in 1:n_train)
    {
      #Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),3)  
      #Xn=(Xn@data)
      xx[k,]=matrix(Xn3[k,,]%*%kp,1,R*P3)
    }
    
    mod=glm(Y[1:n_train] ~ xx-1, family='binomial')
    Unew=array(mod[1]$coefficients,dim=c(P3,R))
    
    # mod=glmnet(xx,Y[1:n_train],alpha=1,lambda=0,
    #            standardize=FALSE,family='binomial',intercept=F,
    #            offset = Z[1:n_train,]%*%(RR))
    # Unew=array(mod$beta,dim=c(P3,R))
    Uold=U3
    Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
    #Unew[which(Unew==0)]=Uold[which(Unew==0)]
    U3=Unew
    
    #---------------
    L=list('mat1'=matrix(U3,P3,R),'mat2'=matrix(U2,P2,R),'mat3'=matrix(U1,P1,R))
    kp=khatri_rao_list(L)
    Bold=B
    B=kp%*%matrix(rep(1),R,1)
    
    #-----------------------
    print(i)
    # L=list('mat1'=matrix(U3,P3,R),'mat2'=matrix(U2,P2,R),'mat3'=matrix(U1,P1,R))
    # kp=khatri_rao_list(L)
    # #loss[i,]=sum((x%*%beta0-xx%*%matrix(G@data,p*p*p,1))^2)+sum(abs(G@data))
    # #e=xx%*%matrix(U[3,,],R*p,1)
    # #loglk=sum(dbinom(as.numeric(y),size=1,prob=exp(e)/(1+exp(e)),log=TRUE))
    # #ll[i,1]=loglk/length(y)+lambda_best*sum(matrix(U[3,,],R*p,1)^2)^0.5
    # #ll[i,1]=logLik(reg)
    # #L=list('mat1'=matrix(U[3,,],p,R),'mat2'=matrix(U[2,,],p,R),'mat2'=matrix(U[1,,],p,R))
    # #kp=khatri_rao_list(L)
    # Bold=B
    # B=kp%*%matrix(rep(1),R,1)
    convergence=0
    if (i>1)
    {
      #if (abs(ll[i,1]-ll[i-1,1])<stop)
      if (sum(abs(Bold-B))/length(B)<stop)
      #if (sum(abs(Bold-B))/length(B)<stop)
      {
        #beta=kp%*%matrix(G@data,p1*p2*p3,1)
        convergence=1
        break
      }
    }
}

#Doing prediction
lg=Xn[(n_train+1):n,]%*%B#+Z[(n_train+1):n,]%*%RR
p=1/(1+exp(-lg))
sum((p>0.5)==Y[(n_train+1):n])/(n-n_train)

setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\data set\\kth")
save(B,file="beta_cp.Rdata")