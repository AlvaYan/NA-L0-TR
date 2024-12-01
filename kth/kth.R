library(rTensor)
library(MASS)
#library(extraDistr)
library(glmnet)
library(RNifti)
library(dplyr)
library(purrr)
library(abind)
rm(list=ls())

#---------------------first get all data------------
# Or you can skip this section and just load data in next section
library(reticulate)
use_python("C:\\Users\\14106\\AppData\\Local\\Programs\\Python\\Python37", required = TRUE)
np <- import("numpy")

aff.path='C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\data set\\kth\\preprocessed'

setwd(aff.path)
file_list <- list.files(path = aff.path, pattern = "*" )
set.seed(1)
file_list=sample(file_list, length(file_list), replace=F)
n=length(file_list)

P1=4#61
P2=4#73
P3=4#61

testing_person <- c('02', '03', '05', '06', '07', '08')
type2 <- c('boxing', 'handclapping', 'handwaving') # 1
X_train <- list()
X_test <- list()
Y_train <- vector()
Y_test <- vector()

XX=array(0,dim=c(n,P1,P2,P3))
# Iterate and categorize
file_path=file_list[1]
file_path="person01_jogging_d1_2.npz"
for (file_path in file_list) {
  filename <- basename(file_path)
  npz_data <- np$load(filename)
  keys <- npz_data$files
  data <- npz_data$get('arr_0')
  
  person <- substr(filename, 7, 8)
  activity_type <- substr(filename, 10, nchar(filename) - 9)
  
  if (person %in% testing_person) {
    X_test[[length(X_test) + 1]] <- data
    Y_test <- c(Y_test, activity_type %in% type2)
  } else {
    X_train[[length(X_train) + 1]] <- data
    Y_train <- c(Y_train, activity_type %in% type2)
  }
}

# Convert to arrays/matrices
n_train=length(X_train)
XX=c(X_train, X_test)
XX=abind(XX, along = 0)
Y=matrix(as.integer(c(Y_train, Y_test)),n,1)
rm(X_train, X_test, Y_train, Y_test)

Xn1=array(rep(0),c(n,P1,P2*P3))
Xn2=array(rep(0),c(n,P2,P1*P3))
Xn3=array(rep(0),c(n,P3,P1*P2))
Xn=array(rep(0),c(n,P1*P2*P3))
for (k in 1:n)
{
  Xn1[k,,]=k_unfold(as.tensor(XX[k,,,]),1)@data
  Xn2[k,,]=k_unfold(as.tensor(XX[k,,,]),2)@data
  Xn3[k,,]=k_unfold(as.tensor(XX[k,,,]),3)@data
  Xn[k,]=vec(as.tensor(XX[k,,,]))
}

setwd('C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\data set\\kth\\Rprocessed')
save(XX, file = "XX.RData")
save(Y, file = "Y.RData")
save(Xn, file = "Xn.RData")
save(Xn1, file = "Xn1.RData")
save(Xn2, file = "Xn2.RData")
save(Xn3, file = "Xn3.RData")
save(Xn3, file = "Xn3.RData")
save(n, file = "n.RData")
save(n_train, file = "n_train.RData")

#-------------------------------------
setwd('C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\data set\\kth\\Rprocessed')
load(file = "XX.RData")
load(file = "Y.RData")
load(file = "Xn.RData")
load(file = "Xn1.RData")
load(file = "Xn2.RData")
load(file = "Xn3.RData")
load(file = "n.RData")
load(file = "n_train.RData")

P1=4#61
P2=4#73
P3=4#61
aff.path='C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\data set\\kth\\preprocessed'

#---------Use simulated data----------------------------------
# P1=4#61
# P2=4#73
# P3=4#61
# XX=array(0,dim=c(n,P1,P2,P3))
# Y=matrix(rep(0),n,1)
# 
# sigma=0.25
# beta0=array(rep((rnorm(P1*P2*P3/8,0,1)),12),c(P1*P2*P3,1))
# Bd0=tucker(as.tensor(array(beta0,c(P1,P2,P3))),ranks=c(P1,P2,P3) ,tol=1e-6,max_iter = 50)
# (G0=Bd0$Z@data)
# Xn=array(rep(0),c(n,P1*P2*P3))
# set.seed(1)
# for (i in 1:n)
# {
#   XX[i,,,]=array(rnorm(P1*P2*P3,0,sigma),c(P1,P2,P3))
#   Xn[i,]=vec(as.tensor(XX[i,,,]))
#   #y=matrix(y_all_bi[r,,],n,1)#Y is n*1. Original tensor without decompositon
#   mu=Xn[i,]%*%beta0#+Z[i,]%*%R0
#   pr=1/(1+exp(-mu))#pr is already prob, why we hist it?should be hist(mu)?
#   #hist(pr)
#   Y[i]=rbinom(1,1,as.numeric(pr))
#   #data0=data.frame(x=x,y=y)
# }
# 
# Xn1=array(rep(0),c(n,P1,P2*P3))
# Xn2=array(rep(0),c(n,P2,P1*P3))
# Xn3=array(rep(0),c(n,P3,P1*P2))
# Xn=array(rep(0),c(n,P1*P2*P3))
# for (k in 1:n)
# {
#   Xn1[k,,]=k_unfold(as.tensor(XX[k,,,]),1)@data
#   Xn2[k,,]=k_unfold(as.tensor(XX[k,,,]),2)@data
#   Xn3[k,,]=k_unfold(as.tensor(XX[k,,,]),3)@data
#   Xn[k,]=vec(as.tensor(XX[k,,,]))
# }
# 
# n_train=floor(n/5*4)

############################################
#be careful of constant term
p1=4
p2=4
p3=4
# U1=array(rnorm(P1*p1,0,1),dim=c(P1,p1))
# U2=array(rnorm(P2*p2,0,1),dim=c(P2,p2))
# U3=array(rnorm(P3*p3,0,1),dim=c(P3,p3))
# G=array(rnorm(p1*p2*p3,0,1),dim=c(p1,p2,p3))
mle=glmnet(Xn[1:n_train,], Y[1:n_train],alpha=1,lambda=0,
           standardize=FALSE,family='binomial',intercept=F)
B=array(mle$beta[1:P1*P2*P3],dim=c(P1*P2*P3,1))
B[which(is.na(B)==TRUE)]=mean(B[which(is.na(B)==F)])
#R=array(mle$beta[(P1*P2*P3+1):length(mle$beta)],dim=c(length(feature)+1,1))
#R[which(is.na(R)==TRUE)]=mean(R[which(is.na(R)==F)])

#beta=array(rnorm(P1*P2*P3,0,1),dim=c(P1*P2*P3,1))

#Choose model parameter
ne=32
m=1
#ne2=0
multi=1
#lambda2=8
#sprs=F
lambda=50/multi
iter=3000
betamat=NULL
#betamat1=NULL
eigvalmat=NULL
#eigvalmat1=NULL
ll=matrix(rep(0),iter,1)
stop=0.1
#beta2=betaols

i=1#i+1
set.seed(1)
for(i in 1:iter){
  ## beta is stored in 2d matrix, therefore to use it as a tensor, we need to fold it first
  Bd=tucker(as.tensor(array(B,c(P1,P2,P3))),ranks=c(p1,p2,p3) ,tol=1e-6,max_iter = 50)
  #Bd=tucker(fold(matrix(B,P2,P1*P3),row_idx=2,col_idx=c(3,1),modes=c(P1,P2,P3)),ranks=c(p1,p2,p3) ,tol=1e-6,max_iter = 50)
  U1=Bd$U[[1]]
  U2=Bd$U[[2]]
  U3=Bd$U[[3]]
  G=Bd$Z@data
  
  #beta_1=as.tensor(matrix(beta,P,P))
  ### apply tucker decomposition to get the core tensor and unitary matrices
  #set.seed(1)
  #beta_1_dec=tucker(beta_1,ranks=c(p,p) ,tol=1e-6,max_iter = 20)
  ## treat the core tnsor as eigen values and unitary matrices as eigenvectors
  #eigval1=beta_1_dec$Z@data
  #eigvec=beta_1_dec$U
  
  eigvalmat=rbind(eigvalmat,G)
  
  if(i>600){
    sign=apply(sign( eigvalmat[(i-400):i,]),2,median)
    eigval_abs=matrix(apply( abs(eigvalmat[(i-400):i,]),2,median),p1*p2*p3,m)
    #sign2=sign(eigval2)
    eigval=eigval_abs*sign
  }else{
    #sign=sign(beta_0_eval)
    #sign2=sign(eigval1)
    eigval=(G)
  }
  #eigvalmat1=rbind(eigvalmat1,as.numeric(eigval))
  
  
  smalllizt <- list('mat1' = matrix(U3,P3,p3),
                    'mat2' = matrix(U2,P2,p2),
                    'mat3' = matrix(U1,P1,p1))
  eigvec_expand=kronecker_list(smalllizt)
  
  
  
  ## set some numerical thresholds
  thred1=1e-6
  thred2=1e+6
  eigval=eigval*(abs(eigval)>thred1)+thred1*(abs(eigval)<=thred1) 
  
  ## generate noises
  noise=NULL
  for(j in 1:ne){
    linrad=rnorm(p1*p2*p3,0,sqrt(lambda/abs(as.numeric(eigval))^4))
    linrad=linrad*(abs(linrad)<thred2)+thred2*(abs(linrad)>=thred2)*sign(linrad)
    lincomb=matrix((linrad),p1*p2*p3,1)
    
    noise_j=matrix(array(eigvec_expand%*%lincomb,dim=c(P1,P2,P3)),1,P1*P2*P3)
    #p=4
    #noise_j=matrix(as.numeric(t(matrix(eigvec_expand%*%lincomb,p,p^2))) ,1,p^3,byrow=TRUE)
    
    noise=rbind(noise,noise_j)
  }
  
  noisex=NULL
  noisey=NULL
  #noisez=NULL
  #zrep=100000
  for(k in 1:multi){
    noisex=rbind(noisex,noise,-noise,noise,-noise)
    noisey=rbind(noisey,matrix(0,ne,m),matrix(0,ne,m),matrix(1,ne,m),matrix(1,ne,m) )
    #noisez=rbind(noisez,matrix(zrep,ne,length(R)),matrix(zrep,ne,length(R)),matrix(-zrep,ne,length(R)),matrix(-zrep,ne,length(R)))
  }
  
  xe=rbind(Xn[1:n_train,],noisex)
  ye=rbind(matrix(Y[1:n_train,],n_train,1),noisey)
  #ze=rbind(Z[1:n_train,],noisez)#matrix(0,length(noisey),length(R)))
  
  order=sample(1:nrow(xe),nrow(xe))
  xe=xe[order,]
  ye=ye[order]
  #ze=ze[order,]
  
  datae=data.frame(xe,ye)
  
  ## do linear regression on vectorized tensor data
  # reg=glmnet(xe,ye,alpha=0,lambda=0,
  #            standardize=FALSE,family='binomial',intercept=F)
  #            #offset = ze%*%R)
  #            #offset = rbind(Z[1:n_train]%*%(R),matrix(1,length(noisey),1)))
  # Bnew=array(reg$beta,dim=c(P1*P2*P3,1))
  
  reg=glm(ye ~ xe - 1, family = binomial())
  Bold=B
  Bnew=array(reg$coefficients,dim=c(P1*P2*P3,1))
  Bnew[which(is.na(Bnew)==TRUE)]=Bold[which(is.na(Bnew)==TRUE)]
  #Rnew[which(Rnew==0)]=Rold[which(Rnew==0)]
  B=Bnew
  ## record the beta estimation
  betamat=rbind(betamat,as.numeric(B))
  #Bold=B
  ## moving average
  if(i>600){
    B=matrix(apply(betamat[(i-400):i,],2,median),P1*P2*P3,m)
  }else{
    if(i==iter){
      B=matrix(apply(betamat[(i-400):i,],2,median),P1*P2*P3,m)
    }
    else{
      B=B
    }
  }
  
  #-----------------------
  if(!i%%400){print(i)}
  converge=0
  if (i>1)
  {
    #if (abs(mean(ll[(i-400):i])-mean(ll[(i-401):(i-1)]))<100)
    if (sum(abs(Bold-B))<0.01)
    {
      converge=1
      break
    }
  }
}


#Doing prediction
lg=Xn[(n_train+1):n,]%*%B#+Z[(n_train+1):n,]%*%R
#B_mle=o[1:(P1*P2*P3)]
#R_mle=o[(P1*P2*P3+1):(length(o)-1)]
p=1/(1+exp(-lg))
sum((p>0.5)==Y[(n_train+1):n])/(n-n_train)

set.seed(1)
glm <- glm(Y[1:n_train] ~ Xn[1:n_train,] - 1, family = binomial())
o=coef((glm))
summary(glm)
lg=Xn[(n_train+1):n,]%*%o#+Z[(n_train+1):n,]%*%R
#B_mle=o[1:(P1*P2*P3)]
#R_mle=o[(P1*P2*P3+1):(length(o)-1)]
p=1/(1+exp(-lg))
sum((p>0.5)==Y[(n_train+1):n])/(n-n_train)
#sum((pred>0.5)==Y[(n_train+1):n,])/(n-n_train)

set.seed(123)
cvglm <- cv.glmnet(Xn[1:n_train,], Y[1:n_train],
                 family = "binomial", alpha = 1,intercept=F)#, lambda = 0.1)
glmcv <- glmnet(Xn[1:n_train,], Y[1:n_train], 
              lambda=cvglm$lambda.min, family = "binomial", intercept=F,
              alpha = 1)#, lambda = 0.1)
oo=as.numeric(glmcv$beta)
pred=predict(glmcv,newx=Xn[(n_train+1):n,], type=c("response"))
sum((pred>0.5)==Y[(n_train+1):n,])/(n-n_train)

setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\data set\\kth")
save(B,file="beta_our.Rdata")
save(o,file="beta_glm.Rdata")
save(oo,file="beta_cvglm.Rdata")
#save(p,file="p_our.Rdata")

