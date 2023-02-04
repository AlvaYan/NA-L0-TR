library(rTensor)
library(MASS)
library(glmnet)
rm(list=ls())
set.seed(2)

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
######################
U=array(rnorm(P*p*2,0,1),dim=c(2,P,p))
G=as.tensor(array(rnorm(p^2,0,1),dim=c(p,p)))
stop=0.01
ite=5000
ll=matrix(rep(0),ite,1)
betamat=NULL
i=1
i=i+1
jump=0
lam=10^(-10)
for (i in 1:ite)
{
  if (i>1){beta=kp1%*%matrix(G@data,p^2,1)}#preserve the beta from last ite
  j=1
  for (j in 1:2)
  {
    xx=matrix(rep(0),n1,P*p)#is the data that can be used as GLM reg to update Uj
    if(j==1)
    {
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
      temp=glm(as.numeric(Y)~xx-1,family='gaussian')
      Uold=matrix(U[1,,],P*p,1)
      U[1,,]=matrix(coef(temp),P,p)
      if (sum(is.na(U[1,,]))>0)
      {
        Unew=matrix(U[1,,],P*p,1)
        Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
        U[1,,]=matrix(Unew,P,p)
      }
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
      temp=glm(as.numeric(Y)~xx-1,family='gaussian')
      Uold=matrix(U[2,,],P*p,1)
      U[2,,]=matrix(coef(temp),P,p)
      if (sum(is.na(U[2,,]))>0)
      {
        Unew=matrix(U[2,,],P*p,1)
        Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
        U[2,,]=matrix(Unew,P,p)
      }
    }
    if (jump==1)
    {
      break
    }
  }
  if (jump==1)
  {
    break
  }
  xx=matrix(rep(0),n1,p*p)#is the data that can be used as GLM reg to update Uj
  L=list('mat1'=matrix(U[2,,],P,p),'mat2'=matrix(U[1,,],P,p))
  kp=kronecker_list(L)
  kp1=kp
  
  k=1
  for (k in 1:n1)
  {
    Xn=vec(as.tensor(array(XX[k,],dim=c(P,P))))
    xx[k,]=t(t(kp)%*%Xn)
  }
  lambdas <- 10^seq(-1, -5, by = -.1)
  lasso_reg <- cv.glmnet(xx, Y, alpha = 1, standardize = FALSE, nfolds = 5, intercept=FALSE)
  
  # Best 
  lambda_best <- lasso_reg$lambda.min 
  lasso_model=glmnet(xx,Y,alpha=1,lambda=lambda_best,standardize=FALSE,family='gaussian',intercept=F)
  Gold=G
  G=as.tensor(array(lasso_model$beta,dim=c(p,p)))
  Gvec=matrix(G@data,1,p*p)
  Goldvec=matrix(Gold@data,1,p*p)
  Gvec[which(is.na(G@data)==TRUE)]=Goldvec[which(is.na(G@data)==TRUE)]
  G=as.tensor(array(Gvec,dim=c(p,p)))
  if (sum(is.na(G@data))>0)
  {
    convergence=0
    jump=1
    ll[i,1]=ll[i-1,1]
    break
  }
  beta=kp%*%matrix(G@data,p^2,1)
  betamat=rbind(betamat,as.numeric(beta))
  rss=sum((Y-xx%*%matrix(G@data,400,1))^2/2)
  ll[i,1]=rss/length(Y)+lambda_best*sum(abs(G@data))
  print(i)
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
  Y_p[i]=as.numeric(substr(files[n1+i],5,6))
}
diff=Y_p-XX_p%*%matrix(beta,P*P,1)
sum(abs(diff))/(n-n1)
setwd("D:\\desktop\\simulation new\\data set\\FGNET\\FGNET\\redo1")
save(beta,file="beta_tr.Rdata")
save(diff,file="diff_tr.Rdata")
