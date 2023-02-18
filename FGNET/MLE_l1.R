library(rTensor)
library(MASS)
#library(extraDistr)
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
n1=floor(n*9/10)#length(files)
#X=array(rep(0),dim=c(n,P,P))
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
###############
#n=200;p=4 ## n by p matrix , x???? n by p^3 matrix
#P=10
#beta0=matrix(rep(rnorm(P^2/2,0,1),2),P^2,1)
#XX=matrix(rnorm(P^2*n,0,1),n,P^2)
#Y=XX%*%beta0+matrix(rnorm(n,0,0.5),n,1)#Y is n*1. Original tensor without decompositon
########################
lambdas <- 10^seq(1, -4, by = -.1)
lasso_reg <- cv.glmnet(XX, Y, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
lambda_best <- lasso_reg$lambda.min 
rg1=glmnet(XX,Y,alpha=1,lambda=lambda_best,standardize=FALSE,family='gaussian')
beta=matrix(rg1$beta)
#return(list(beta,convergence,ll[i]))
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
diff=Y_p-predict(rg1,newx=XX_p)
sum(abs(diff))/(n-n1)
setwd("D:\\desktop\\simulation new\\data set\\FGNET\\FGNET\\redo1")
save(beta,file="beta_l1MLE.Rdata")
save(diff,file="diff_l1MLE.Rdata")

# > sum(abs(diff))/(n-n1)
# [1] 4.297829