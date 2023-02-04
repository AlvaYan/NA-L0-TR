library(rTensor)
library(MASS)
library(glmnet)
rm(list=ls())
set.seed(2)

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

###############
lambdas <- 10^seq(1, -4, by = -.1)
lasso_reg <- cv.glmnet(XX, Y, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
lambda_best <- lasso_reg$lambda.min 
rg1=glmnet(XX,Y,alpha=1,lambda=lambda_best,standardize=FALSE,family='gaussian')

#Doing prediction
setwd(aff.path)
XX_p=matrix(rep(0),n-n1,P*P)
Y_p=matrix(rep(0),n-n1,1)
for (i in 1:(n-n1))
{
  #X[i,,]=t(as.matrix(read.csv(files[i],header = F)))
  XX_p[i,]=matrix(t(as.matrix(read.csv(files[n1+i],header = F))),1,P^2)
  name=strsplit(files[n1+i], split=".",fixed=T)[[1]][1]
  sp=unlist(strsplit(name,split="(?=[+-])",perl=T))
  Y_p[i]=as.numeric(paste(sp[4],sp[5],collapse="",sep=""))
}
diff=Y_p-predict(rg1,newx=XX_p)
sum(abs(diff))/(n-n1)

save(beta,file="beta_olsreg.Rdata")
save(diff,file="diff_olsreg.Rdata")
