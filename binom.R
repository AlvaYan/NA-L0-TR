###########  binom
library(rTensor)
library(MASS)
library(glmnet)
library(e1071)
rm(list=ls())
set.seed(1)

n=300;p=4 ## n by p data matrix
beta0=matrix(rep((rnorm(p^3/8,0,1)),8),p^3,1)
stop=0.01

our=function(x,y,beta0,beta,betaols,data0)
{
  beta_0=fold(matrix(beta0,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_0_dec=tucker(beta_0,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)
  beta_0_eval=as.numeric(beta_0_dec$Z@data)
  
  beta_ols=fold(matrix(betaols,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_ols_dec=tucker(beta_ols,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)
  beta_ols_eval=as.numeric(beta_ols_dec$Z@data)
  
  #Choose model parameter
  ne=floor(p^3-2)
  m=1
  ne2=0
  multi=1
  lambda=50/multi
  iter=5000
  betamat=NULL
  betamat1=NULL
  eigvalmat=NULL
  eigvalmat1=NULL
  ll=matrix(rep(0),iter,1)
  ll_bar=matrix(rep(0),iter,1)
  i=1
  for(i in 1:iter){
    noise=NULL
    
    ## beta is stored in 2d matrix, therefore to use it as a tensor, we need to fold it first
    beta_1=fold(matrix(beta,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
    ## apply tucker decomposition to get the core tensor and unitary matrices
    beta_1_dec=tucker(beta_1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 20)
    ## treat the core tensor as eigen values and unitary matrices as eigenvectors
    eigval1=beta_1_dec$Z@data
    eigvec=beta_1_dec$U
    
    eigvalmat=rbind(eigvalmat,eigval1)
    
    if(i>600){
      sign=apply(sign( eigvalmat[(i-600):i,]),2,median)
      eigval2=matrix(apply( abs(eigvalmat[(i-600):i,]),2,median),p^3,m)
      sign2=sign(eigval2)
      eigval=eigval2*sign
    }else{
      sign=sign(beta_0_eval)
      sign2=sign(eigval1)
      eigval=(eigval1)
    }
    eigvalmat1=rbind(eigvalmat1,as.numeric(eigval))
    
    
    ## unfold the core tensor into a vector, make sure index matches with that used for unfloding beta 
    eigvec_expand=NULL
    for(ii in 1:p^3){
      eigvec_ii=NULL
      for(jj in 1:p){
        
        eigvec_ii=cbind(eigvec_ii,eigvec[[1]][index[ii,1],]%o%eigvec[[2]][index[ii,2],]*eigvec[[3]][index[ii,3],jj])
      }
      eigvec_expand=rbind(eigvec_expand,as.numeric(eigvec_ii))
    }

    ## set some numerical thresholds
    thred1=1e-6
    thred2=1e+6
    eigval=eigval*(abs(eigval)>thred1)+thred1*(abs(eigval)<=thred1) 
    
    ## generate noises
    for(j in 1:ne){
      linrad=rnorm(p^3,0,sqrt(lambda/abs(as.numeric(eigval))^4))
      linrad=linrad*(abs(linrad)<thred2)+thred2*(abs(linrad)>=thred2)*sign(linrad)
      lincomb=matrix((linrad),p^3,1)
      
      noise_j=matrix(as.numeric(t(matrix(eigvec_expand%*%lincomb,p,p^2))) ,1,p^3,byrow=TRUE)
      noise=rbind(noise,noise_j)
    }
    
    noisex=NULL
    noisey=NULL
    for(k in 1:multi){
      noisex=rbind(noisex,noise,-noise,noise,-noise)
      noisey=rbind(noisey,matrix(0,ne,m),matrix(0,ne,m),matrix(1,ne,m),matrix(1,ne,m) )
      
    }
    
    xe=rbind(x,noisex)
    ye=rbind(y,noisey)
    
    order=sample(1:nrow(xe),nrow(xe))
    xe=xe[order,]
    ye=ye[order]
    
    datae=data.frame(xe,ye)
    rg=glm(ye ~ xe-1, data=datae,family='binomial')
    beta1=rg[1]$coefficients
    
    betamat=rbind(betamat,as.numeric(beta1))
    ll[i]=logLik(rg)
    if(i>600){
      beta=matrix(apply(betamat[(i-600):i,],2,median),p^3,m)
      ll_bar[i]=median(ll[(i-600):i,])
    }else{
      if(i==iter){
        beta=matrix(apply(betamat[(i-600):i,],2,median),p^3,m)
        ll_bar[i]=median(ll[(i-600):i,])
      }else{
        beta=beta1
        ll_bar[i]=ll[i,]
      }
    }
    betamat1=rbind(betamat1,as.numeric(beta))
    if(!i%%50){print(i)}
    convergence=0
    if (i>4000)
    {
      if (abs(ll_bar[i]-ll_bar[i-1])<stop)
      #if (sum(abs(betamat1[i,]-betamat1[i-1,]))<0.01)
      {
        convergence=1
        break
      }
    }
  }
  return(list(beta,convergence,ll[i]))
}
tuc=function (U,G,x,y)#U is array 3*p*p,G is tensor p*p*p
{
  ite=5000
  ll=matrix(rep(0),ite,1)#llh function
  betamat=NULL
  i=1
  for (i in 1:ite)
  {
    #update each U(egvector), according to Tucker Tensor Regression and Neuroimaging Analysis
    for (j in 1:3)
    {
      xx=matrix(rep(0),n,p*p)#is the data that can be used as GLM reg to update Uj
      if(j==1)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[3,,],p,p),'mat2'=matrix(U[2,,],p,p))
        kp=kronecker_list(L)
        Gn=k_unfold(G,1)  
        Gn=(Gn@data)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),1)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp%*%t(Gn),1,p*p)
        }
        datanew=data.frame(xx,y)
        U[1,,]=glm(y ~ xx-1, data=datanew, family='binomial')[1]$coefficients
      }
      else if (j==2)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[3,,],p,p),'mat2'=matrix(U[1,,],p,p))
        kp=kronecker_list(L)
        Gn=k_unfold(G,2)  
        Gn=(Gn@data)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),2)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp%*%t(Gn),1,p*p)
        }
        datanew=data.frame(xx,y)
        U[2,,]=glm(y ~ xx-1,data=datanew,family='binomial')[1]$coefficients
      }
      else if (j==3)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[2,,],p,p),'mat2'=matrix(U[1,,],p,p))
        kp=kronecker_list(L)
        Gn=k_unfold(G,3)  
        Gn=(Gn@data)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),3)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp%*%t(Gn),1,p*p)
        }
        datanew=data.frame(xx,y)
        U[3,,]=glm(y ~ xx-1, data=datanew,family='binomial')[1]$coefficients
      }
    }
    #using tucker decom to optimize core tensor
    xx=matrix(rep(0),n,p*p*p)#is the data that can be used as GLM reg to update Uj
    L=list('mat1'=matrix(U[3,,],p,p),'mat2'=matrix(U[2,,],p,p),'mat3'=matrix(U[1,,],p,p))
    kp=kronecker_list(L)
    k=1
    for (k in 1:n)
    {
      Xn=vec(as.tensor(array(x[k,],dim=c(p,p,p))))
      xx[k,]=t(t(kp)%*%Xn)
    }
    datanew=data.frame(xx,y)
    reg=glm(y ~ xx-1, data=datanew,family='binomial')
    G=as.tensor(array(reg$coefficients,dim=c(p,p,p)))
    
    beta=kp%*%matrix(G@data,p^3,1)
    betamat=rbind(betamat, as.numeric(beta))
    print(i)
    ll[i,1]=logLik(reg)
    convergence=0
    if (i>1)
    {
      if (sum(abs(betamat[i,]-betamat[i-1,]))<stop)
      #if ((ll[i,1]-ll[i-1,1])<stop)
      {
        beta=kp%*%matrix(G@data,p^3,1)
        ll=ll[1:i]
        convergence=1
        break
      }
    }
  }
  return(list(beta,convergence,ll[i]))
}

tuc_reg=function (U,G,x,y)#U is array 3*p*p,G is tensor p*p*p
{
  ite=5000
  ll=matrix(rep(0),ite,1)#llh function
  betamat=NULL
  i=1
  jump=0
  for (i in 1:ite)
  {
    #update each U(egvector), according to Tucker Tensor Regression and Neuroimaging Analysis
    if (i>1){beta=kp1%*%matrix(G@data,p^3,1)}#preserve the beta from last ite
    j=1
    for (j in 1:3)
    {
      xx=matrix(rep(0),n,p*p)#is the data that can be used as GLM reg to update Uj
      if(j==1)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[3,,],p,p),'mat2'=matrix(U[2,,],p,p))
        kp=kronecker_list(L)
        Gn=k_unfold(G,1)  
        Gn=(Gn@data)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),1)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp%*%t(Gn),1,p*p)
        }
        datanew=data.frame(xx,y)
        temp=glm(as.numeric(y)~xx-1,family='binomial')
        Uold=matrix(U[1,,],p*p,1)
        U[1,,]=matrix(coef(temp),p,p)
        if (sum(is.na(U[1,,]))>0)
        {
          Unew=matrix(U[1,,],p*p,1)
          Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
          U[1,,]=matrix(Unew,p,p)
        }
      }
      else if (j==2)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[3,,],p,p),'mat2'=matrix(U[1,,],p,p))
        kp=kronecker_list(L)
        Gn=k_unfold(G,2)  
        Gn=(Gn@data)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),2)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp%*%t(Gn),1,p*p)
        }
        datanew=data.frame(xx,y)
        temp=glm(as.numeric(y)~xx-1,family='binomial')
        Uold=matrix(U[2,,],p*p,1)
        U[2,,]=matrix(coef(temp),p,p)
        if (sum(is.na(U[2,,]))>0)
        {
          Unew=matrix(U[2,,],p*p,1)
          Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
          U[2,,]=matrix(Unew,p,p)
        }
      }
      else if (j==3)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[2,,],p,p),'mat2'=matrix(U[1,,],p,p))
        kp=kronecker_list(L)
        Gn=k_unfold(G,3)  
        Gn=(Gn@data)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),3)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp%*%t(Gn),1,p*p)
        }
        datanew=data.frame(xx,y)
        temp=glm(as.numeric(y)~xx-1,family='binomial')
        Uold=matrix(U[3,,],p*p,1)
        U[3,,]=matrix(coef(temp),p,p)
        if (sum(is.na(U[3,,]))>0)
        {
          Unew=matrix(U[3,,],p*p,1)
          Unew[which(is.na(Unew)==TRUE)]=Uold[which(is.na(Unew)==TRUE)]
          U[3,,]=matrix(Unew,p,p)
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
    #using tucker decom to optimize core tensor
    xx=matrix(rep(0),n,p*p*p)#is the data that can be used as GLM reg to update Uj
    L=list('mat1'=matrix(U[3,,],p,p),'mat2'=matrix(U[2,,],p,p),'mat3'=matrix(U[1,,],p,p))
    kp=kronecker_list(L)
    kp1=kp
    k=1
    for (k in 1:n)
    {
      Xn=vec(as.tensor(array(x[k,],dim=c(p,p,p))))
      xx[k,]=t(t(kp)%*%Xn)
    }
    
    lambdas <- 10^seq(1, -5, by = -.1)
    lasso_reg <- cv.glmnet(xx, y, alpha = 1, lambda = lambdas, standardize = F, nfolds = 5, family='binomial',intercept=FALSE)
    lambda_best <- lasso_reg$lambda.min 
    datanew=data.frame(xx,y)
    lasso_model=glmnet(xx,y,alpha=1,lambda=lambda_best,standardize=FALSE,family='binomial',intercept=F)
    
    Gold=G
    G=as.tensor(array(lasso_model$beta,dim=c(p,p,p)))
    Gvec=matrix(G@data,1,p*p*p)
    Goldvec=matrix(Gold@data,1,p*p*p)
    Gvec[which(is.na(G@data)==TRUE)]=Goldvec[which(is.na(G@data)==TRUE)]
    G=as.tensor(array(Gvec,dim=c(p,p,p)))
    
    if (sum(is.na(G@data))>0)
    {
      convergence=0
      jump=1
      ll[i,1]=ll[i-1,1]
      break
    }
    print(i)
    beta=kp%*%matrix(G@data,p^3,1)
    betamat=rbind(betamat,as.numeric(beta))
    
    loglk=sum(dbinom(as.numeric(y),size=1,prob=exp(xx%*%G@data)/(1+exp(xx%*%G@data)),log=TRUE))
    ll[i,1]=loglk/length(y)+lambda_best*sum(abs(G@data))
    convergence=0
    if (i>1)
    {
      if (sum(abs(betamat[i,]-betamat[i-1,]))<stop)
      #if ((ll[i,1]-ll[i-1,1])<stop)
      {
        beta=kp%*%matrix(G@data,p^3,1)
        ll=ll[1:i,1]
        convergence=1
        break
      }
    }
  }
  return(list(beta,convergence,ll[i]))
}
cp_reg=function (U,x,y)#U is array 3*R*p
{
  ite=5000
  ll=matrix(rep(0),ite,1)#llh function
  betamat=NULL
  i=1
  for (i in 1:ite)
  {
    #update each U(egvector), according to Tucker Tensor Regression and Neuroimaging Analysis
    for (j in 1:3)
    {
      xx=matrix(rep(0),n,p*R)#is the data that can be used as GLM reg to update Uj
      if(j==1)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[3,,],p,R),'mat2'=matrix(U[2,,],p,R))
        kp=khatri_rao_list(L)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),1)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp,1,R*p)
        }
        lambdas <- 10^seq(-1, -5, by = -.1)
        lasso_reg <- cv.glmnet(xx, y, alpha = 0, lambda = lambdas, standardize = F, nfolds = 5, family='binomial',intercept=FALSE)
        lambda_best <- lasso_reg$lambda.min
        lasso_model=glmnet(xx,y,alpha=0,lambda=lambda_best,standardize=FALSE,family='binomial',intercept=F)
        U[1,,]=matrix(lasso_model$beta,p,R)
      }
      else if (j==2)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[3,,],p,R),'mat2'=matrix(U[1,,],p,R))
        kp=khatri_rao_list(L)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),2)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp,1,R*p)
        }
        lambdas <- 10^seq(-1, -5, by = -.1)
        lasso_reg <- cv.glmnet(xx, y, alpha = 0, lambda = lambdas, standardize = F, nfolds = 5, family='binomial',intercept=FALSE)
        lambda_best <- lasso_reg$lambda.min 
        lasso_model=glmnet(xx,y,alpha=0,lambda=lambda_best,standardize=FALSE,family='binomial',intercept=F)
        U[2,,]=matrix(lasso_model$beta,p,R)
      }
      else if (j==3)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[2,,],p,R),'mat2'=matrix(U[1,,],p,R))
        kp=khatri_rao_list(L)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),3)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp,1,R*p)
        }
        lambdas <- 10^seq(-1, -5, by = -.1)
        lasso_reg <- cv.glmnet(xx, y, alpha = 0, lambda = lambdas, standardize = F, nfolds = 5, family='binomial',intercept=FALSE)
        lambda_best <- lasso_reg$lambda.min 
        lasso_model=glmnet(xx,y,alpha=0,lambda=lambda_best,standardize=FALSE,family='binomial',intercept=F)
        U[3,,]=matrix(lasso_model$beta,p,R)
      }
    }
    #using tucker decom to optimize core tensor
    print(i)
    L=list('mat1'=matrix(U[3,,],p,R),'mat2'=matrix(U[2,,],p,R),'mat2'=matrix(U[1,,],p,R))
    kp=khatri_rao_list(L)
    beta=kp%*%matrix(rep(1),R,1)
    betamat=rbind(betamat,as.numeric(beta))
    e=xx%*%matrix(U[3,,],R*p,1)
    loglk=sum(dbinom(as.numeric(y),size=1,prob=exp(e)/(1+exp(e)),log=TRUE))
    ll[i,1]=loglk/length(y)+lambda_best*sum(matrix(U[3,,],R*p,1)^2)^0.5
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
  return(list(beta,convergence,ll[i,1]))
}
cp=function (U,x,y)#U is array 3*R*p
{
  ite=5000
  ll=matrix(rep(0),ite,1)#llh function
  betamat=NULL
  i=1
  i=i+1
  lambda_best=0.001
  for (i in 1:ite)
  {
    #update each U(egvector), according to Tucker Tensor Regression and Neuroimaging Analysis
    for (j in 1:3)
    {
      xx=matrix(rep(0),n,p*R)#is the data that can be used as GLM reg to update Uj
      if(j==1)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[3,,],p,R),'mat2'=matrix(U[2,,],p,R))
        kp=khatri_rao_list(L)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),1)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp,1,R*p)
        }
        datanew=data.frame(xx,y)
        U[1,,]=glm(y ~ xx-1, data=datanew, family='binomial')[1]$coefficients
      }
      else if (j==2)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[3,,],p,R),'mat2'=matrix(U[1,,],p,R))
        kp=khatri_rao_list(L)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),2)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp,1,R*p)
        }
        datanew=data.frame(xx,y)
        U[2,,]=glm(y ~ xx-1, data=datanew, family='binomial')[1]$coefficients
      }
      else if (j==3)
      {
        #using tucker decom to turn optimization of Uj into a GLM problem
        L=list('mat1'=matrix(U[2,,],p,R),'mat2'=matrix(U[1,,],p,R))
        kp=khatri_rao_list(L)
        for (k in 1:n)
        {
          Xn=k_unfold(as.tensor(array(x[k,],dim=c(p,p,p))),3)  
          Xn=(Xn@data)
          xx[k,]=matrix(Xn%*%kp,1,R*p)
        }
        datanew=data.frame(xx,y)
        reg=glm(y ~ xx-1, data=datanew, family='binomial')
        U[3,,]=reg[1]$coefficients
      }
    }
    #using tucker decom to optimize core tensor
    print(i)
    L=list('mat1'=matrix(U[3,,],p,R),'mat2'=matrix(U[2,,],p,R),'mat2'=matrix(U[1,,],p,R))
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
  return(list(beta,convergence,ll[i,1]))
}
rep=200
R=6

DV_bi=matrix(rep(0),rep,8)
LL_bi=matrix(rep(0),rep,7)
MSE_bi=matrix(rep(0),rep,8)
cMSE_bi=matrix(rep(0),rep,8)
CR_bi=matrix(rep(0),rep,9)
beta_bi_all=array(rep(0),dim=c(rep,7,p^3))
cvg_bi=matrix(rep(0),rep,5)
time=matrix(rep(0),rep,8)

r=1
sigma=0.25
for (r in 1:rep)
{
  #generate data
  x=matrix(rnorm(p^3*n,0,sigma),n,p^3)
  mu=x%*%beta0
  pr=1/(1+exp(-mu))#pr is already prob, why we hist it?should be hist(mu)?
  y=matrix((rbinom(n,1,as.numeric(pr))),n)
  data0=data.frame(x=x,y=y)
  
  #get the OLS estimator for beta witho~t regularization
  start=Sys.time()
  rg1=glm(y ~ x-1, data=data0,family='binomial')
  end=Sys.time()
  t6=end-start
  beta=betaols=rg1$coefficients
  start=Sys.time()
  rg2=svm(y~.,data=data0)
  end=Sys.time()
  t8=end-start
  
  start=Sys.time()
  lambdas <- 10^seq(1, -3, by = -.1)
  lasso_reg <- cv.glmnet(x, y, alpha = 1, lambda = lambdas, standardize = FALSE, nfolds = 5, family='binomial',intercept=FALSE)
  lambda_best <- lasso_reg$lambda.min
  rg3=glmnet(x,y,alpha=1,lambda=lambda_best,standardize=F,family='binomial',intercept=F)
  end=Sys.time()
  t7=end-start
  
  start=Sys.time()
  result1=our(x,y,beta0=beta0,beta=beta,betaols=betaols,data0)
  end=Sys.time()
  t1=end-start
  cvg_bi[r,1]=result1[[2]]
  
  #Li(2018) with regularization
  start=Sys.time()
  U=array(rnorm(3*p*p,0,1),dim=c(3,p,p))
  G=as.tensor(array(rnorm(p^3,0,1),dim=c(p,p,p)))
  result2=tuc(U,G,x,y)
  end=Sys.time()
  t2=end-start
  cvg_bi[r,2]=result2[[2]]
  start=Sys.time()
  result3=tuc_reg(U,G,x,y)
  end=Sys.time()
  t3=end-start
  cvg_bi[r,3]=result3[[2]]
  
  #Zhou(2013) with and without regularization
  start=Sys.time()
  U=array(rnorm(3*p*R,0,1),dim=c(3,p,R))
  result4=cp(U,x,y)
  end=Sys.time()
  t4=end-start
  cvg_bi[r,4]=result4[[2]]
  start=Sys.time()
  result5=cp_reg(U,x,y)
  end=Sys.time()
  t5=end-start
  cvg_bi[r,5]=result5[[2]]
  
  beta1=matrix(result1[[1]],p^3,1)
  beta2=matrix(result2[[1]],p^3,1)
  beta3=matrix(result3[[1]],p^3,1)
  beta4=matrix(result4[[1]],p^3,1)
  beta5=matrix(result5[[1]],p^3,1)
  beta6=matrix(betaols,p^3,1)
  beta7=matrix(rg3$beta,p^3,1)
  
  #record time
  time[r,1]=t1
  time[r,2]=t2
  time[r,3]=t3
  time[r,4]=t4
  time[r,5]=t5
  time[r,6]=t6
  time[r,7]=t7
  
  #store estimated beta
  beta_bi_all[r,1,]=matrix(beta1,1,p^3)
  beta_bi_all[r,2,]=matrix(beta2,1,p^3)
  beta_bi_all[r,3,]=matrix(beta3,1,p^3)
  beta_bi_all[r,4,]=matrix(beta4,1,p^3)
  beta_bi_all[r,5,]=matrix(beta5,1,p^3)
  beta_bi_all[r,6,]=matrix(beta6,1,p^3)
  beta_bi_all[r,7,]=matrix(beta7,1,p^3)
  betamle_bi_all[r,]=matrix(betaols,1,p^3)
  
  #sum of absolute error of parameter
  MSE_bi[r,1]=sum((beta1-beta0)^2)/p^3
  MSE_bi[r,2]=sum((beta2-beta0)^2)/p^3
  MSE_bi[r,3]=sum((beta3-beta0)^2)/p^3
  MSE_bi[r,4]=sum((beta4-beta0)^2)/p^3
  MSE_bi[r,5]=sum((beta5-beta0)^2)/p^3
  MSE_bi[r,6]=sum((beta6-beta0)^2)/p^3
  MSE_bi[r,7]=sum((beta7-beta0)^2)/p^3
  
  beta_1=fold(matrix(beta1,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_2=fold(matrix(beta2,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_3=fold(matrix(beta3,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_4=fold(matrix(beta4,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_5=fold(matrix(beta5,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_6=fold(matrix(beta6,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_7=fold(matrix(beta7,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_ols=fold(matrix(betaols,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_0=fold(matrix(as.numeric(beta0),p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  
  
  #Compute log-likelihood
  LL_bi[r,1]=result1[[3]]
  LL_bi[r,2]=result2[[3]]
  LL_bi[r,3]=result3[[3]]
  LL_bi[r,4]=result4[[3]]
  LL_bi[r,5]=result5[[3]]
  LL_bi[r,6]=logLik(rg1)#log lh for mle
  
  #Compute deviance
  DV_bi[r,7]=summary(rg1)$null.deviance
  DV_bi[r,6]=summary(rg1)$deviance
  DV_bi[r,1]=DV_bi[r,6]+2*LL_bi[r,6]-2*LL_bi[r,1]
  DV_bi[r,2]=DV_bi[r,6]+2*LL_bi[r,6]-2*LL_bi[r,2]
  DV_bi[r,3]=DV_bi[r,6]+2*LL_bi[r,6]-2*LL_bi[r,3]
  DV_bi[r,4]=DV_bi[r,6]+2*LL_bi[r,6]-2*LL_bi[r,4]
  DV_bi[r,5]=DV_bi[r,6]+2*LL_bi[r,6]-2*LL_bi[r,5]
  
  #Compute misclassification rate
  n1=200
  x_p=matrix(rnorm(p^3*n1,0,sigma),n1,p^3)
  mu=x_p%*%beta0
  pr=1/(1+exp(-mu))
  y_p=matrix((rbinom(n1,1,as.numeric(pr))),n1)
  mu1=x_p%*%beta1
  mu2=x_p%*%beta2
  mu3=x_p%*%beta3
  mu4=x_p%*%beta4
  mu5=x_p%*%beta5
  mu6=x_p%*%betaols
  mu7=x_p%*%beta7
  mu0=x_p%*%beta0
  pr1=1/(1+exp(-mu1))#pr is already prob, why we hist it?should be hist(mu)?
  pr2=1/(1+exp(-mu2))
  pr3=1/(1+exp(-mu3))
  pr4=1/(1+exp(-mu4))
  pr5=1/(1+exp(-mu5))
  pr6=1/(1+exp(-mu6))#pr is already prob, why we hist it?should be hist(mu)?
  pr7=1/(1+exp(-mu7))
  pr0=1/(1+exp(-mu0))
  y_1=as.numeric(pr1>0.5)
  y_2=as.numeric(pr2>0.5)
  y_3=as.numeric(pr3>0.5)
  y_4=as.numeric(pr4>0.5)
  y_5=as.numeric(pr5>0.5)
  y_6=as.numeric(pr6>0.5)
  y_7=as.numeric(pr7>0.5)
  y_8=predict(rg2,data.frame(x=x_p))
  y0=as.numeric(pr0>0.5)
  
  CR_bi[r,1]=sum((y_1-as.numeric(y_p))==0)/n1
  CR_bi[r,2]=sum((y_2-as.numeric(y_p))==0)/n1
  CR_bi[r,3]=sum((y_3-as.numeric(y_p))==0)/n1
  CR_bi[r,4]=sum((y_4-as.numeric(y_p))==0)/n1
  CR_bi[r,5]=sum((y_5-as.numeric(y_p))==0)/n1
  CR_bi[r,6]=sum((y_6-as.numeric(y_p))==0)/n1
  CR_bi[r,7]=sum((y_7-as.numeric(y_p))==0)/n1
  CR_bi[r,8]=sum((as.numeric(y_8>0.5))==(as.numeric(y_p)))/n1
  CR_bi[r,9]=sum((y0-as.numeric(y_p))==0)/n1
  
  print(r)
}

#save data after regression
save(DV_bi,file="DV_bi.Rdata")
save(LL_bi,file="LL_bi.Rdata")
save(MSE_bi,file="MSE_bi.Rdata")
save(cMSE_bi,file="cMSE_bi.Rdata")
save(CR_bi,file="CR_bi.Rdata")
save(beta_bi_all,file="beta_bi_all.Rdata")
save(cvg_bi,file="cvg_bi.Rdata")
save(time,file="time_bi.Rdata")

