#########################################################  linear
library(rTensor)
library(MASS)
library(glmnet)
rm(list=ls())
set.seed(1)

n=300;p=4 ## n by p matrix 
setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\linear\\new6\\simulated data\\")
beta0=matrix(rep((rnorm(p^3/2,0,1)),2),p^3,1)
beta_0=as.tensor(array(beta0,dim=c(p,p,p)))
tucker(as.tensor(array(beta0,dim=c(p,p,p))),ranks=c(p,p,p),tol=1e-6,max_iter = 50)$Z@data

stop=0.01

our=function(x,y,data0, n.ne)
{
  beta=betaols=glm(y ~ x-1, data=data0,family='gaussian')[1]$coefficients
  beta_ols=fold(matrix(betaols,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
  beta_ols_dec=tucker(beta_ols,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)
  beta_ols_eval=as.numeric(beta_ols_dec$Z@data)
  
  index=matrix(0,p^3,3)
  for( i in 1:p^3){
    index[,1]= rep(1:p,p^2)
    index[,2]= rep(rep(1:p,rep(p,p)),p)
    index[,3]= rep(1:p,rep(p^2,p))
  }
  
  #Choose model parameter
  ne=n.ne # floor(p^3-2)
  m=1
  ne2=0
  multi=1
  lambda2=8
  sprs=F
  lambda=10**-2/multi # 
  iter=30000
  betamat=NULL
  betamat1=NULL
  eigvalmat=NULL
  eigvalmat1=NULL
  ll=matrix(rep(0),iter,1)
  ll_bar=matrix(rep(0),iter,1)
  beta2=betaols
  
  for(i in 1:iter){
    noise=NULL
    
    ## beta is stored in 2d matrix, therefore to use it as a tensor, we need to fold it first
    beta_1=fold(matrix(beta,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
    beta_1_dec=tucker(beta_1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 20)
    eigval1=beta_1_dec$Z@data
    eigvec=beta_1_dec$U
    
    eigvalmat=rbind(eigvalmat,eigval1)
    
    if(i>600){
      sign=apply(sign( eigvalmat[(i-600):i,]),2,median)
      eigval2=matrix(apply( abs(eigvalmat[(i-600):i,]),2,median),p^3,m)
      eigval=eigval2*sign
    }
    else{
      eigval=(eigval1)
    }
    eigvalmat1=rbind(eigvalmat1,as.numeric(eigval))
    
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
      noisex=rbind(noisex,noise)
      noisey=rbind(noisey,matrix(0,ne,m))
    }
    
    xe=rbind(x,noisex)
    ye=rbind(y,noisey)
    
    order=sample(1:nrow(xe),nrow(xe))
    xe=xe[order,]
    ye=ye[order]
    
    datae=data.frame(xe,ye)
    
    ## do linear regression on vectorized tensor data
    rg=lm(ye ~ xe-1, data=datae)
    beta1=rg[1]$coefficients
    
    RSS <- sum(residuals(rg)^2)
    ll[i] = RSS
    
    ## record the beta estimation
    betamat=rbind(betamat,as.numeric(beta1))
    betaold=beta
    ## moving average
    if(i>600){
      beta=matrix(apply(betamat[(i-600):i,],2,median),p^3,m)
      ll_bar[i]=median(ll[(i-600):i,])
    }
    else{
      if(i==iter){
        beta=matrix(apply(betamat[(i-600):i,],2,median),p^3,m)
        ll_bar[i]=median(ll[(i-600):i,])
      }
      else{
        beta=beta1
        ll_bar[i]=ll[i,]
      }
    }
    
    
    betamat1=rbind(betamat1,as.numeric(beta))
    
    if(!i%%400){print(i)}
    converge=0
    if (i>4000)
    {
      if (abs(ll_bar[i]-ll_bar[i-1])<stop)
      {
        converge=1
        break
      }
    }
  }
  return(list(beta,converge,ll[1:i],ll_bar[1:i],betamat,betamat1,i))
}


rep=1
R=6
ne_candi = c(8, 16, 24, 32, 40, 48, 56, 64)
MSE=matrix(rep(0),rep,length(ne_candi))#mean bias
PMAE=matrix(rep(0),rep,length(ne_candi))#prediction MSE
beta_all=array(rep(0),dim=c(rep,length(ne_candi),p^3))
time=matrix(rep(0),rep,length(ne_candi))

r=1
sigma=0.5
for (r in 1:rep)
{
  x=matrix(rnorm(p^3*n,0,1),n,p^3)
  y=x%*%beta0+matrix(rnorm(n,0,sd=sigma),n,1)
  data0=data.frame(x,y)
  x_p=matrix(rnorm(p^3*n,0,1),n,p^3)
  y_p=x_p%*%beta0+matrix(rnorm(n,0,sigma),n,1)
  
  for (i in 1:length(ne_candi)){
    
    start=Sys.time()
    result1=our(x,y,data0, ne_candi[i])
    end=Sys.time()
    t1=end-start
    time[r,i]=t1
    beta1=matrix(result1[[1]],p^3,1)
    
    #store estimated beta
    beta_all[r,i,]=matrix(beta1,1,p^3)
    
    #PMSE
    PMAE[r,i]=(sum(abs(x_p%*%beta1-y_p))/n)
    
    #sum of absolute error of parameter
    MSE[r,i]=sum((beta1-beta0)^2)/p^3
  }
  print(r)
}

# load("ne0_MSE.Rdata")
# load("ne0_PMAE.Rdata")
# load("ne0_beta_all.Rdata")

BMSE = matrix(rep(0),length(ne_candi), 1)
for (i in 1:length(ne_candi)){
  for (j in 1:rep){
    BMSE[i, 1] = BMSE[i, 1] + mean((beta0 - beta_all[j,i,])^2)^0.5
  }
  BMSE[i, 1] = BMSE[i, 1]/rep
}


num_zeros=matrix(rep(0),rep,length(ne_candi))
for (i in 1:rep){
  for (j in 1:length(ne_candi)){
    beta1=(beta_all[i,j,])
    beta1=as.tensor(array(beta1,dim=c(p,p,p)))
    egv1=tucker(beta1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
    thrsh=0.005
    num_zeros[i,j]=length(which((abs(egv1))<thrsh))
  }
}

sort(abs(egv1))
print(round((aPMAE=colMeans(PMAE)),4))
print(BMSE)
print(colMeans(num_zeros))

ne0.linear = c(round((aPMAE=colMeans(PMAE)),4), BMSE, colMeans(num_zeros))
save(ne0.linear,file="ne0_linear.Rdata")

save(MSE,file="ne0_MSE.Rdata")
save(PMAE,file="ne0_PMAE.Rdata")
save(beta_all,file="ne0_beta_all.Rdata")


#--- plot them -------------------------------
# Load required library
library(ggplot2)

true_zeros <- ne_candi  
experimented_zeros <- colMeans(num_zeros)

# Create a data frame
data <- data.frame(
  TrueZeros = true_zeros,
  ExperimentedZeros = experimented_zeros
)

# Create the plot
plot <- ggplot(data, aes(x = TrueZeros, y = ExperimentedZeros)) +
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +  # Add a regression line
  labs(
    x = "Size of ne",
    y = "Avg. Number of Zeros",
    title = ""  #True vs. Avg. Number of Zeros in Core Tensor
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 64), ylim = c(0, 64)) +  # Ensure equal axis lengths and set limits
  theme_minimal(base_size = 16) +  # Use a minimal theme with a larger base font size
  theme(
    axis.title = element_text(size = 22),  # Larger font size for axis titles
    axis.text = element_text(size = 20),  # Larger font size for axis text
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Centered title with larger font
  )

# Save the plot (optional, for high-quality output)
ggsave("ne0.pdf", plot = plot, width = 8, height = 6)

# Display the plot
print(plot)
