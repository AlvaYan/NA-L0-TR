#########################################################  linear
library(rTensor)
library(MASS)
#library(extraDistr)
library(glmnet)
rm(list=ls())
set.seed(1)

n=300;p=4 ## n by p matrix , x???? n by p^3 matrix

setwd("C:\\Users\\alvay\\Downloads\\linearDV\\")
load("true_beta.Rdata") #beta0
#beta_0=as.tensor(array(beta0,dim=c(p,p,p)))
beta0=matrix(rep((rnorm(p^3/2,0,1)),2),p^3,1)

stop=0.01

#set up storage matrix

n_smooth = 600
our=function(x,y,data0, n_smooth)
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
  ne=32 #floor(p^3-2)
  m=1
  ne2=0
  multi=1
  lambda2=8
  sprs=F
  lambda=50/multi
  iter=30000
  betamat=NULL
  betamat1=NULL
  eigvalmat=NULL
  eigvalmat1=NULL
  ll=matrix(rep(0),iter,1)
  ll_bar=matrix(rep(0),iter,1)
  beta2=betaols
  i=1
  for(i in 1:iter){
    noise=NULL
    
    ## beta is stored in 2d matrix, therefore to use it as a tensor, we need to fold it first
    beta_1=fold(matrix(beta,p,p^2),row_idx=2,col_idx=c(3,1),modes=c(p,p,p))
    ## apply tucker decomposition to get the core tensor and unitary matrices
    #set.seed(1)
    beta_1_dec=tucker(beta_1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 20)
    ## treat the core tnsor as eigen values and unitary matrices as eigenvectors
    eigval1=beta_1_dec$Z@data
    eigvec=beta_1_dec$U
    
    eigvalmat=rbind(eigvalmat,eigval1)
    
    if(i>n_smooth){
      sign=apply(sign( eigvalmat[(i-n_smooth):i,]),2,median)
      eigval2=matrix(apply( abs(eigvalmat[(i-n_smooth):i,]),2,median),p^3,m)
      #sign2=sign(eigval2)
      eigval=eigval2*sign
    }else{
      #sign=sign(beta_0_eval)
      #sign2=sign(eigval1)
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
    
    ## record the beta estimation
    betamat=rbind(betamat,as.numeric(beta1))
    betaold=beta
    
    # record deviance
    # null_model <- lm(ye ~ 1)
    RSS <- sum(residuals(rg)^2)
    ll[i] = RSS
    
    ## moving average
    if(i>n_smooth){
      beta=matrix(apply(betamat[(i-n_smooth):i,],2,median),p^3,m)
      ll_bar[i]=median(ll[(i-n_smooth):i,])
    }
    else{
      if(i==iter){
        beta=matrix(apply(betamat[(i-n_smooth):i,],2,median),p^3,m)
        ll_bar[i]=median(ll[(i-n_smooth):i,])
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

rep = 1
MSE=matrix(rep(0),rep,3)#mean bias
PMAE=matrix(rep(0),rep,3)#prediction MSE
beta_all=array(rep(0),dim=c(rep,3,p^3))

r=1
sigma=0.5
#setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\linear\\new6\\simulated data\\")
x=matrix(rnorm(p^3*n,0,1),n,p^3)
y=x%*%beta0+matrix(rnorm(n,0,sd=sigma),n,1)
data0=data.frame(x,y)
x_p=matrix(rnorm(p^3*n,0,1),n,p^3)
y_p=x_p%*%beta0+matrix(rnorm(n,0,sigma),n,1)
# name = paste0("x", r, ".Rdata")  # "x2.Rdata"
# load(name) # x
# name = paste0("y", r, ".Rdata")  # "y2.Rdata"
# load(name) # y
# name = paste0("xp", r, ".Rdata")  # "x2.Rdata"
# load(name)
# name = paste0("yp", r, ".Rdata")  # "y2.Rdata"
# load(name)
# data0=data.frame(x,y)

rg = lm(y ~ x-1, data=data0)
#do regression using 4 types of method
start=Sys.time()
result1=our(x,y,data0, 600)
result2=our(x,y,data0, 300)
result3=our(x,y,data0, 100)
end=Sys.time()
t1=end-start

#calculate MSE, PMSE, check number of zeros in core tensor
beta1=matrix(result1[[1]],p^3,1)
beta_1=as.tensor(array(beta1,dim=c(p,p,p)))
tucker(beta_1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
beta_all[r,1,]=matrix(beta1,1,p^3)
(PMAE[r,1]=mean(abs(x_p%*%beta1-y_p)))
(MSE[r,1]=sum((beta1-beta0)^2)/p^3)
# check performance for ols case
betaols=glm(y ~ x-1, data=data0,family='gaussian')[1]$coefficients
beta_ols=as.tensor(array(betaols,dim=c(p,p,p)))
tucker(beta_ols,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
sum((betaols-beta0)^2)/p^3
#check performance for our method with result from first iteration
betaite1 = matrix(result1[[6]][1,],p^3,1)
beta_ite1=as.tensor(array(betaite1,dim=c(p,p,p)))
mean(abs(x_p%*%betaite1-y_p))
tucker(beta_ite1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
sum((betaite1-beta0)^2)/p^3


beta2=matrix(result2[[1]],p^3,1)
beta_all[r,2,]=matrix(beta2,1,p^3)
PMAE[r,2]=mean(abs(x_p%*%beta2-y_p))
MSE[r,2]=sum((beta2-beta0)^2)/p^3

beta3=matrix(result3[[1]],p^3,1)
beta_all[r,3,]=matrix(beta3,1,p^3)
PMAE[r,3]=mean(abs(x_p%*%beta3-y_p))
MSE[r,3]=sum((beta3-beta0)^2)/p^3

# log likelihood
LL1=result1[[3]]
LL_bar1=result1[[4]]
LL2=result2[[3]]
LL_bar2=result2[[4]]
LL3=result3[[3]]
LL_bar3=result3[[4]]

#save data after regression
#setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\linear\\new6\\simulated data\\")
save(LL1,file="LL1.Rdata")
save(LL_bar1,file="LL_bar1.Rdata")
save(LL2,file="LL2.Rdata")
save(LL_bar2,file="LL_bar2.Rdata")
save(LL3,file="LL3.Rdata")
save(LL_bar3,file="LL_bar3.Rdata")

save(MSE,file="DV_MSE.Rdata")
save(PMAE,file="DV_PMAE.Rdata")
save(beta_all,file="DV_beta_all.Rdata")

round((aPMAE=colMeans(PMAE)),4)
round((aMSE=colMeans(MSE)),4)

#----------residual deviance plt--------------------------------------------
# Load ggplot2 package
library(ggplot2)

# Sample data with different numbers of iterations (replace with your actual data)
iterations_1 <- 1:length(LL_bar1)
iterations_2 <- 1:length(LL_bar2)
iterations_3 <- 1:length(LL_bar3)

# Combine data into a single data frame, including the specific iteration number for each sequence
data <- data.frame(
  Iteration = c(iterations_1, iterations_2, iterations_3),
  Residual_Deviance = c(LL_bar1, LL_bar2, LL_bar3),
  Sequence = rep(c("m=600", "m=300", "m=100"),
                 times = c(length(iterations_1), length(iterations_2), length(iterations_3)))
)

# Create the plot with different lengths for each sequence
plot <- ggplot(data, aes(x = Iteration, y = Residual_Deviance, color = Sequence)) +
  geom_line(size = 1) +  # Line plot for each sequence
  labs(
    title = "",  # Residual Deviance Convergence for Linear TR
    x = "Iteration Number",
    y = "Residual Deviance"
  ) +
  scale_color_manual(values = c("blue", "red", "green")) +  # Custom colors for clarity
  theme_minimal(base_size = 14) +  # Clean theme for publication
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and increase title size
    axis.title = element_text(size = 22),               # Larger axis title size for readability
    axis.text = element_text(size = 20),                # Larger axis text size for readability
    legend.title = element_blank(),                     # Remove legend title
    legend.text = element_text(size = 14),              # Increase legend text size
    legend.position = "top"                             # Position legend at the top
  )


plot
# Save the plot as PNG (high resolution)
ggsave("DV_multiple_windows_linear.png", plot, width = 7, height = 5, dpi = 300)

#------ spaggetti plot -------------------
# Example: Generating a matrix of parameter estimates (replace this with your data)
set.seed(123)
N <- dim(result1[[6]])[1] # Number of iterations
P <- dim(result1[[6]])[2]   # Number of parameters
parameter_estimates <- result1[[6]]
#parameter_estimates[,1]

# Convert the matrix to a data frame in long format for plotting
library(ggplot2)
library(reshape2)

df <- as.data.frame(parameter_estimates)
df$Iteration <- 1:N
df_long <- melt(df, id.vars = "Iteration", variable.name = "Parameter", value.name = "Estimate")

# Determine a suitable y-axis range to emphasize close values
y_min <- min(df_long$Estimate) - 0.1  # Slightly below the min
y_max <- max(df_long$Estimate) + 0.1  # Slightly above the max

# Create the spaghetti plot with adjusted y-axis limits
plot <- ggplot(df_long, aes(x = Iteration, y = Estimate, color = Parameter, group = Parameter)) +
  geom_line(alpha = 0.7) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(
    title = "",  # Plot of Parameter Estimates Over Iterations for Linear TR
    x = "Iteration",
    y = "Parameter Estimate"
  ) +
  theme_minimal(base_size = 16) +  # Increase base font size for better readability
  theme(
    legend.position = "none",                              # Remove legend
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),  # Larger title font size
    axis.title = element_text(face = "bold", size = 22),   # Larger axis title font size
    axis.text = element_text(size = 20)                   # Larger tick label font size
  )

# Save the plot as a high-resolution PNG
ggsave("spaghetti_parameters_linear.png", plot = plot, width = 8, height = 6, dpi = 300)

plot  # Display the plot

#------------ trace plot
set.seed(123)
N <- dim(result1[[6]])[1] # Number of iterations
P <- dim(result1[[6]])[2]   # Number of parameters
parameter_estimates <- result1[[6]]

# Calculate differences between consecutive iterations for each parameter
parameter_diffs <- apply(parameter_estimates, 2, diff)
parameter_diffs <- as.data.frame(parameter_diffs)
parameter_diffs$Iteration <- 2:N

# Reshape data for plotting
library(ggplot2)
library(reshape2)
df_long <- melt(parameter_diffs, id.vars = "Iteration", variable.name = "Parameter", value.name = "Estimate_Diff")

# Create the plot
plot <- ggplot(df_long, aes(x = Iteration, y = Estimate_Diff, color = Parameter, group = Parameter)) +
  geom_line(alpha = 0.7) +
  labs(
    title = "",  # Convergence Behavior of Parameter Estimates for Linear TR
    x = "Iteration",
    y = "Change in Parameter Estimate"
  ) +
  theme_minimal(base_size = 16) +  # Increase base font size for overall readability
  theme(
    legend.position = "none",                                # Remove legend
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),  # Larger title font
    axis.title = element_text(face = "bold", size = 22),     # Larger axis titles
    axis.text = element_text(size = 20)                     # Larger axis tick labels
  )

# Save the plot
ggsave("convergence_parameters_linear.png", plot = plot, width = 8, height = 6, dpi = 300)

plot  # Display the plot
