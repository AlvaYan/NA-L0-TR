#########################################################  linear
library(rTensor)
library(MASS)
#library(extraDistr)
library(glmnet)
library(BoomSpikeSlab)
library(Boom)
library(coda)
rm(list=ls())
set.seed(1)

n=200;p=4 ## n by p matrix , x???? n by p^3 matrix

setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\poi\\new8\\simulated data\\")
load("true_beta.Rdata") #beta0

# Function for k-fold cross-validation with parameter tuning for poisson.spike
cv_poisson_spike <- function(X, y, k = 5, n.ite = 1000) {
  # Simplified hyperparameter grid with three values for expected.model.size
  param_grid <- expand.grid(
    expected.model.size = c(5, 10, 15),   # Three values for expected model size
    prior.information.weight = c(0.1, 0.5, 0.9) # Influence of prior on coefficients
  )
  
  # Set up cross-validation folds
  set.seed(123) # Ensure reproducibility
  n <- nrow(X)
  folds <- sample(rep(1:k, length.out = n))  # Randomly assign data to k folds
  
  # Prepare a data frame to store CV results
  cv_results <- data.frame(param_grid, MAE = rep(NA, nrow(param_grid)))
  
  # Loop over each parameter combination
  for (i in 1:nrow(param_grid)) {
    # Extract the hyperparameters for this iteration
    expected_model_size <- param_grid$expected.model.size[i]
    prior_info_weight <- param_grid$prior.information.weight[i]
    
    # Store RMSE for each fold
    fold_mae <- numeric(k)
    
    for (j in 1:k) {
      # Split into training and validation sets
      train_indices <- which(folds != j)
      test_indices <- which(folds == j)
      
      X_train <- X[train_indices, ]
      y_train <- y[train_indices]
      X_test <- X[test_indices, ]
      y_test <- y[test_indices]
      
      # Define the custom prior for this fold and parameter combination
      prior <- PoissonZellnerPrior(
        predictors = X_train,
        counts = y_train,
        expected.model.size = expected_model_size,
        prior.information.weight = prior_info_weight
      )
      
      # Fit the Poisson spike-and-slab model on the training data
      model <- poisson.spike(
        y_train ~ X_train - 1,
        niter = n.ite,         # Number of MCMC iterations
        prior = prior           # Use the custom prior object
      )
      
      # Make predictions on the validation set
      beta_samples <- as.matrix(model$beta)
      last_samples <- tail(beta_samples, as.integer(n.ite * 0.4))
      thinned_samples <- last_samples[seq(1, nrow(last_samples), by = 2), ]
      beta_means <- colMeans(thinned_samples)
      
      log_lambda_pred <- X_test %*% beta_means       # Linear predictor for the Poisson rate
      lambda_pred <- exp(log_lambda_pred)            # Predicted Poisson rate
      
      # Calculate RMSE for the validation fold
      fold_mae[j] <- mean(abs(y_test-as.numeric(lambda_pred)))
    }
    
    # Store the average RMSE across folds for the current hyperparameter combination
    cv_results$MAE[i] <- mean(fold_mae)
  }
  
  # Return the cross-validation results
  list(
    results = cv_results,
    best_params = cv_results[which.min(cv_results$MAE), ]
  )
}


rep=200
R=6
MSE_poi=matrix(rep(0),rep,1)#mean bias
PMAE_poi=matrix(rep(0),rep,1)#classification rate
beta_poi_all=array(rep(0),dim=c(rep,1,p^3))
time=matrix(rep(0),rep,1)

r=1
r=r+1
for (r in 1:rep)
{
  setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\poi\\new8\\simulated data\\")
  name = paste0("x", r, ".Rdata")  # "x2.Rdata"
  load(name) # x
  
  name = paste0("y", r, ".Rdata")  # "y2.Rdata"
  load(name) # y
  
  # data0=data.frame(x,y)
  #get the OLS estimator for beta without regularization
  start=Sys.time()
  cv_results <- cv_poisson_spike(X = x, y = y, n.ite = 1000, k = 5)
  
  # Extract the best hyperparameters from cross-validation results
  best_params <- cv_results$best_params
  expected_model_size <- best_params$expected.model.size
  prior_info_weight <- best_params$prior.information.weight
  
  # Define the final prior using the best parameters
  final_prior <- PoissonZellnerPrior(
    predictors = x,                      # Full predictor matrix
    counts = y,                      # Full response vector
    expected.model.size = expected_model_size,
    prior.information.weight = prior_info_weight
  )
  
  # Fit the final lm.spike model with the selected prior
  rg1 <- poisson.spike(
    y ~ x - 1,
    niter = 10000,              # Use more iterations for the final model for better convergence
    prior = final_prior
  )
  end=Sys.time()
  
  # # trace plot for convergence
  # i = 2
  # plot(rg1$beta[, i], type = "l", main = paste("Trace plot for beta[", i, "]"),
  #      xlab = "Iteration", ylab = paste("beta[", i, "]"))
  # # effectiveSize, close to real size is better
  # mcmc_samples = as.matrix(rg1$beta)[, i]
  # mcmc_obj <- as.mcmc(mcmc_samples)
  # ess <- effectiveSize(mcmc_obj)
  # print(paste("Effective Sample Size (ESS):", round(ess, 2)))
  # # |Z| < 2 indicate stationary
  # geweke_diag <- geweke.diag(mcmc_obj)
  # print(geweke_diag)
  # # test for stationarity and accuracy of the estimated mean of the parameter
  # heidel_diag <- heidel.diag(mcmc_obj)
  # print(heidel_diag)
  
  t1=end-start
  
  beta_samples <- as.matrix(rg1$beta)
  last_4000_samples <- tail(beta_samples, 4000)
  thinned_samples <- last_4000_samples[seq(1, nrow(last_4000_samples), by = 2), ]
  beta <- colMeans(thinned_samples)
  # res = summary(rg1)
  # beta = res$coefficients[,1]
  # beta <- beta[order(as.numeric(gsub("x", "", names(beta))))]
  
  #record time
  time[r,1]=t1
  beta1=matrix(beta,p^3,1)
  #store log likelihood
  #LL[r,1]=result1[[3]]
  #store estimated beta
  beta_poi_all[r,1,]=matrix(beta1,1,p^3)
  #RMSE
  #predict(rg1, newdata = x)
  MSE_poi[r,1]=sum((beta1-beta0)^2)/p^3
  
  #PMSE
  setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\poi\\new8\\simulated data\\")
  name = paste0("xp", r, ".Rdata")  # "x2.Rdata"
  load(name)
  name = paste0("yp", r, ".Rdata")  # "y2.Rdata"
  load(name)
  
  mu1=x_p%*%beta1
  pr1=exp(mu1)
  PMAE_poi[r,1]=mean(abs(pr1-as.numeric(y_p)))
  
  print(r)
}

#save data after regression
save(MSE_poi,file="MSE_poi.Rdata")
save(PMAE_poi,file="PMAE_poi.Rdata")
save(beta_poi_all,file="beta_poi_all.Rdata")
#save(betamle_poi_all,file="betamle_poi_all.Rdata")
save(time,file="time_poi.Rdata")

# round((MSE=colMeans(MSE_poi)),2)
# round((cMSE=colMeans(cMSE_poi)),3)
# round((RMSPE=colMeans(PMAE_poi)),5)
# round((DV_poi=colMeans(DV_poi)),3)


BMSE = 0
for (i in 1:rep){
  BMSE = BMSE + mean((beta0 - beta_poi_all[i,,])^2)^0.5
}

num_zeros=matrix(0,rep,1)
for (i in 1:rep){
  beta1=(beta_poi_all[i,1,])
  #value unadjusted order
  beta1=as.tensor(array(beta1,dim=c(p,p,p)))
  
  egv1=tucker(beta1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
  
  #Calulate false postive and false negative for non-zero core-tensor
  thrsh=0.05
  num_zeros[i,1]=length(which((abs(egv1))<thrsh))
}

print(round((aPMAE=colMeans(PMAE_poi)),4))
print(BMSE/rep)
print(mean(num_zeros))

res.poi = c(round((aPMAE=colMeans(PMAE_poi)),4), BMSE/rep, mean(num_zeros))
save(res.poi,file="spike_poi.Rdata")

#-- reload data and analysis --------
setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\poi\\new8\\simulated data\\")

load("PMAE_poi.Rdata")
load("beta_poi_all.Rdata")

BMSE = 0
for (i in 1:rep){
  BMSE = BMSE + mean((beta0 - beta_poi_all[i,,])^2)^0.5
}

num_zeros=matrix(0,rep,1)
for (i in 1:rep){
  beta1=(beta_poi_all[i,1,])
  #value unadjusted order
  beta1=as.tensor(array(beta1,dim=c(p,p,p)))
  
  egv1=tucker(beta1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
  
  #Calulate false postive and false negative for non-zero core-tensor
  thrsh=0.05
  num_zeros[i,1]=length(which((abs(egv1))<thrsh))
}

print(round((aPMAE_poi=colMeans(PMAE_poi)),4))
column_sd <- apply(PMAE_poi, 2, sd)
print(column_sd)
print(BMSE/rep)
print(mean(num_zeros))
