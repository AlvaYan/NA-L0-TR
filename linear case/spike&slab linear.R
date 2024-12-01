#########################################################  linear
library(rTensor)
library(MASS)
library(glmnet)
library(BoomSpikeSlab)
library(Boom)
library(coda)
rm(list=ls())
set.seed(1)

n=300;p=4 ## n by p matrix

setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\linear\\new6\\simulated data\\")
load("true_beta.Rdata") #beta0

# Function for k-fold cross-validation with parameter tuning for lm.spike
cv_lm_spike <- function(X, y, n.ite = 1000, k = 5) {
  # Define hyperparameter grid based on prior explanations
  param_grid <- expand.grid(
    expected.model.size = c(3, 11, 19),      # Expected number of predictors
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
      prior <- SpikeSlabPrior(
        x = X_train,
        y = y_train,
        expected.model.size = expected_model_size,
        prior.information.weight = prior_info_weight
      )
      
      # Fit the model on the training data
      model <- lm.spike(
        y_train ~ X_train - 1,
        niter = n.ite,                    # Number of MCMC iterations
        prior = prior
      )
      
      # Make predictions on the validation set
      beta_samples <- as.matrix(model$beta)
      last_samples <- tail(beta_samples, as.integer(n.ite * 0.4))
      thinned_samples <- last_samples[seq(1, nrow(last_samples), by = 2), ]
      beta_means <- colMeans(thinned_samples)
      predictions <- X_test %*% beta_means
      
      # Calculate RMSE for the validation fold
      fold_mae[j] <- mean(abs(y_test - predictions))
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
RMSE=matrix(rep(0),rep,1)
MSE=matrix(rep(0),rep,1)#mean bias
PMAE=matrix(rep(0),rep,1)#prediction MSE
beta_all=array(rep(0),dim=c(rep,1,p^3))
#betaols_all=matrix(rep(0),rep,p^3)
time=matrix(rep(0),rep,1)

r=1
for (r in 1:rep)
{
  setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\linear\\new6\\simulated data\\")
  name = paste0("x", r, ".Rdata")  # "x2.Rdata"
  load(name) # x
  
  name = paste0("y", r, ".Rdata")  # "y2.Rdata"
  load(name) # y
  
  start=Sys.time()
  
  cv_results <- cv_lm_spike(X = x, y = y, n.ite = 1000, k = 5)
  
  # Extract the best hyperparameters from cross-validation results
  best_params <- cv_results$best_params
  expected_model_size <- best_params$expected.model.size
  prior_info_weight <- best_params$prior.information.weight
  
  # Define the final prior using the best parameters
  final_prior <- SpikeSlabPrior(
    x = x,                      # Full predictor matrix
    y = y,                      # Full response vector
    expected.model.size = expected_model_size,
    prior.information.weight = prior_info_weight
  )
  
  # Fit the final lm.spike model with the selected prior
  rg1 <- lm.spike(
    y ~ x - 1,
    niter = 10000,              # Use more iterations for the final model for better convergence
    prior = final_prior
  )
  end=Sys.time()
  
  t1=end-start
  beta_samples <- as.matrix(rg1$beta)
  last_4000_samples <- tail(beta_samples, 4000)
  thinned_samples <- last_4000_samples[seq(1, nrow(last_4000_samples), by = 2), ]
  beta <- colMeans(thinned_samples)
  
  #record time
  time[r,1]=t1
  beta1=matrix(beta,p^3,1)
  beta_all[r,1,]=matrix(beta1,1,p^3)
  RMSE[r,1]=(sum((x%*%beta1-y)^2)/n)^0.5
  setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\linear\\new6\\simulated data\\")
  name = paste0("xp", r, ".Rdata")  # "x2.Rdata"
  load(name)
  name = paste0("yp", r, ".Rdata")  # "y2.Rdata"
  load(name)
  PMAE[r,1]=mean(abs(x_p%*%beta1-y_p))
  MSE[r,1]=sum((beta1-beta0)^2)/p^3
  
  print(r)
}

#save data after regression
save(PMAE,file="PMAE.Rdata")
save(beta_all,file="beta_all.Rdata")
save(time,file="time_linear.Rdata")


BMSE = 0
for (i in 1:rep){
  BMSE = BMSE + mean((beta0 - beta_all[i,,])^2)^0.5
}

num_zeros=matrix(0,rep,1)
for (i in 1:rep){
  beta1=(beta_all[i,1,])
  #value unadjusted order
  beta1=as.tensor(array(beta1,dim=c(p,p,p)))
  
  egv1=tucker(beta1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
  
  thrsh=0.005
  num_zeros[i,1]=length(which((abs(egv1))<thrsh))
}

print(round((aPMAE=colMeans(PMAE)),4))
print(BMSE/rep)
print(mean(num_zeros))

res.linear = c(round((aPMAE=colMeans(PMAE)),4), BMSE/rep, mean(num_zeros))
save(res.linear,file="spike_linear.Rdata")

#-- reload data and analysis --------
setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\linear\\new6\\simulated data\\")

load("PMAE.Rdata")
load("beta_all.Rdata")
load("time_linear.Rdata")

BMSE = 0
for (i in 1:rep){
  BMSE = BMSE + mean((beta0 - beta_all[i,,])^2)^0.5
}

num_zeros=matrix(0,rep,1)
for (i in 1:rep){
  beta1=(beta_all[i,1,])
  beta1=as.tensor(array(beta1,dim=c(p,p,p)))
  
  egv1=tucker(beta1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
  
  thrsh=0.005
  num_zeros[i,1]=length(which((abs(egv1))<thrsh))
}

print(round((aPMAE=colMeans(PMAE)),4))
column_sd <- apply(PMAE, 2, sd)
print(column_sd)
print(BMSE/rep)
print(mean(num_zeros))
