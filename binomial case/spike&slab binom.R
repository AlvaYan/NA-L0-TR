#########################################################  linear
library(rTensor)
library(MASS)
library(glmnet)
library(BoomSpikeSlab)
library(Boom)
library(coda)
rm(list=ls())
set.seed(1)

n=300;p=4 ## n by p matrix , x???? n by p^3 matrix


setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\binom\\new16\\simulated data\\")
load("true_beta.Rdata") #beta0

cv_logit_spike <- function(X, y, n.ite = 1000, k = 5) {
  # Define hyperparameter grid based on prior explanations
  param_grid <- expand.grid(
    expected.model.size = c(3, 11, 32),      # Expected number of predictors
    prior.information.weight = c(0.1, 0.5, 0.9) # Influence of prior on coefficients
  )
  
  # Set up cross-validation folds
  set.seed(123) # Ensure reproducibility
  n <- nrow(X)
  folds <- sample(rep(1:k, length.out = n))  # Randomly assign data to k folds
  
  # Prepare a data frame to store CV results
  cv_results <- data.frame(param_grid, LogLoss = rep(NA, nrow(param_grid)))
  
  # Loop over each parameter combination
  i=1
  for (i in 1:nrow(param_grid)) {
    # Extract the hyperparameters for this iteration
    expected_model_size <- param_grid$expected.model.size[i]
    prior_info_weight <- param_grid$prior.information.weight[i]
    
    # Store RMSE for each fold
    fold_logloss <- numeric(k)
    
    j=1
    for (j in 1:k) {
      # Split into training and validation sets
      train_indices <- which(folds != j)
      test_indices <- which(folds == j)
      
      X_train <- X[train_indices, ]
      y_train <- y[train_indices]
      X_test <- X[test_indices, ]
      y_test <- y[test_indices]
      
      # Define the custom prior for this fold and parameter combination
      prior <- LogitZellnerPrior(
        predictors = X_train,
        successes = y_train,
        expected.model.size = expected_model_size,
        prior.information.weight = prior_info_weight
      )
      
      # Fit the model on the training data
      model <- logit.spike(
        y_train ~ X_train - 1,
        niter = n.ite,                    # Number of MCMC iterations
        prior = prior
      )
      
      # Make predictions on the validation set
      beta_samples <- as.matrix(model$beta)
      last_samples <- tail(beta_samples, as.integer(n.ite * 0.4))
      thinned_samples <- last_samples[seq(1, nrow(last_samples), by = 2), ]
      beta_means <- colMeans(thinned_samples)
      
      logits <- X_test %*% beta_means
      predictions <- 1 / (1 + exp(-logits))  # Convert logits to probabilities
      
      # Calculate Log-Loss for the validation fold
      fold_logloss[j] <- -mean(y_test * log(predictions) + (1 - y_test) * log(1 - predictions))
      
    }
    
    # Store the average RMSE across folds for the current hyperparameter combination
    cv_results$LogLoss[i] <- mean(fold_logloss)
  }
  
  # Return the cross-validation results
  list(
    results = cv_results,
    best_params = cv_results[which.min(cv_results$LogLoss), ]
  )
}


rep=200
R=6

MSE_bi=matrix(rep(0),rep,1)
CR_bi=matrix(rep(0),rep,1)
beta_bi_all=array(rep(0),dim=c(rep,1,p^3))

time=matrix(rep(0),rep,1)

r=1
r=r+1
for (r in 1:rep)
{
  setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\binom\\new16\\simulated data\\")
  name = paste0("x", r, ".Rdata")  # "x2.Rdata"
  load(name) # x
  
  name = paste0("y", r, ".Rdata")  # "y2.Rdata"
  load(name) # y
  
  # data0=data.frame(x,y)
  #get the OLS estimator for beta without regularization
  
  start=Sys.time()
  cv_results <- cv_logit_spike(X = x, y = y, n.ite = 1000, k = 5)
  
  # Extract the best hyperparameters from cross-validation results
  best_params <- cv_results$best_params
  expected_model_size <- best_params$expected.model.size
  prior_info_weight <- best_params$prior.information.weight
  
  # Define the final prior using the best parameters
  final_prior <- LogitZellnerPrior(
    predictors = x,                      # Full predictor matrix
    successes = y,                      # Full response vector
    expected.model.size = expected_model_size,
    prior.information.weight = prior_info_weight
  )
  
  # Fit the final lm.spike model with the selected prior
  rg1 <- logit.spike(
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
  
  #record time
  time[r,1]=t1
  beta1=matrix(beta,p^3,1)
  time[r,1]=t1
  beta_bi_all[r,1,]=matrix(beta1,1,p^3)
  MSE_bi[r,1]=sum((beta1-beta0)^2)/p^3
  
  # #RMSE
  # mu=x%*%beta
  # pr=1/(1+exp(-mu))
  # y_hat=as.numeric(pr>0.5)
  # RMSE[r,1]=sum((y_hat-as.numeric(y))==0)/n
  
  #PMSE
  setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\binom\\new16\\simulated data\\")
  name = paste0("xp", r, ".Rdata")  # "x2.Rdata"
  load(name)
  name = paste0("yp", r, ".Rdata")  # "y2.Rdata"
  load(name)
  
  mu1=x_p%*%beta1
  pr1=1/(1+exp(-mu1))
  y_1=as.numeric(pr1>0.5)
  CR_bi[r,1]=mean((y_1-as.numeric(y_p))==0)
  
  print(r)
}


#save data after regression
save(MSE_bi,file="MSE_bi.Rdata")
save(CR_bi,file="CR_bi.Rdata")
save(beta_bi_all,file="beta_bi_all.Rdata")
#save(betamle_bi_all,file="betamle_bi_all.Rdata")
save(time,file="time_bi.Rdata")

# round((aMSE=colMeans(MSE_bi)),4)
# round((acMSE=colMeans(cMSE_bi)),4)
# round((aMCR=1-colMeans(CR_bi)),4)
# round((aDV_bi=colMeans(DV_bi)),4)

BMSE = 0
for (i in 1:rep){
  BMSE = BMSE + mean((beta0 - beta_bi_all[i,,])^2)^0.5
}

num_zeros=matrix(0,rep,1)
for (i in 1:rep){
  beta1=(beta_bi_all[i,1,])
  #value unadjusted order
  beta1=as.tensor(array(beta1,dim=c(p,p,p)))
  
  egv1=tucker(beta1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
  
  #Calulate false postive and false negative for non-zero core-tensor
  thrsh=0.05
  num_zeros[i,1]=length(which((abs(egv1))<thrsh))
}

print(round((aMCR=1-colMeans(CR_bi)),4))
print(BMSE/rep)
print(mean(num_zeros))

res.binom = c(round((aMCR=1-colMeans(CR_bi)),4), BMSE/rep, mean(num_zeros))
save(res.binom, file="spike_binom.Rdata")

#-- reload data and analysis --------
setwd("C:\\TianLaptopData\\LaptopData\\desktop\\simulation new\\adressing cmts, rd1\\binom\\new16\\simulated data\\")


#load("MSE_bi.Rdata")
load("CR_bi.Rdata")
load("beta_bi_all.Rdata")

BMSE = 0
for (i in 1:rep){
  BMSE = BMSE + mean((beta0 - beta_bi_all[i,,])^2)^0.5
}

num_zeros=matrix(0,rep,1)
for (i in 1:rep){
  beta1=(beta_bi_all[i,1,])
  #value unadjusted order
  beta1=as.tensor(array(beta1,dim=c(p,p,p)))
  
  egv1=tucker(beta1,ranks=c(p,p,p) ,tol=1e-6,max_iter = 50)$Z@data
  
  #Calulate false postive and false negative for non-zero core-tensor
  thrsh=0.05
  num_zeros[i,1]=length(which((abs(egv1))<thrsh))
}

MCR = 1-(CR_bi)
print(round(colMeans(MCR),4))
column_sd <- apply(MCR, 2, sd)
print(column_sd)
print(BMSE/rep)
print(mean(num_zeros))
