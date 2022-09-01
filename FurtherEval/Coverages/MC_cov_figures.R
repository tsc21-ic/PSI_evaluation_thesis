## Load required packages
# install.packages("pacman")
pacman::p_load(selectiveInference, prospectr, knockoff, mvtnorm, glmnet, stabs, MASS, dplyr,ggthemes, kableExtra, ggplot2, tibble, tidyverse,
               viridis, pbapply, data.table, progress, randomcoloR, cowplot, ggpubr,tcltk, ggh4x, patchwork, profvis, fGarch, mvtnorm, beepr)
source("fun.R")

set.seed(2)
SampleNormal = function(beta, X, sigma) X%*%beta + rnorm(dim(X)[1], sd = sigma)

## Simulation setup
iters = 1000 ## number of replications, R to record
times = 5 ## times to run MC
n = 100; p = 3; rho = 0.1
sigma = 1; f=0.5; level=0.9; alpha=0.1
beta = c(1,0.5,0.2)
estimateVar=FALSE   ## Assume known variance


## Change to B = 200 (first figure) ###################
B = 200

## Change to B = 1000 (second figure) ###################
# B = 1000  ## number condition on selection samples to simulate


###################
## Time the running
pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = iters*times,
                      complete = "=",
                      incomplete = "-",
                      current = ">",   
                      clear = FALSE,   
                      width = 100)

## Generate design matrix, X
Sigma_X = matrix(nrow=p, ncol=p)
for (i in 1:p){
  for (j in 1:p){
    Sigma_X[i, j] = rho^(abs(i - j))
  }
}
##########################################
my_results = list()
for(time in 1:times){
  ## To store MC data
  beta_poly = c()
  coverage_poly_theory=c()
  coverage_poly_ori = c()
  coverage_poly_ori_R= c()
  coverage_poly_ori_R_more = c()
  beta_obs_full = c()
  iter_count=c()
  
  for (i in 1:iters){
    pb$tick()
    X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
    y = SampleNormal(beta, X, sigma = sigma)
    
    ## For estimated sigma
    est_sig = sigma_hat(y, X)
    
    ############# MC ##########################
    ## Carry out LASSO to determine sub-model
    initial = Selection(y, X, type="lasso.fixed.lambda", p = p, first=TRUE, signs=FALSE) ## returns lambda, coefficient vector, and variables
    vars_initial = initial$vars
    vars_m = which(vars_initial)     ## selected variables' number for the sub-model
    X_m = X[,vars_m, drop=FALSE]
    print(paste0("Selected variables are: ", vars_m))
    
    # Actual values of submodel coefficients
    mu = X %*% beta
    eta = X_m %*% (solve(t(X_m) %*% X_m))
    beta_obs = t(eta) %*% mu
    
    est_sigma = sigma_hat(y, X)   ## estimate variance case
    
    for (j in 1:length(vars_m)){
      cat(sprintf("Inference for variable %i in iteration %i \n",vars_m[j], i))
      
      ## Calculate residual
      eta_j = eta[,j, drop=FALSE]   ## want to test the variable j in X_m
      denom = sum(t(eta_j) %*% eta_j)
      resid = (diag(n) - ((eta_j/denom) %*% t(eta_j))) %*% y  ## residual
      
      ## Actual projection (partial regression) parameter value
      beta_obs_j = beta_obs[j]
      
      ## target parameter estimate
      betahat_obs = t(eta_j) %*% y
      
      if(estimateVar){
        Sigma2 = sqrt(sum((est_sigma$sigmahat^2)*(t(eta_j) %*% eta_j)))   ## variance of estimate
      }else{
        Sigma2 = sqrt(sum((sigma^2)*(t(eta_j) %*% eta_j)))
      }
      
      ## Fix the actual parameter value
      pivot = function(param=beta_obs_j, randomization=FALSE, morerandomization=FALSE){
        betahat_sim_all = numeric(B)
        for(b in 1:B){
          condition = FALSE
          while(condition == FALSE){
            ## Generate beta_hat distribution
            betahat_sim = rnorm(1, mean = param, sd = Sigma2)
            y_sim = resid + ((eta_j %*% betahat_sim) / denom)
            if(randomization){
              y_sim = y_sim + rnorm(n, sd=est_sig$sigmahat)
            }
            if(morerandomization){
              y_sim = y_sim + rnorm(n, sd=2.5*est_sig$sigmahat)
            }
            selected_vars = Selection(y_sim, X, type="lasso.fixed.lambda", p = p, first=FALSE, lam=initial$lambda, signs = FALSE)   ## Use the same lambda since it is a function of X and errors only
            if(all(selected_vars == vars_initial) & !is.na(selected_vars[vars_m[j]])){
              condition=TRUE
            }
          }
          betahat_sim_all[b] = betahat_sim
        }
        lower_R = quantile(betahat_sim_all, probs = alpha/2)
        upper_R = quantile(betahat_sim_all, probs = 1- alpha/2)
        coverage_MC = (lower_R < betahat_obs)*(upper_R > betahat_obs)    ## check whether greater than the parameter value or not: if 0 => do not cover
        return(coverage_MC)
      }
      coverage_poly_ori = c(coverage_poly_ori, pivot(randomization=FALSE))
      coverage_poly_ori_R = c(coverage_poly_ori_R, pivot(randomization=TRUE))
      coverage_poly_ori_R_more = c(coverage_poly_ori_R_more, pivot(morerandomization=TRUE))
    }
    
    ######################## Theoretical comparison ##################
    ## actual polyhedral lemma
    theory_res = fixedLassoInf(X, y, beta = initial$coefs, lambda = initial$lambda, sigma = sigma, intercept = FALSE, type="partial")
    ## Extract the CIs from theoretical
    CI = theory_res$ci
    ## Store coverages yes or no for the iteration i
    coverage_poly_theory = c(coverage_poly_theory, (CI[, 1] < beta_obs) * (CI[, 2] > beta_obs))
    beta_poly = c(beta_poly, abs(beta[vars_m]))
    beta_obs_full = c(beta_obs_full, beta_obs)
    iter_count = c(iter_count, rep(i,length(beta_obs)))
  }
  my_results[[time]]=list(beta_poly=beta_poly,
                          coverage_poly_theory=coverage_poly_theory,
                          coverage_poly_ori=coverage_poly_ori,
                          coverage_poly_ori_R=coverage_poly_ori_R,
                          coverage_poly_ori_R_more=coverage_poly_ori_R_more,
                          beta_obs_full = beta_obs_full,
                          iter = iter_count)
}

# save(my_results, file= "MC_R_1000_B_200.Rdata")
# save(my_results, file= "MC_R_1000_B_1000.Rdata")