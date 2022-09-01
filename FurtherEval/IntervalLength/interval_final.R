pacman::p_load(selectiveInference, prospectr, knockoff, mvtnorm, glmnet, stabs, MASS, dplyr,ggthemes, kableExtra, ggplot2, tibble, tidyverse,
               viridis, pbapply, data.table, progress, randomcoloR, cowplot, ggpubr,tcltk, ggh4x, patchwork, profvis, fGarch, mvtnorm, beepr)
source("fun.R")
require(R.utils)

## To store MC data
beta_poly = c()
beta_obs_all = c()
coverage_poly_ori = c()
length_poly_ori = c()
## To store Theoretical data
coverage_poly_theory = c()
length_poly_theory = c()
## To store MC Randomization data
coverage_poly_ori_R= c()
length_poly_ori_R = c()

B = 200 ## samples
SampleNormal = function(beta, X, sigma) X%*%beta + rnorm(dim(X)[1], sd = sigma)
n = 100; p = 3; rho = 0.1
sigma = 1; f=0.5; level=0.9; alpha=0.1
beta = c(1,0.5,0.2)
estimateVar=FALSE
iters = 300  ## number of coverages to record

Sigma_X = matrix(nrow=p, ncol=p)
for (i in 1:p){
  for (j in 1:p){
    Sigma_X[i, j] = rho^(abs(i - j))
  }
}

pb = tkProgressBar(title = "Tk progress bar",      # Window title
                   label = "Percentage completed", # Window label
                   min = 0,      # Minimum value of the bar
                   max = B, # Maximum value of the bar
                   initial = 0,  # Initial value of the bar
                   width = 300)  # Width of the window

for (i in 1:iters){
  set.seed(i+6000)   # different seed
  X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
  y = SampleNormal(beta, X, sigma = sigma)    ## Actual sigma = 5
  
  ## Estimated sigma
  est_sig = sigma_hat(y, X)
  
  ############# MC ##########################
  ## Carry out LASSO to determine sub-model
  initial = Selection(y, X, type="lasso.fixed.lambda", p = p, first=TRUE, signs=FALSE)      ## returns lambda, coefficient vector, and variables
  vars_initial = initial$vars
  vars_m = which(vars_initial)     ## selected variables' number for the sub-model
  X_m = X[,vars_m, drop=FALSE]
  if(length(vars_m) == 0){ next }
  
  # Actual values of sub-model coefficients
  mu = X %*% beta
  eta = X_m %*% (solve(t(X_m) %*% X_m))
  beta_obs = t(eta) %*% mu
  est_sigma = sigma_hat(y, X)   ## estimate variance case

  
  for (j in 1:length(vars_m)){
  
    cat(sprintf("Inference for variable %i in iteration %i \n",vars_m[j], i))
    ## Calculate residual
    eta_j = eta[,j, drop=FALSE]   ## want to test the variable j in X_m
    denom = sum(t(eta_j) %*% eta_j)
    resid = (diag(n) - ((eta_j/denom) %*% t(eta_j))) %*% y  ## residual: projection to orthogonal of eta (dim: 200 x 1)
    
    ## Actual projection (partial regression) parameter value
    beta_obs_j = beta_obs[j]
    
    ## target parameter estimate
    betahat_obs = t(eta_j) %*% y
    
    Sigma2 = sqrt(sum((sigma^2)*(t(eta_j) %*% eta_j)))   ## variance of estimate
    
    pivot = function(param, lower_bound=FALSE, randomization=FALSE){
      coverages = numeric(B)
      for(b in 1:B){
        condition = FALSE
        counter = 1
        while(condition == FALSE){
          ## Generate beta_hat distribution
          betahat_sim = rnorm(1, mean = param, sd = Sigma2)
          y_sim = resid + ((eta_j %*% betahat_sim) / denom)
          if(randomization){
            y_sim = y_sim + rnorm(n, sd=2.5*est_sig$sigmahat)
          }
          selected_vars = Selection(y_sim, X, type="lasso.fixed.lambda", p = p, first=FALSE, lam=initial$lambda, signs = FALSE)   ## Use the same lambda since it is a function of X and errors only
          if(all(selected_vars == vars_initial) & !is.na(selected_vars[vars_m[j]])){
            condition=TRUE
            mes = "OK"

          }
          
          counter = counter + 1
          if((counter > 100)){
            print("here")
            condition=TRUE
            mes = "nothing"
          }
        }
        
        if(lower_bound & (mes=="OK")){
          coverages[b] = betahat_sim <= betahat_obs
        }else if(lower_bound & (mes == "nothing")){
          coverages[b] = NaN
        }
        
        if(!lower_bound & (mes=="OK")){
          coverages[b] = betahat_sim > betahat_obs
        }else if(!lower_bound & (mes == "nothing")){
          # print("here upper")
          coverages[b] = NaN
        }
        pctg <- paste(round(b/B *100, 0), "% completed")
        setTkProgressBar(pb, b, label = pctg)
      }
      return(mean(coverages, na.rm=TRUE)) ## cdf or pvalues
    }
    
    L = betahat_obs
    params_grid = c(L - 10*Sigma2, L + 10*Sigma2)  ## best estimate of the bounds
    
    if(sign(params_grid[1]) == -1){
      params_grid[1] = params_grid[1]
    }else{
      params_grid[1] = -params_grid[1]
    }
    
    if(sign(params_grid[2]) == 1){
      params_grid[2] = params_grid[2]
    }else{
      params_grid[2] = -params_grid[2]
    }
    
    # cat("Parameters: ",params_grid)
    
    CI_bisect = function(pivot,alpha,prec=1e-2, lower_bound=FALSE, randomization=FALSE){
      ## try initial value
      trial = mean(params_grid)
      if(lower_bound){
        cdf = pivot(trial, lower_bound, randomization=randomization)
        gtrial = cdf-(1-alpha/2)
        ## if NaN is returned, we assume that the CDF = 0
        if(is.nan(cdf)){ gtrial = 0 - (1-alpha/2) }
      }else{
        pvalue = pivot(trial, randomization=randomization)
        gtrial = pvalue-(1-alpha/2)
        ## if NaN is returned, we assume that the pvalue = 0
        if(is.nan(pvalue)){ gtrial = 0 - (1-alpha/2) }
      }
      
      if(gtrial == 0) return(gtrial)
      
      condition = ifelse(lower_bound, gtrial < 0, gtrial > 0)
      if(condition){
        
        ## Left boundary
        a = params_grid[1]
        
        if(lower_bound){
          cdf = pivot(a, lower_bound, randomization=randomization)
          ga = cdf-(1-alpha/2)
          if(is.nan(cdf)){ ga = 0 - (1-alpha/2) }
        }else{
          pvalue = pivot(a, randomization=randomization)
          ga = pvalue-(1-alpha/2)
          ## if NaN is returned, we assume that the pvalue = 0
          if(is.nan(pvalue)){ ga = 0 - (1-alpha/2) }
        }
        
        condition = ifelse(lower_bound, ga<0, ga >0)
        while(condition){ ## check if boundary works
          a <- a*2 ## decrease boundary if it does not work yet
          if(lower_bound){
            cdf = pivot(a, lower_bound, randomization=randomization)
            ga <- cdf-(1-alpha/2)
            if(is.nan(cdf)){ ga = 0 - (1-alpha/2) }
            if(ga > 0){ condition=FALSE}
          }else{
            pvalue = pivot(a, randomization=randomization)
            ga <- pvalue-(1- alpha/2)
            if(is.nan(pvalue)){ ga = 0 - (1-alpha/2) }
            if(ga < 0){ condition=FALSE}
          }
          if (a < -1e4){
            print("Error: g does not converge to 0 as x -> -Inf")
            return(-Inf)
          }
        }
        
        ## Right boundary
        b = trial
        gb = gtrial
        
      }else{
        
        a = trial
        ga = gtrial
        
        b = params_grid[2]
        
        if(lower_bound){
          cdf = pivot(b, lower_bound, randomization=randomization)
          gb = cdf-(1-alpha/2)
          if(is.nan(cdf)){ gb = 0 - (1-alpha/2) }
        }else{
          pvalue = pivot(b, randomization=randomization)
          gb = pvalue-(1-alpha/2)
          ## if NaN is returned, we assume that the pvalue = 0
          if(is.nan(pvalue)){ gb = 0 - (1-alpha/2) }
        }
        
        condition = ifelse(lower_bound, gb>0, gb<0)
        
        while(condition){ ## check if boundary works
          b <- b*2                  ## increase boundary if it does not work yet
          if(lower_bound){
            cdf = pivot(b, lower_bound, randomization=randomization)
            gb <- cdf-(1- alpha/2)
            if(is.nan(cdf)){ gb = 0 - (1-alpha/2) }
            if(gb<0){ condition = FALSE}
          }else{
            # print("hi6")
            pvalue = pivot(b, randomization=randomization)
            gb <- pvalue-(1- alpha/2)
            if(is.nan(pvalue)){ gb = 0 - (1-alpha/2) }
            if(gb>0){ condition = FALSE}
          }
          
          if (b > 1e4){
            print("Error: g does not converge to 0 as x -> Inf")
            return(Inf)
          }
        }
      }

      ####################  start bisectioning ###############
      count = 1
      while(abs(b-a)>prec | count == 20){
        c = (a+b)/2
        if(lower_bound){
          cdf = pivot(c, lower_bound=lower_bound, randomization=randomization)
          gc = cdf - (1-alpha/2)
          ## if NaN is returned, we assume that the CDF = 0
          if(is.nan(cdf)){ gc = 0 - (1-alpha/2) }
        }else{
          pvalue = pivot(c , randomization=randomization)
          gc = pvalue - (1-alpha/2)
          if(is.nan(pvalue)){ gc = 0 - (1-alpha/2) }
        }
        ## If exact convergence return the root, c
        if (gc==0) return(c)
        
        if (gc*ga<0){
          b = c
          gb = gc
        }else{
          a = c
          ga = gc
        }
        cat(sprintf("Value of a = %.3f at iteration %s \n",a,count))
        cat(sprintf("Value of b = %.3f at iteration %s \n",b,count))
        cat(sprintf("Error = %.3f at iteration %s \n",(b-a), count))
        count = count + 1
      }
      return((a+b)/2)
    }
    print("LOWER")
    lower = CI_bisect(pivot, alpha=alpha, lower_bound=TRUE)    
    print("UPPER")
    upper = CI_bisect(pivot, alpha=alpha, lower_bound=FALSE)   
    lower_R = CI_bisect(pivot, alpha=alpha, lower_bound=TRUE, randomization=TRUE)   
    upper_R = CI_bisect(pivot, alpha=alpha, lower_bound=FALSE, randomization=TRUE)
    
    coverage_poly_ori = c(coverage_poly_ori, (lower < beta_obs_j)*(upper > beta_obs_j))    ## check whether greater than the parameter value or not: if 0 => do not cover
    coverage_poly_ori_R = c(coverage_poly_ori_R, (lower_R < beta_obs_j)*(upper_R > beta_obs_j))    ## check whether greater than the parameter value or not: if 0 => do not cover
    
    length_poly_ori = c(length_poly_ori, upper - lower)   ## calculate the length
    length_poly_ori_R = c(length_poly_ori_R, upper_R - lower_R)   ## calculate the length
  }
  
  ######################## Theoretical comparison ##################
    ### actual polyhedral lemma
    if(estimateVar){
      theory_res = fixedLassoInf(X, y, beta = initial$coefs, lambda = initial$lam, sigma = est_sigma$sigmahat,intercept = FALSE, type="partial")
    }else{
      theory_res = fixedLassoInf(X, y, beta = initial$coefs, lambda = initial$lam, sigma = sigma, intercept = FALSE, type="partial")
    }

  ## Extract the CIs from theoretical
  CI = theory_res$ci
  
  ## Store coverages yes or no for the iteration i
  coverage_poly_theory = c(coverage_poly_theory, (CI[, 1] < beta_obs) * (CI[, 2] > beta_obs))
  ## Store Length for the iteration i
  length_poly_theory = c(length_poly_theory, CI[, 2] - CI[, 1])
  ## Make CI for the iteration i
  beta_poly = c(beta_poly, abs(beta[vars_m]))
}

my_results=list(beta_poly=beta_poly,
          coverage_poly_theory=coverage_poly_theory,
          coverage_poly_ori=coverage_poly_ori,
          coverage_poly_ori_R=coverage_poly_ori_R,
          length_poly_theory=length_poly_theory,
          length_poly_ori=length_poly_ori,
          length_poly_ori_R=length_poly_ori_R
)

# save(my_results, file= "0901_length10502.Rdata")
