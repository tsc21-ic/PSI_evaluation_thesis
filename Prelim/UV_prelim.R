pacman::p_load(selectiveInference, prospectr, knockoff, mvtnorm, glmnet, stabs, MASS, dplyr,ggthemes, kableExtra, ggplot2, tibble, tidyverse,
               viridis, pbapply, data.table, progress, randomcoloR, cowplot, ggpubr,tcltk, ggh4x, patchwork, profvis, fGarch, mvtnorm, beepr)
source('fun.R')

# Generative model (t-distribution)
SampleNormal = function(beta, X, sigma) X%*%beta + rnorm(dim(X)[1], sd = sigma)

B = 100  ## samples
n = 100; p = 3; rho = 0.1
sigma = 1; f=0.5; level=0.9; alpha=0.1
beta = c(1,0.8,0.2)
estimateVar=FALSE
coefficient=FALSE
projection=TRUE
iters = 200   ## number of coverages to record


pb = tkProgressBar(title = "Tk progress bar",      # Window title
                   label = "Percentage completed", # Window label
                   min = 0,      # Minimum value of the bar
                   max = iters, # Maximum value of the bar
                   initial = 0,  # Initial value of the bar
                   width = 300)  # Width of the window

## For data splitting
n1 = floor(f*n)
n2 = n - n1

## For randomization
gamma = sqrt(1/f - 1)
k = sqrt((1 + gamma^(-2)))

Sigma_X = matrix(nrow=p, ncol=p)
for (i in 1:p){
  for (j in 1:p){
    Sigma_X[i, j] = rho^(abs(i - j))
  }
}

## To store R data
COV_R_HD = c()
beta_R = c()
LENGTH_R = c()

## To store DS data
COV_DS_HD = c()
beta_DS = c()
LENGTH_DS = c()

initial_all=c()

for (i in 1:iters){
  set.seed(i+2000)
  pctg <- paste(round(i/iters *100, 0), "% completed")
  setTkProgressBar(pb, iters, label = pctg)
  
  X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
  y = SampleNormal(beta, X, sigma = sigma)    ## Actual sigma = 5
  
  ## Estimated sigma
  est_sig = sigma_hat(y, X)
  
  ## Carry out LASSO to determine sub-model
  initial = Selection(y, X, type="lasso.fixed.lambda", p = p, first=TRUE, signs=FALSE)      ## returns lambda, coefficient vector, and variables
  initial_all = c(initial_all , initial)
  ###################### Data splitting ########################
  ind = I(X,n,n1)  ## indices for n1
  X1 = X[ind, ]
  y1 = y[ind]
  # S_DS = Selection(y1, X1, type = SelType, p = p)  ## TRUE => selected variable , FALSE => not selected variable
  S_DS = Selection(y1, X1, type="lasso.fixed.lambda", p = p, first=TRUE, signs=FALSE)$vars
  
  ## "Hold-out" (Correct) approach
  if(coefficient){
    X_inf = X[-ind , ]
    y_inf = y[-ind]
    model = lm(y_inf ~ X_inf - 1)
    if(estimateVar){
      CIs_DS_HD = confint.truesigma(model, level = level, truesigma = est_sig$sigmahat, df=est_sig$df, estimate = TRUE)[S_DS, ]
    }else{
      CIs_DS_HD = confint.truesigma(model, level = level, truesigma = sigma)[S_DS, ]
    }
  }
  
  ## Projection parameter inference
  if(projection){
    if(sum(S_DS) != 0){
      X_inf = X[-ind, S_DS]   ## s is data dependent
      y_inf = y[-ind]
      beta_DS_HD = solve(t(X_inf)%*%X_inf)%*%t(X_inf)%*%X[-ind, ]%*%beta    ## actual partial target coefficient upon model selection beta_M
      model = lm(y_inf ~ X_inf - 1)
      if(estimateVar){
        CIs_DS_HD = confint.truesigma(model, level=level, truesigma = est_sig$sigmahat, df=est_sig$df, estimate = TRUE)
      }else{
        CIs_DS_HD = confint.truesigma(model, level=level, truesigma = sigma)  ## true sigma
      }
    }
  }
  
  #################### Randomisation ####################

  ## Selection
  if(estimateVar){
    w = rnorm(n, sd = est_sig$sigmahat)
  }else{
    w = rnorm(n, sd = sigma)
  }
  u = y + gamma*w
  S_R = Selection(u, X, type="lasso.fixed.lambda", p = p, lam = initial$lambda) ## Try using the same lambda by negahban
  
  ## Inference
  if(coefficient){
    X_inf = X
    v = y - w/gamma
    model = lm(v ~ X_inf - 1)
    if(estimateVar){
      CIs_R_HD = confint.truesigma(model, level = level, truesigma = k*est_sig$sigmahat, df=est_sig$df, estimate = TRUE)[S_R, ]
    }else{
      CIs_R_HD = confint.truesigma(model, level = level, truesigma = k*sigma)[S_R, ]   ## randomisation confidence intervals
    }
  }
  
  ## Projection inference
  if(projection){
    if (sum(S_R) != 0){
      X_inf = X[ , S_R]
      v = y - w/gamma
      beta_R_HD = solve(t(X_inf)%*%X_inf)%*%t(X_inf)%*%X%*%beta
      model = lm(v ~ X_inf - 1)
      if(estimateVar){
        CIs_R_HD = confint.truesigma(model, level= level, truesigma = k*est_sig$sigmahat, df=est_sig$df, estimate = TRUE)
      }else{
        CIs_R_HD = confint.truesigma(model, level = level, truesigma = k*sigma)   ## randomisation confidence intervals
      }
    }
  }
  
  beta_DS = c(beta_DS, abs(beta[S_DS]))
  beta_R = c(beta_R, abs(beta[S_R]))
  
  ################### RESULTS ##########################################
  
  
  ## Data Splitting
  if (sum(S_DS) == 1){
    # COV_DS_FV = c(COV_DS_FV, (CIs_DS_FV[1] < beta[S_DS])*(CIs_DS_FV[2] > beta[S_DS]))   ## TRUE => CIs covers selected variables' coefficient, FALSE => CIs does not cover
    if(coefficient){ COV_DS_HD = c(COV_DS_HD, (CIs_DS_HD[1] < beta[S_DS])*(CIs_DS_HD[2] > beta[S_DS])) }
    if(projection){ COV_DS_HD = c(COV_DS_HD, (CIs_DS_HD[1] < beta_DS_HD)*(CIs_DS_HD[2] > beta_DS_HD)) }
    LENGTH_DS = c(LENGTH_DS, CIs_DS_HD[2] - CIs_DS_HD[1])    ## length of confidence interval
  }
  
  if (sum(S_DS) > 1){
    # COV_DS_FV = c(COV_DS_FV, (CIs_DS_FV[ , 1] < beta[S_DS])*(CIs_DS_FV[ , 2] > beta[S_DS]))
    if(coefficient){ COV_DS_HD = c(COV_DS_HD, (CIs_DS_HD[ , 1] < beta[S_DS])*(CIs_DS_HD[ , 2] > beta[S_DS])) }
    if(projection){ COV_DS_HD = c(COV_DS_HD, (CIs_DS_HD[ , 1] < beta_DS_HD)*(CIs_DS_HD[ , 2] > beta_DS_HD)) }
    LENGTH_DS = c(LENGTH_DS, CIs_DS_HD[ , 2] - CIs_DS_HD[ , 1])
  }
  
  ## Randomisation
  if (sum(S_R) == 1){
    if(coefficient){ COV_R_HD = c(COV_R_HD, (CIs_R_HD[1] < beta[S_R])*(CIs_R_HD[2] > beta[S_R])) }
    if(projection){ COV_R_HD = c(COV_R_HD, (CIs_R_HD[1] < beta_R_HD)*(CIs_R_HD[2] > beta_R_HD)) }
    LENGTH_R = c(LENGTH_R,  CIs_R_HD[2] - CIs_R_HD[1])
  }
  
  if (sum(S_R) > 1){
    if(coefficient){ COV_R_HD = c(COV_R_HD, (CIs_R_HD[ , 1] < beta[S_R])*(CIs_R_HD[ , 2] > beta[S_R])) }
    if(projection){ COV_R_HD = c(COV_R_HD, (CIs_R_HD[ , 1] < beta_R_HD)*(CIs_R_HD[ , 2] > beta_R_HD)) }
    LENGTH_R = c(LENGTH_R,  CIs_R_HD[ , 2] - CIs_R_HD[ , 1])
  }
}

my_results_DS_UV =list(beta_DS=beta_DS,
                       beta_R=beta_R,
                       COV_DS_HD=COV_DS_HD,
                       COV_R_HD=COV_R_HD,
                       LENGTH_DS=LENGTH_DS,
                       LENGTH_R=LENGTH_R)

# save(my_results_DS_UV, file="DS_UV_length_10802.Rdata")
