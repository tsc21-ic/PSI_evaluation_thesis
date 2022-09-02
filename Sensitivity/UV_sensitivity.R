pacman::p_load(selectiveInference, prospectr, knockoff, mvtnorm, glmnet, stabs, MASS, dplyr,ggthemes, kableExtra, ggplot2, tibble, tidyverse,
               viridis, data.table, randomcoloR, cowplot, ggpubr, ggh4x, patchwork, profvis, fGarch)
source("fun.R")

################## Objective: compare with data splitting ######################
SampleNormal = function(beta, X, sigma) X%*%beta + rnorm(dim(X)[1], sd = sigma)
SampleT = function(beta, X, df) X%*%beta + rt(dim(X)[1], df = df)
SampleSkewNormal = function(beta, X, xi) X%*%beta + rsnorm(dim(X)[1], mean = 0, sd = 1, xi = xi)  ## here, we fix sd = 1, mainly look at how skewness affects performance

## Function
COV2 = function(n, p, rho, B, SelType, f, level, SampleType, sigma=NULL, df=NULL, xi = NULL,estimateVar = FALSE,coefficient=FALSE, projection=FALSE){
  # n;p;rho=0;B=2; SelType; f=0.5; level;sigma=5; estimateVar = TRUE; coefficient=TRUE; projection=FALSE
  n1 = floor(f*n)
  n2 = n - n1
  gamma = sqrt(1/f - 1)
  k = sqrt((1 + gamma^(-2)))
  
  Sigma_X = matrix(nrow = p, ncol = p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma_X[i, j] = rho^(abs(i - j))
    }
  }
  
  beta0_DS = c()
  beta0_R = c()
  COV_DS_HD = c()
  COV_R_HD = c()
  LENGTH_DS = c()
  LENGTH_R = c()
  
  pb = txtProgressBar(min = 0, max = B, style = 3)
  for (i in 1:B){
    
    ########################### DATA ###########################
    set.seed(i)
    
    beta0 = c(1, -1, 0.5, -0.5, 0.2, -0.2, rep(0, p - 6))    ## First 6 betas are non-zero
    
    X = rmvnorm(n, mean = rep(0, p), sigma = Sigma_X)
    
    if(SampleType == "normal"){
      y = SampleNormal(beta0, X, sigma = sigma)    ## Actual sigma = 5
    }else if(SampleType == "t"){
      y = SampleT(beta0, X, df = df)
    }else if(SampleType == "sn"){
      y = SampleSkewNormal(beta0, X, xi = xi)
    }
    
    ## Estimated sigma
    est_sig = sigma_hat(y, X)
    
    ##################### Data splitting ########################
    ind = I(X,n,n1)  ## indices for n1
    X1 = X[ind, ]
    y1 = y[ind]
    S_DS = Selection(y1, X1, type = SelType, p = p)  ## TRUE => selected variable , FALSE => not selected variable
    
    ## Selected coefficients inference
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
        beta_DS_HD = solve(t(X_inf)%*%X_inf)%*%t(X_inf)%*%X[-ind, ]%*%beta0    ## actual partial target coefficient upon model selection beta_M
        # beta_DS_HD_hat = solve(t(X_inf)%*%X_inf)%*%t(X_inf)%*%y_inf    ## beta^hat_M
        model = lm(y_inf ~ X_inf - 1)
        if(estimateVar){
          CIs_DS_HD = confint.truesigma(model, level=level, truesigma = est_sig$sigmahat, df=est_sig$df, estimate = TRUE)
        }else{
          CIs_DS_HD = confint.truesigma(model, level=level, truesigma = sigma)
        }
      }
    }
    
    #################### Randomisation ####################
    ## "Hold-out" (Correct) approach
    
    ## Selection stage
    if(estimateVar){
      w = rnorm(n, sd = est_sig$sigmahat)
    }else{
      w = rnorm(n, sd = sigma)
    }
    u = y + gamma*w
    S_R = Selection(u, X, type = SelType, p = p)  
    
    ## Selected coefficients inference
    if(coefficient){
      X_inf = X
      v = y - w/gamma
      model = lm(v ~ X_inf - 1)
      if(estimateVar){
        CIs_R_HD = confint.truesigma(model, level = level, truesigma = k*est_sig$sigmahat, df=est_sig$df, estimate = TRUE)[S_R, ] 
      }else{
        CIs_R_HD = confint.truesigma(model, level = level, truesigma = k*sigma)[S_R, ]   ## randomization confidence intervals
      }
    }
    
    ## Projection coefficients inference
    if(projection){
      if(sum(S_R) != 0){
        X_inf = X[ , S_R]
        v = y - w/gamma
        beta_R_HD = solve(t(X_inf)%*%X_inf)%*%t(X_inf)%*%X%*%beta0   ## actual partial target coefficient upon model selection beta_M
        model = lm(v ~ X_inf - 1)
        if(estimateVar){
          CIs_R_HD = confint.truesigma(model, level= level, truesigma = k*est_sig$sigmahat, df=est_sig$df, estimate = TRUE)
        }else{
          CIs_R_HD = confint.truesigma(model, level = level, truesigma = k*sigma)   ## randomization confidence intervals  
        }
      }
    }
    
    ########################################################################
    beta0_DS = c(beta0_DS, abs(beta0[S_DS]))     ## absolute value of betas selected by data splitting
    beta0_R = c(beta0_R, abs(beta0[S_R]))        ## absolute value of betas selected by randomization
    ########################################################################
    ## RESULTS
    
    ## Data Splitting
    if (sum(S_DS) == 1){
      if(coefficient){ COV_DS_HD = c(COV_DS_HD, (CIs_DS_HD[1] < beta0[S_DS])*(CIs_DS_HD[2] > beta0[S_DS])) }
      if(projection){ COV_DS_HD = c(COV_DS_HD, (CIs_DS_HD[1] < beta_DS_HD)*(CIs_DS_HD[2] > beta_DS_HD)) }
      LENGTH_DS = c(LENGTH_DS, CIs_DS_HD[2] - CIs_DS_HD[1])    ## length of confidence interval
    }
    
    if (sum(S_DS) > 1){
      if(coefficient){ COV_DS_HD = c(COV_DS_HD, (CIs_DS_HD[ , 1] < beta0[S_DS])*(CIs_DS_HD[ , 2] > beta0[S_DS])) }
      if(projection){ COV_DS_HD = c(COV_DS_HD, (CIs_DS_HD[ , 1] < beta_DS_HD)*(CIs_DS_HD[ , 2] > beta_DS_HD)) }
      LENGTH_DS = c(LENGTH_DS, CIs_DS_HD[ , 2] - CIs_DS_HD[ , 1])
    }
    
    ## Randomization
    if (sum(S_R) == 1){
      if(coefficient){ COV_R_HD = c(COV_R_HD, (CIs_R_HD[1] < beta0[S_R])*(CIs_R_HD[2] > beta0[S_R])) }
      if(projection){ COV_R_HD = c(COV_R_HD, (CIs_R_HD[1] < beta_R_HD)*(CIs_R_HD[2] > beta_R_HD)) }
      LENGTH_R = c(LENGTH_R,  CIs_R_HD[2] - CIs_R_HD[1]) 
    }
    
    if (sum(S_R) > 1){
      if(coefficient){ COV_R_HD = c(COV_R_HD, (CIs_R_HD[ , 1] < beta0[S_R])*(CIs_R_HD[ , 2] > beta0[S_R])) }
      if(projection){ COV_R_HD = c(COV_R_HD, (CIs_R_HD[ , 1] < beta_R_HD)*(CIs_R_HD[ , 2] > beta_R_HD)) }
      LENGTH_R = c(LENGTH_R,  CIs_R_HD[ , 2] - CIs_R_HD[ , 1])
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(list(beta0_DS=beta0_DS, beta0_R=beta0_R, COV_DS_HD=COV_DS_HD, COV_R_HD=COV_R_HD, LENGTH_DS=LENGTH_DS, LENGTH_R=LENGTH_R))
}

############### To run results ###########################################
run_result_cov2 = function(n, p, rhos, B, SelType, fs, level, SampleType, sigma=NULL, df=NULL, xi = NULL, estimateVar = FALSE, coefficient=FALSE, projection=FALSE, save=FALSE, load=FALSE){
  count = 1
  for(j in rhos){
    for(k in fs){
      
      ############################# Assign object  ###########################
      ## Coefficient
      if(coefficient){
        
        if(SampleType == "t"){
          result_obj = paste0("cov", count, '_df_', df)  ## e.g. cov1_df_1
        }
        if(SampleType == "sn"){
          result_obj = paste0("cov", count, '_sn_', xi)  ## e.g. cov1_sn_1
        }
        
        if(SampleType == "normal"){
          if(estimateVar){
            result_obj = paste0("cov", count, '_est_sigma_', sigma)  ## e.g. cov1_est_sigma_1
          }else{
            result_obj = paste0("cov", count, '_fixed_sigma_', sigma)  ## e.g. cov1_fixed_sigma_1
          }
        }
        
      }
      
      ## Projection
      if(projection){
        
        if(SampleType == "t"){
          result_obj = paste0("cov_pj", count, '_df_', df)  ## e.g. cov_pj1_df_1
        }
        if(SampleType == "sn"){
          result_obj = paste0("cov_pj", count, '_sn_', xi)  ## e.g. cov_pj1_sn_1
        }
        
        if(SampleType == "normal"){
          if(estimateVar){
            result_obj = paste0("cov_pj", count, '_est_sigma_', sigma)  ## e.g. cov_pj1_est_sigma_1
          }else{
            result_obj = paste0("cov_pj", count, '_fixed_sigma_', sigma)  ## e.g. cov_pj1_fixed_sigma_1
          }
        }
        
      }
      
      ############################### Assign object ###############################
      if(!save & !load){
        assign(result_obj, COV2(n, p = p, rho = j, B, SelType, f = k, 
                                level = level, SampleType = SampleType, sigma=sigma, df=df, xi = xi,
                                estimateVar = estimateVar,
                                coefficient = coefficient, projection=projection), envir = parent.frame())
      }
      
      ############################### Save object ###############################
      if(save){
        if(estimateVar){
          filename = paste0(result_obj, '.Rdata')
        }else{
          filename = paste0(result_obj, '.Rdata')
        }
        save(list=result_obj, file=filename)
      }
      
      ############################### Load object ###############################
      if(load){
        if(estimateVar){
          filename = paste0(result_obj, '.Rdata')
        }else{
          filename = paste0(result_obj, '.Rdata')
        }
        load(file=filename, verbose=TRUE, envir = .GlobalEnv)
      }
      count = count + 1
    }
  }
}


##### Run results #### 
#################### t-distributed errors ##############################
dfs = c(1,2,5,10,Inf)   ## Inf correspond to normal case
n = 200; p = 30; rhos = c(0, 0.5) ; SelType = "lasso"; fs = c(1/2, 3/4); level = 0.9

## Selected coefficient
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[1], estimateVar = TRUE, coefficient = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[2], estimateVar = TRUE, coefficient = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[3], estimateVar = TRUE, coefficient = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[4], estimateVar = TRUE, coefficient = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[5], estimateVar = TRUE, coefficient = TRUE)

## Save results
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[1], estimateVar = TRUE, coefficient = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[2], estimateVar = TRUE, coefficient = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[3], estimateVar = TRUE, coefficient = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[4], estimateVar = TRUE, coefficient = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[5], estimateVar = TRUE, coefficient = TRUE, save=TRUE)


## Projection
n = 100; p = 150; rhos = c(0, 0.5) ; SelType = "stabs.lasso"; fs = c(1/2, 3/4); level = 0.9

run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[1], estimateVar = TRUE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[2], estimateVar = TRUE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[3], estimateVar = TRUE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[4], estimateVar = TRUE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[5], estimateVar = TRUE, projection = TRUE)

## Save results
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[1], estimateVar = TRUE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[2], estimateVar = TRUE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[3], estimateVar = TRUE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[4], estimateVar = TRUE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "t", df = dfs[5], estimateVar = TRUE, projection = TRUE, save=TRUE)

#################### Skewed normal distribution ####################
xis = c(0.1, 0.7, 1, 3, 10)   ## 1 correspond to normal case
n = 200; p = 30; rhos = c(0, 0.5) ; SelType = "lasso"; fs = c(1/2, 3/4); level = 0.9

## Selected coefficient
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[1], estimateVar = TRUE, coefficient = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[2], estimateVar = TRUE, coefficient = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[3], estimateVar = TRUE, coefficient = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[4], estimateVar = TRUE, coefficient = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[5], estimateVar = TRUE, coefficient = TRUE)

## Save results
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[1], estimateVar = TRUE, coefficient = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[2], estimateVar = TRUE, coefficient = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[3], estimateVar = TRUE, coefficient = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[4], estimateVar = TRUE, coefficient = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[5], estimateVar = TRUE, coefficient = TRUE, save=TRUE)

## Projection
n = 100; p = 150; rhos = c(0, 0.5) ; SelType = "stabs.lasso"; fs = c(1/2, 3/4); level = 0.9

run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[1], estimateVar = TRUE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[2], estimateVar = TRUE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[3], estimateVar = TRUE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[4], estimateVar = TRUE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[5], estimateVar = TRUE, projection = TRUE)

## Save results
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[1], estimateVar = TRUE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[2], estimateVar = TRUE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[3], estimateVar = TRUE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[4], estimateVar = TRUE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "sn", xi = xis[5], estimateVar = TRUE, projection = TRUE, save=TRUE)


#################### Variance ####################
## Coefficient
# true sigma = 5, true sigma = 1 (we know (U, V) is dependent on the sigma value)
n = 200; p = 30; rhos = c(0, 0.5) ; SelType = "lasso"; fs = c(1/2, 3/4); level = 0.9
## sigma = 5
run_result_cov2(n,p,rhos, B = 2000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = FALSE, coefficient=TRUE)
run_result_cov2(n,p,rhos, B = 2000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = TRUE, coefficient=TRUE)

## Save results
run_result_cov2(n,p,rhos, B = 2000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = FALSE, coefficient=TRUE, save=TRUE)
run_result_cov2(n,p,rhos, B = 2000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = TRUE, coefficient=TRUE, save=TRUE)

## sigma = 1
n = 200; p = 30; rhos = c(0, 0.5) ; SelType = "lasso"; fs = c(1/2, 3/4); level = 0.9
run_result_cov2(n,p,rhos, B = 2000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = FALSE, coefficient=TRUE)
run_result_cov2(n,p,rhos, B = 2000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = TRUE, coefficient=TRUE)

## Save results
run_result_cov2(n,p,rhos, B = 2000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = FALSE, coefficient=TRUE, save=TRUE)
run_result_cov2(n,p,rhos, B = 2000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = TRUE, coefficient=TRUE, save=TRUE)


## Projection
## n = 100, p = 150 (Need to use projection method and LASSO stability selection: less prone to wrong estimates, and identifiability issues)
## Run higher iterations because harder to select active variables now
n = 100; p = 150; rhos = c(0, 0.5) ; SelType = "stabs.lasso"; fs = c(1/2, 3/4); level = 0.9
## sigma = 5
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = FALSE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = TRUE, projection = TRUE)

run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = FALSE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = TRUE, projection = TRUE, save=TRUE)

## sigma = 1
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = FALSE, projection = TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = TRUE, projection = TRUE)

run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = FALSE, projection = TRUE, save=TRUE)
run_result_cov2(n,p,rhos,B = 2000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = TRUE, projection = TRUE, save=TRUE)
