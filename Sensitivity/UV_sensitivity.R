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
run_result_cov2(n,p,rhos, B = 1000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = FALSE, coefficient=TRUE)
run_result_cov2(n,p,rhos, B = 1000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = TRUE, coefficient=TRUE)

## Save results
run_result_cov2(n,p,rhos, B = 1000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = FALSE, coefficient=TRUE, save=TRUE)
run_result_cov2(n,p,rhos, B = 1000, SelType, fs, level, SampleType = "normal", sigma=5, estimateVar = TRUE, coefficient=TRUE, save=TRUE)

## sigma = 1
n = 200; p = 30; rhos = c(0, 0.5) ; SelType = "lasso"; fs = c(1/2, 3/4); level = 0.9
run_result_cov2(n,p,rhos, B = 1000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = FALSE, coefficient=TRUE)
run_result_cov2(n,p,rhos, B = 1000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = TRUE, coefficient=TRUE)

## Save results
run_result_cov2(n,p,rhos, B = 1000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = FALSE, coefficient=TRUE, save=TRUE)
run_result_cov2(n,p,rhos, B = 1000, SelType, fs, level, SampleType = "normal", sigma=1, estimateVar = TRUE, coefficient=TRUE, save=TRUE)


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





###################### Alternatively: load prepared data files ###########################################
# ## t-distributed
# dir = "~/github/Part_1_data"
# dfs = c(1,2,5,10,Inf)   ## Inf correspond to normal case
# xis = c(0.1, 0.7, 1, 3, 10)   ## 1 correspond to normal case
# rhos = c(0, 0.5) ; fs = c(1/2, 3/4); level = 0.9
# 
# for(i in dfs){
#   count = 1
#   for(j in rhos){
#     for(k in fs){
#       load(verbose=TRUE, file = paste0(dir, "/cov", count, "_df" , "_", i, ".Rdata"))
#       load(verbose=TRUE, file = paste0(dir, "/cov_pj", count, "_df_", i, ".Rdata"))
#       count=count+1
#     }
#   }
# }
# ## skewed normal
# for(i in xis){
#   count = 1
#   for(j in rhos){
#     for(k in fs){
#       load(verbose=TRUE, file = paste0(dir, "/cov", count, "_sn_", i, ".Rdata"))
#       load(verbose=TRUE, file = paste0(dir, "/cov_pj", count, "_sn_", i, ".Rdata"))
#       count=count+1
#     }
#   }
# }
# ## sigma
# for(i in c("fixed", "est")){
#   count = 1
#   for(j in rhos){
#     for(k in fs){
#       load(verbose=TRUE, file = paste0(dir, "/cov", count, "_" ,i, "_sigma_", 1, ".Rdata"))
#       load(verbose=TRUE, file = paste0(dir, "/cov_pj", count, "_" ,i, "_sigma_", 1, ".Rdata"))
#       load(verbose=TRUE, file = paste0(dir, "/cov", count, "_" ,i, "_sigma_", 5, ".Rdata"))
#       load(verbose=TRUE, file = paste0(dir, "/cov_pj", count, "_" ,i, "_sigma_", 5, ".Rdata"))
# 
#       count=count+1
#     }
#   }
# }

# ##### There are total 112 lists of data
# data.frame(do.call("rbind", lapply(ls(), function(x) {
#   obj = get(x)
#   if (class(obj) == "list")
#     c(name =x)
# })))
# 

#####################  Generate full data frame of results ##############################################################################
generatedf = function(rhos, fs, SampleType, sigma=NULL, df=NULL, xi = NULL, estimateVar=FALSE, coefficient=FALSE, projection=FALSE){
  list_cov = list()
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
      
      result = eval(parse(text = result_obj))
      list_cov = append(list_cov, result)
      count = count + 1
    }
  }
  
  list_cov_final=with(stack(list_cov), split(values, ind))   ## stack all combination of outputs into one as summary
  return(list_cov_final)
}


############# t-results ############# 
## n > p case
cov_t1_df = generatedf(rhos, fs,SampleType = "t", df=1, estimateVar=TRUE, coefficient = TRUE)
cov_t2_df = generatedf(rhos, fs,SampleType = "t", df=2, estimateVar=TRUE, coefficient = TRUE)
cov_t5_df = generatedf(rhos, fs,SampleType = "t", df=5, estimateVar=TRUE, coefficient = TRUE)
cov_t10_df = generatedf(rhos, fs,SampleType = "t", df=10, estimateVar=TRUE, coefficient = TRUE)
cov_tInf_df = generatedf(rhos, fs,SampleType = "t", df=Inf, estimateVar=TRUE, coefficient = TRUE)
t_dfs = list(cov_t1_df, cov_t2_df, cov_t5_df, cov_t10_df, cov_tInf_df)


## n < p case
cov_pj_t1_df = generatedf(rhos, fs, SampleType = "t", df=1, estimateVar=TRUE, projection = TRUE)
cov_pj_t2_df = generatedf(rhos, fs, SampleType = "t", df=2, estimateVar=TRUE, projection = TRUE)
cov_pj_t5_df = generatedf(rhos, fs, SampleType = "t", df=5, estimateVar=TRUE, projection = TRUE)
cov_pj_t10_df = generatedf(rhos, fs, SampleType = "t", df=10, estimateVar=TRUE, projection = TRUE)
cov_pj_tInf_df = generatedf(rhos, fs, SampleType = "t", df=Inf, estimateVar=TRUE, projection = TRUE)
t_pj_dfs = list(cov_pj_t1_df, cov_pj_t2_df, cov_pj_t5_df, cov_pj_t10_df, cov_pj_tInf_df)


############# Skewed normal ############# 
## n > p case
cov_sn1_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[1], estimateVar=TRUE, coefficient = TRUE)
cov_sn2_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[2], estimateVar=TRUE, coefficient = TRUE)
cov_sn3_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[3], estimateVar=TRUE, coefficient = TRUE)
cov_sn4_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[4], estimateVar=TRUE, coefficient = TRUE)
cov_sn5_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[5], estimateVar=TRUE, coefficient = TRUE)
sn_dfs = list(cov_sn1_df, cov_sn2_df, cov_sn3_df, cov_sn4_df, cov_sn5_df)

## n < p case
cov_pj_sn1_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[1], estimateVar=TRUE, projection = TRUE)
cov_pj_sn2_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[2], estimateVar=TRUE, projection = TRUE)
cov_pj_sn3_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[3], estimateVar=TRUE, projection = TRUE)
cov_pj_sn4_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[4], estimateVar=TRUE, projection = TRUE)
cov_pj_sn5_df = generatedf(rhos, fs,SampleType = "sn", xi=xis[5], estimateVar=TRUE, projection = TRUE)
sn_pj_dfs = list(cov_pj_sn1_df, cov_pj_sn2_df, cov_pj_sn3_df, cov_pj_sn4_df, cov_pj_sn5_df)


############# sigma = 5 #############
cov_fixed_df = generatedf(rhos, fs,  SampleType = "normal",sigma=5, estimateVar=FALSE,coefficient = TRUE)
cov_est_df = generatedf(rhos, fs, SampleType = "normal", sigma=5, estimateVar=TRUE,coefficient = TRUE)
normal_low_d_dfs = list(cov_est_df, cov_fixed_df)  ## KEEP THIS ORDER

cov_fixed_df = generatedf(rhos, fs, SampleType = "normal", sigma=5, estimateVar=FALSE,projection=TRUE)
cov_est_df = generatedf(rhos, fs, SampleType = "normal", sigma=5, estimateVar=TRUE,projection=TRUE)
normal_high_d_dfs = list(cov_est_df, cov_fixed_df)  ## KEEP THIS ORDER

############# sigma = 1 #############
cov_fixed_sigma1_df = generatedf(rhos, fs,  SampleType = "normal",sigma=1, estimateVar=FALSE,coefficient = TRUE)
cov_est_sigma1_df = generatedf(rhos, fs, SampleType = "normal", sigma=1, estimateVar=TRUE,coefficient = TRUE)
normal_low_d_sigma1_dfs = list(cov_est_sigma1_df, cov_fixed_sigma1_df)  ## KEEP THIS ORDER

cov_fixed_sigma1_df = generatedf(rhos, fs, sigma=1, SampleType = "normal", estimateVar=FALSE,projection=TRUE)
cov_est_sigma1_df = generatedf(rhos, fs, sigma=1, SampleType = "normal", estimateVar=TRUE,projection=TRUE)
normal_high_d_sigma1_dfs = list(cov_est_sigma1_df, cov_fixed_sigma1_df)  ## KEEP THIS ORDER


sigma_status = c("est", "fixed")

############## Generate output function ####################
## Function to take a list (to plot boxplots for the lengths)
generateResult  = function(...){
  
  args = sys.call()
  all_args = paste0(lapply(args[-1], as.character))
  
  abs_beta = factor(c(0, 0.2, 0.5, 1))
  
  ## Combine results for all beta
  all_beta = function(input, split){
    df = data.frame(beta=factor(), length=numeric(), type=NULL)
    if(split == "R"){
      LENGTH = input$LENGTH_R
      beta0 = input$beta0_R
    }
    
    if(split == "DS"){
      LENGTH = input$LENGTH_DS
      beta0 = input$beta0_DS
    }
    
    count <<- count + 1
    for(i in abs_beta){
      if(grepl("t", all_args)){
        tmp_df = try({rbind.data.frame(df, data.frame(beta=i,length=LENGTH[which(beta0 == i)], type=dfs[count]))}, silent=TRUE)
        if('try-error' %in% class(tmp_df)){ next } else{df = tmp_df}
      }
      else if(grepl("sn", all_args)){
        tmp_df = try({rbind.data.frame(df, data.frame(beta=i,length=LENGTH[which(beta0 == i)], type=xis[count]))}, silent=TRUE)
        if('try-error' %in% class(tmp_df)){ next } else{df = tmp_df}
      }
      else{
        tmp_df = try({rbind.data.frame(df, data.frame(beta=i,length=LENGTH[which(beta0 == i)], type=sigma_status[count]))}, silent=TRUE)
        if('try-error' %in% class(tmp_df)){ next } else{df = tmp_df}
      }
    }
    return(df)
  }
  
  count = 0
  df_list_R = lapply(..., all_beta, split = "R")
  df_R = rbindlist(df_list_R) %>% mutate("split" = "(U, V)")
  
  count = 0
  df_list_DS = lapply(..., all_beta, split = "DS")
  df_DS = rbindlist(df_list_DS) %>% mutate("split" = "DS")
  
  
  ## For plotting
  df_all = rbind(df_R, df_DS) %>% mutate(split = factor(split,levels=c("DS", "(U, V)")))
  if(grepl("t", all_args)){
    df_all = df_all %>% mutate(type = factor(type, level=dfs))
  }
  else if(grepl("sn", all_args)){
    df_all = df_all %>% mutate(type = factor(type,levels=xis))
  }
  
  var_names = c(`0` = paste0("bgroup('|', beta[i], '|')", "== 0"),
                `0.2` = paste0("bgroup('|', beta[i], '|')", "== 0.2"),
                `0.5` = paste0("bgroup('|', beta[i], '|')", "== 0.5"),
                `1` = paste0("bgroup('|', beta[i], '|')", "== 1"))
  
  tmp = df_all %>% split(., .$beta) %>% map(function(x){
    ggplot(data = x , aes(x=type, y=length))+
      geom_boxplot(position=position_dodge(0.7), linetype = "dashed")+
      stat_boxplot(aes(ymin = ..lower.., ymax = ..upper.., fill=type),outlier.shape = 1) +
      stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.2)+
      stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width=0.2) +
      stat_summary(fun = "mean", geom = "point", shape=18, size=2, color="cornsilk2")+
      theme_bw() +
      facet_nested(. ~ beta + split, labeller = labeller(beta = as_labeller(var_names, label_parsed)), scales = "free_y")+
      {if(grepl("normal", all_args)){
        scale_fill_discrete(labels = c(bquote(Estimated~sigma^2), bquote(Fixed~sigma^2)))
      }else if(grepl("t", all_args)){
        scale_fill_discrete(labels = c(dfs[1:(length(dfs)-1)],bquote(infinity)))
      }
      }+
      {if(grepl("t", all_args)) coord_cartesian(ylim = c(0,3))} +
      theme(strip.text.x = element_text(margin = margin(b = 0, t = 0), face =  "bold.italic"),
            axis.title = element_blank(),
            legend.position = "none",
            panel.spacing=unit(0.1,"lines"),
            panel.spacing.y = unit(0.9, "lines"),
            strip.background = element_blank(),
            legend.spacing.x = unit(1.0, 'cm'),
            axis.ticks.x=element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank())+
      { if(grepl("t", all_args)){ labs(fill=expression(italic(k)))
      }else if(grepl("sn", all_args)){ labs(fill=expression(paste(xi, ' ', sep = ' ')))
      }else{ guides(fill=guide_legend(title=""))}
      }
  })
  
  out_p = wrap_plots(tmp, nrow = 1) + plot_layout(guides = "collect")
  
  ## For table presentation
  ## Randomization
  test = function(input, split){
    if(split == "R"){
      LENGTH = input$LENGTH_R
      beta0 = input$beta0_R
      coverage = input$COV_R_HD
    }
    
    if(split == "DS"){
      LENGTH = input$LENGTH_DS
      beta0 = input$beta0_DS
      coverage = input$COV_DS_HD
    }
    
    table = data.frame(beta=beta0, length=LENGTH, coverage=coverage)
    
    ## For normal case
    args2 = sys.call()
    all_args2 = paste0(lapply(args2[-1], as.character))
    if(grepl("fixed", all_args2[1])){
      title_normal = "fixed"
    }else{
      title_normal = "est"
    }
    
    
    count <<- count + 1
    if(grepl("t", all_args)){
      table = table %>% group_by(beta) %>% 
        summarise(!! paste("length.df",dfs[count], sep = ".") := format(round(median(length),2), nsmall=2),
                  !! paste("coverage.df",dfs[count], sep = ".") := format(round(mean(coverage)*100, 2), nsmall = 1))
    }else if(grepl("sn", all_args)){
      table = table %>% group_by(beta) %>% 
        summarise(!! paste("length.xi",xis[count], sep = ".") := format(round(median(length),2), nsmall=2),
                  !! paste("coverage.xi",xis[count], sep = ".") := format(round(mean(coverage)*100, 2), nsmall = 1))
    }else{
      table = table %>% group_by(beta) %>% 
        summarise(!! paste("length.", title_normal, sep = ".") := format(round(median(length),2), nsmall=2),
                  !! paste("coverage.", title_normal, sep = ".") := format(round(mean(coverage)*100,2), nsmall = 1))
    }
    return(table)
  }
  
  
  count = 0
  table_DS = lapply(...,test, split="DS")
  count = 0
  table_R = lapply(...,test, split="R")
  
  table_DS = Reduce(function(...) merge(..., by="beta"),table_DS)
  table_R = Reduce(function(...) merge(..., by="beta"),table_R)
  
  ## Header names
  if(grepl("t", all_args)){
    header_names=c("")
    count = 0
    for(i in dfs){
      count = count + 1
      header_names = c(header_names, setNames(2,as.character(i)))
      colNames = c("$|{\\beta}_i|$", rep(c("Length", "Coverage (%)"), length(dfs)))
    }
    header_names2 = c("", setNames(2*length(dfs), paste0("Degrees of freedom, k")))
  }else if(grepl("sn", all_args)){
    header_names=c("")
    count = 0
    for(i in xis){
      count = count + 1
      header_names = c(header_names, setNames(2,as.character(i)))
      header_names2 = c("", setNames(2*length(dfs), paste0("Skewness, xi")))
      colNames = c("$|{\\beta}_i|$", rep(c("Length", "Coverage (%)"), length(xis)))
    }
  }else{
    header_names = c("",setNames(2, "Estimated $\\sigma^{2}$"),setNames(2, "Fixed $\\sigma^{2}$"))
    colNames = c("$|{\\beta}_i|$", rep(c("Length", "Coverage (%)"), 2))
  }
  
  
  table_DS_final = kbl(table_DS,align = 'c', format="latex", booktabs = T, escape = F, caption = "DS", 
                       col.names = colNames) %>%
    add_header_above(header_names, bold = T) %>%
    # {if(!grep("normal", all_args)) add_header_above(header_names2, bold = T)} %>%
    column_spec(1, width="5em") %>%
    column_spec(2:5, width = "5em") %>%
    kable_styling(position = "center", latex_options = c("hold_position", "scale_down"), font_size = 12)
  
  table_R_final = kbl(table_R,align = 'c', booktabs = T, format="latex", escape = F, caption = "(U, V)", 
                      col.names = colNames) %>%
    add_header_above(header_names, bold = T) %>%
    # {if(!grep("normal", all_args)) add_header_above(header_names2, bold = T)} %>%
    column_spec(1, width="5em") %>%
    column_spec(2:5, width = "5em") %>%
    kable_styling(position = "center", latex_options = c("hold_position", "scale_down"), font_size = 12)
  
  return(list(out_p=out_p, table_DS_final = table_DS_final, table_R_final = table_R_final))
}



################# Final results #################
############# t #############
t_final_results = generateResult(t_dfs)
# t_final_results[[1]]
t_final_results[[2]] ## DS
t_final_results[[3]] ## (U, V)
# ggsave(plot = t_final_results[[1]], filename = "t_plot.pdf", width = 210, height = 297, units = "mm")


t_pj_final_results = generateResult(t_pj_dfs)
# t_pj_final_results[[1]]
t_pj_final_results[[2]] ## DS
t_pj_final_results[[3]] ## (U, V)
# ggsave(plot = t_pj_final_results[[1]], filename = "t_pj_plot.pdf", width = 210, height = 297, units = "mm")


t_ps = ggarrange(t_final_results[[1]],  t_pj_final_results[[1]],
                 labels = c("A", "B"), nrow = 2, common.legend = TRUE, legend="bottom")
t_ps
# ggsave(plot=t_ps, filename = "t_ps.pdf", width = 210, height = 297, units = "mm")

############# Skewed normal #############
sn_final_results = generateResult(sn_dfs)
# sn_final_results[[1]]
sn_final_results[[2]] ## DS
sn_final_results[[3]] ## (U, V)
# ggsave(plot = sn_final_results[[1]], filename = "sn_plot.pdf", width = 210, height = 297, units = "mm")


sn_pj_final_results = generateResult(sn_pj_dfs)
# sn_pj_final_results[[1]]
sn_pj_final_results[[2]] ## DS
sn_pj_final_results[[3]] ## (U, V)
# ggsave(plot = sn_pj_final_results[[1]], filename = "sn_pj_plot.pdf", width = 210, height = 297, units = "mm")

sn_ps = ggarrange(sn_final_results[[1]],  sn_pj_final_results[[1]],
                  labels = c("A", "B"), nrow = 2, common.legend = TRUE, legend="bottom")
sn_ps
# ggsave(plot=sn_ps, filename = "sn_ps.pdf", width = 210, height = 297, units = "mm")


############# sigma = 5 #############
normal_low_d_final_results = generateResult(normal_low_d_dfs)
# normal_low_d_final_results[[1]]
normal_low_d_final_results[[2]] ## DS
normal_low_d_final_results[[3]] ## (U, V)

normal_high_d_final_results = generateResult(normal_high_d_dfs)
# normal_high_d_final_results[[1]]
normal_high_d_final_results[[2]] ## DS
normal_high_d_final_results[[3]] ## (U, V)

normal_high_sigma_ps = ggarrange(normal_low_d_final_results[[1]],  normal_high_d_final_results[[1]],
                                 labels = c("A", "B"), nrow = 2, common.legend = TRUE, legend="bottom")
normal_high_sigma_ps
# ggsave(plot=normal_high_sigma_ps, filename = "normal_high_sigma_ps.pdf", width = 210, height = 297, units = "mm")


############# sigma = 1 #############
normal_low_d_sigma1_final_results = generateResult(normal_low_d_sigma1_dfs)
# normal_low_d_sigma1_final_results[[1]]
normal_low_d_sigma1_final_results[[2]] ## DS
normal_low_d_sigma1_final_results[[3]] ## (U, V)

normal_high_d_sigma1_final_results = generateResult(normal_high_d_sigma1_dfs)
# normal_high_d_sigma1_final_results[[1]]
normal_high_d_sigma1_final_results[[2]] ## DS  
normal_high_d_sigma1_final_results[[3]] ## (U, V)

normal_low_sigma_ps = ggarrange(normal_low_d_sigma1_final_results[[1]],  normal_high_d_sigma1_final_results[[1]],
                                labels = c("A", "B"), nrow = 2, common.legend = TRUE, legend="bottom")
normal_low_sigma_ps
# ggsave(plot=normal_low_sigma_ps, filename = "normal_low_sigma_ps.pdf", width = 210, height = 297, units = "mm")

