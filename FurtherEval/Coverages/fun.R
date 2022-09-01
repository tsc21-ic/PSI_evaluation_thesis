## Standardize X: center and scale
std <- function(X){
  n <- nrow(X)
  s <- apply(X, 2, sd) * sqrt((n-1)/n)   ## (n-1)/n 
  val <- scale(X, scale=s)
  dimnames(val) <- dimnames(X)
  val
}

# Selection algorithms
Selection = function(y, X, type, p, lam = NULL, first=FALSE, signs=FALSE){
  if (type == "stabs.lasso"){
    p = dim(X)[2]
    sel = selected(stabsel(X, y, fitfun = lars.lasso, PFER = 3, cutoff = 0.7))
    out = rep(F, p)
    out[sel] = TRUE
    # out = c(TRUE, out)
    return(out)
  } 
  else if (type == "lasso"){
    fit = cv.glmnet(X, y, intercept=FALSE)   ## Use cross-validation to select the best lambda
    sel = as.matrix(coef(fit, s=fit$lambda.min)[-1]) != 0
    # sel = as.integer(sub('.', '', rownames(tst)[tst]))
    out = rep(F, p)
    out[sel] = TRUE
    return(out)
    # out = c(TRUE, out)     ## TRUE => selected variable , FALSE => not selected variable
  } else if (type == "lasso.fixed.lambda"){
    fit = glmnet(X, y, intercept = FALSE, standardize = TRUE)
    if(first){
      n = length(y)
      lam = sqrt(2*log(p)/n)*n   ## universal lambda parameter (liu et.al 2018)
    }
    coefs = coef(fit, s = lam/n)[-1]
    sel = which(coefs != 0)
    out = rep(F, p)
    out[sel] = TRUE      ## TRUE => selected variable , FALSE => not selected variable
    if(first) return(list(lambda=lam,coefs=coefs, vars=out))
    ifelse(signs, return(list(coefs=coefs, vars=out)), return(out))
  }
}


## Normal Confidence interval when true sigma is applied
confint.truesigma <- function (object, parm, level = 0.95, truesigma, df=NULL, estimate=FALSE)
{
  cf = coef(object)
  pnames = names(cf)
  if(missing(parm)) parm = pnames
  else if(is.numeric(parm)) parm = pnames[parm]
  a = (1 - level)/2
  a = c(a, 1 - a)
  pct = format(a, 3)
  if(estimate){
    fac = qt(a, df=df)
  }else{
    fac = qnorm(a)
  }
  ci = array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  invXtX = solve(crossprod(model.matrix(object)), tol = 1e-50)
  manual = invXtX * (truesigma^2)
  ses = sqrt(diag(manual))[parm]
  ci[] = cf[parm] + ses %o% fac
  return(ci)
}




# Estimator of the standard deviation
sigma_hat = function(y, X){
  if(is.null(dim(X))){
    n = length(X)
    p = 1
  }else{
    n = dim(X)[1]
    p = dim(X)[2]
  }
  if (p < n/4){
    # X_full = cbind(rep(1, n), X)
    ## no intercept term since our we assume no intercept in X
    X_full = X
    # beta_hat_full = solve(t(X_full)%*%X_full, tol=1e-51)%*%t(X_full)%*%y
    beta_hat_full = solve(t(X_full)%*%X_full)%*%t(X_full)%*%y
    df = n - p
    sigmahat = sqrt(sum((y - X_full%*%beta_hat_full)^2)/df)
    return(list(sigmahat=sigmahat, df=df))
  }else{
    est = estimateSigma2(X, y, intercept=FALSE)
    return(list(sigmahat=est$sigmahat, df=est$df))
  }
}

## estimateSigma temporary function
estimateSigma2 <- function(x, y, intercept=TRUE, standardize=TRUE){
  # X; y; intercept=FALSE; standardize=TRUE
  # checkargs.xy(x,rep(0,nrow(x)))
  if(nrow(x)<10) stop("Number of observations must be at least 10 to run estimateSigma")
  cvfit=cv.glmnet(x,y,intercept=intercept,standardize=standardize)
  lamhat=cvfit$lambda.min
  fit=glmnet(x,y,standardize=standardize, intercept=intercept)
  yhat=predict(fit,x,s=lamhat)
  nz=sum(predict(fit,s=lamhat, type="coef")[-1, ] !=0)
  df = max(length(y)- nz, 1)  ## prevent from dividing by 0
  sigma=sqrt(sum((y-yhat)^2)/(df))
  return(list(sigmahat=sigma, df=df))
}



## Indices for data splitting n1
I = function(X, n, n1){
  n1_idx = sample(1:n, size = n1, replace = FALSE)
  return(n1_idx)
}




