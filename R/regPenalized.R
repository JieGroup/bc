#' Penalized high dimensional regresssion
#'
#' This function uses a particular penalized regression "MCP", "SCAD", or "lasso"
#' and cross validation to find the regression variables and cofficients,
#' and generate prediction (on a testing design matrix), selected variable support, 
#' a solution path. The function relies on ncvreg package. Intercept is considered. 
#'
#' @param dat List containing design Matrix X and observations Vector y 
#' @param hd_method Strong ("MCP", "SCAD", "lasso") indicating the method in use
#' @param dfmax Integer of maximum number of variables to select, 
#' default to the number of variables
#' @param X_test Matrix of design (for testing)
#' @return pre_opt Vector of predicted observations of size 'number of test size' x 1
#' @return supp_opt Vector of indices of the selected variables 
#' @return Beta Matrix of estimated coefficients along solution path, 
#' of size ('number of variables' + 1) x 'number of lambdas'
#' @return fit Object of from package 'ncvreg' that can be used to predict
#' @export
#' @examples
#' n <- 150
#' p <- 200
#' X <- genDesignMat(n,p,0.5)
#' beta <- rep(0, p)
#' beta[c(99,199)] = c(10, 5)
#' mu <- X %*% as.matrix(beta, ncol=1)
#' y <- mu + stats::rnorm(n)
#' dat <- list(X = X, y= y)
#' res <- regPenalized(dat, c("MCP", "SCAD", "lasso"))

regPenalized <- function(dat, hd_method, dfmax = NULL, X_test = NULL){
  
  X <- dat$X 
  y <- dat$y
  p <- ncol(X)
  if (is.null(dfmax)){
    dfmax <- p
  }

  #fit and choose the best by default 10 fold 
  cvfit <- ncvreg::cv.ncvreg(X, y, penalty = hd_method, dfmax = dfmax)  
  #plot(cvfit)
  #summary(cvfit)
  lambda_opt <- cvfit$lambda.min
  fit <- cvfit$fit
  #summary(fit)
  
  
  if (!is.null(X_test)){
    pre_opt <- stats::predict(fit, X_test, type="response", lambda_opt)
  }else{
    pre_opt <- NULL
  }
  Beta <- fit$beta
  #beta <- fit$beta[,cvfit$min] #first row is the intercept 
  supp_opt <- which(fit$beta[2:(p+1),cvfit$min]!=0)
  
  res <- list(pre_opt = pre_opt, supp_opt = supp_opt, Beta = Beta, fit)

}
