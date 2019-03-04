#' Variable selection using different criteria 
#'
#' This function uses specified criteria to select the optimal subset from a list of subsets, 
#' given design matrix X and observation y.
#'
#' @param candidates List of vectors each representing the indices of significant variables
#' @param dat List containing design matrix X and observation y
#' @param criteria Vector of strings ('GIC2, GICn, Cp, AIC, BIC, BC') indicating criteria to use
#' @param penaltyBC Vector of non-default penalty values in using BC, default to NULL so that
#' only the suggested value n^(1/3) will be considered
#' @param adaptiveBC Boolean indicating whether to use adaptive BC, default to FALSE so that
#' only the suggested value n^(1/3) will be considered
#' @param methodPI String ('BC', 'Drop', or 'Both') indicating the method to calculate BC, default to 'BC'
#' @return List (cri_opt) of indices of candidate subsets, each corresponding to a method (cri_name), 
#' along with parametricness index PI for each type of BC
#' @export
#' @examples
#' n <- 150
#' p <- 200
#' X <- genDesignMat(n,p,0.5)
#' beta <- rep(0, p)
#' beta[c(99,199)] = c(10, 5)
#' mu <- X %*% as.matrix(beta, ncol=1)
#' y <- mu + stats::rnorm(n)
#' candidates <- vector('list', 0)
#' candidates[[1]] <- c(1, 99)
#' candidates[[2]] <- c(99, 199) 
#' candidates[[3]] <- c(1, 99, 199)
#' dat <- list(X = X, y= y)
#' res <- varSelection(candidates, dat)

varSelection <- function(candidates, 
                         dat, 
                         criteria=c('BC'), 
                         penaltyBC=NULL, 
                         adaptiveBC=FALSE, 
                         methodPI='BC'){

  X <- dat$X
  y <- dat$y
  nCand <- length(candidates)
  n <- nrow(X)

  # obtain the "consistent" estimator of sigma2 in order to cal. gic2 & gicn
  # only calculable when the number of distinct variables is less than sample size
  vars <- c()
  for (j in 1:nCand){
    vars <- c(vars, candidates[[j]])
  }
  vars <- unique(vars)

  if (length(y) >  length(vars)) { 
    wellCondition <- TRUE
    ls <- stats::lsfit(X[, vars], y, intercept = TRUE)
    sig2 <- mean((ls$residuals)^2) * n / (n-length(vars))    
    
  }else{
    wellCondition <- FALSE
    sig2 <- NA
  }

  # penalty param for BC
  if (is.null(penaltyBC)) { 
    WN <- c(n^(1/3))
  }else { 
    WN <- penaltyBC 
  }

  if (adaptiveBC && wellCondition){
    WN <- c(WN, chooseWNAdaptively(y, X, candidates, sig2, 5))
  }
  # else{
  #   WN <- c(WN, NA)
  # }
  nWN <- length(WN)
  
  # pre-allocate 
  dim <- rep(NA, nCand)
  mse <- rep(NA, nCand)
  aic <- rep(NA, nCand)
  bic <- rep(NA, nCand)
  gic2 <- rep(NA, nCand)
  gicn <- rep(NA, nCand)
  bc <- vector('list', nWN)
  for (i in 1:nWN) { bc[[i]] <- rep(NA, nCand) }
  #cv1 <- rep(NA, nCand)
  #cvd <- rep(NA, nCand)
  for (i in 1:nCand){
    dim[i] <- length(candidates[[i]])
  }
  # the parametricness index, where the last one is defined by Liu and Yang
  PI <- rep(NA, nWN+1) 

  #compute mse
  for (k in 1:nCand){
    if (length(candidates[[k]]) >= (length(y)-1)){
      mse[k] <- Inf
    }else{
      ls <- stats::lsfit(X[,candidates[[k]]], y, intercept = TRUE)
      mse[k] <- mean((ls$residuals)^2)
    }
  }

  # result using method AIC, calculated regardless of input param
  aic <- log(mse) + 2 * dim / n
  aic_opt <- which.min(aic) 

  # result using method BIC, calculated regardless of input param
  bic <- log(mse) + log(n) * dim / n  
  bic_opt <- which.min(bic)   

  bc_opt <- rep(NA, nWN)

  if (wellCondition){
    # result using method GIC2, calculated regardless of input param
    gic2 <- mse + 2 * sig2 * dim / n
    gic2_opt <- which.min(gic2) 
    
    # result using method GICn, calculated regardless of input param
    # gicn uses lambda=log(n) to resemble BIC
    gicn <- mse + log(n) * sig2 * dim / n 
    gicn_opt <- which.min(gicn) 

    # gic2/gicn if well defined and aic/bic otherwise
    surro_gic2 <- gic2_opt
    surro_gicn <- gicn_opt

  }else{ 
    # if not in wellCondition, use AIC, BIC based BC, and return NA for GIC, GICn
    gic2_opt <- NA
    gicn_opt <- NA
    surro_gic2 <- aic_opt
    surro_gicn <- bic_opt
  }

  # result using method BC
  for (i in 1:nWN){
    for (k in 1:nCand){ 
      if (wellCondition){
        bc[[i]][k] <- mse[k] + 2 * sig2 * sum(1/(1:dim[k])) * WN[i] / n
        
        # filter out all the candidates that have dimension larger than gic2's
        if (dim[k] > dim[which.min(gic2)])  { bc[[i]][k] <- Inf }
        
      }else{
        bc[[i]][k] <- log(mse[k]) + 2 * sum(1/(1:dim[k])) * WN[i] / n
        if (dim[k] > dim[which.min(aic)])  { bc[[i]][k] <- Inf }
      }
    }
    bc_opt[i] <- which.min(bc[[i]])
  }
 
  for (i in 1:nWN) { 
    # Note that 1. the candidates (not user-specified but) selected from solution path 
    # may be just the true model 
    # 2. d_BC is not necessarily beweteen d_AIC d_BIC
    if ((surro_gic2 == surro_gicn) && (surro_gic2  == bc_opt[i])) {
      PI[i] <- 1
    }else{
      PI[i] <- abs( dim[bc_opt[i]] - dim[surro_gic2] ) / 
      ( abs( dim[bc_opt[i]] - dim[surro_gic2] ) + abs( dim[bc_opt[i]] - dim[surro_gicn] ) )
    }
  }

  if (wellCondition && (methodPI == 'Drop' || methodPI == 'Both')){ #only in this case is sig2 defined 
    PI[nWN+1] <- getLiuPI( mse, dim, which.min(gicn), sig2, n )
  }else{
    PI[nWN+1] <- NA
  }

  # use string matching to decide the final return 
  # index of subset selected by each criterion
  cri_opt <- c()
  # name of each criterion, for BC of different penalties it looks as 'BC', 'BC 2', ...
  cri_name <- c()
  if (sum(grepl('GIC2', criteria)) > 0) {
    if (!wellCondition) {stop('Due to small sample size, noise variance cannot be estimated for GIC2')}
    cri_name <- c(cri_name, 'GIC2')
    cri_opt <- c(cri_opt, gic2_opt)
  }
  if (sum(grepl('Cp', criteria)) > 0) {#same as GIC2
    if (!wellCondition) {stop('Due to small sample size, noise variance cannot be estimated for Cp')}
    cri_name <- c(cri_name, 'Cp')
    cri_opt <- c(cri_opt, gic2_opt)
  } 
  if (sum(grepl('GICn', criteria)) > 0) {
    if (!wellCondition) {stop('Due to small sample size, noise variance cannot be estimated for GICn')}
    cri_name <- c(cri_name, 'GICn')
    cri_opt <- c(cri_opt, gicn_opt)
  }
  if (sum(grepl('AIC', criteria)) > 0) {
    cri_name <- c(cri_name, 'AIC')
    cri_opt <- c(cri_opt, aic_opt)
  }
  if (sum(grepl('BIC', criteria)) > 0) {
    cri_name <- c(cri_name, 'BIC')
    cri_opt <- c(cri_opt, bic_opt)
  }
  if (sum(grepl('BC', criteria)) > 0) {
    cri_name <- c(cri_name, 'BC')
    if (length(bc_opt) > 1){
      for (j in 2:length(bc_opt)) {
        cri_name <- c(cri_name, paste('BC', j))
      }
    }
    if (adaptiveBC) { cri_name[length(cri_name)] <- 'adaBC' }
    cri_opt <- c(cri_opt, bc_opt)
  }

  res <- list()
  res$cri_name <- cri_name
  res$cri_opt <- cri_opt

  if (methodPI == 'BC') {
    res$PI <- PI[1:nWN]
  }else{
    if (!wellCondition) {stop('Due to small sample size, noise variance cannot be estimated for calculating PI')}
    if (methodPI == 'Both'){ 
      res$PI <- PI[1:(nWN+1)]
    }else{
      res$PI <- PI[nWN+1]
    }
  }

  res
  
}