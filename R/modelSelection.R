#' Model selection using different criteria 
#'
#' This function uses specified criteria to select the optimal subset from a list of subsets, 
#' given design matrix X and observation y.
#'
#' @param loglik Vector of (unnormalized) log-likelihood values evaluated for each model
#' evaluted at the maximum likelihood estimator
#' @param dim Vector of dimension (or number of free unknown parameters) of each model
#' @param n Integer of sample size
#' @param criteria Vector of strings ('AIC, BIC, BC') indicating criteria to use
#' @param penaltyBC Vector of non-default penalty values in using BC, default to NULL so that
#' only the suggested value n^(1/3) will be considered
#' @return Vector (cri_opt) of indices of candidate models
#' @return Vector (cri_name) of the criteria used
#' @return Vector (PI) of parametricness index, default to one value that corresponds to 
#' the default BC penalty
#' @export
#' @examples
#' res <- modelSelection(loglik = c(4.1, 5.2, 6.3), dim = c(1, 2, 3), n = 100, 
#' criteria='AIC, BIC, BC', penaltyBC=NULL)

modelSelection <- function(loglik, 
                         dim, 
                         n,
                         criteria, 
                         penaltyBC=NULL){

  nCand <- length(loglik)

  # penalty param for BC
  if (is.null(penaltyBC)) { 
    WN <- c(n^(1/3))
  }else { 
    WN <- penaltyBC 
  }
 
  nWN <- length(WN)
  
  # pre-allocate 
  aic <- rep(NA, nCand)
  bic <- rep(NA, nCand)
  bc <- vector('list', nWN)
  for (i in 1:nWN) { bc[[i]] <- rep(NA, nCand) }

  # the parametricness index based on BC
  PI <- rep(NA, nWN) 

  # result using method AIC, calculated regardless of input param
  aic <- -2 * loglik + 2 * dim
  aic_opt <- which.min(aic) 

  # result using method BIC, calculated regardless of input param
  bic <- -2 * loglik + log(n) * dim
  bic_opt <- which.min(bic) 

  bc_opt <- rep(NA, nWN)

  # result using method BC
  for (i in 1:nWN){
    for (k in 1:nCand){ 
        bc[[i]][k] <- -2 * loglik[k] + 2 * sum(1/(1:dim[k])) * WN[i]
        if (dim[k] > dim[aic_opt])  { bc[[i]][k] <- Inf }
    }
    bc_opt[i] <- which.min(bc[[i]])
  }
 
  for (i in 1:nWN) { 
    # note that d_BC is not necessarily beweteen d_AIC d_BIC
    if ((aic_opt == bic_opt) && (aic_opt  == bc_opt[i])) {
      PI[i] <- 1
    }else{
      PI[i] <- abs( dim[bc_opt[i]] - dim[aic_opt] ) / 
      ( abs( dim[bc_opt[i]] - dim[aic_opt] ) + abs( dim[bc_opt[i]] - dim[bic_opt] ) )
    }
  }

  # use string matching to decide the final return 
  # index of subset selected by each criterion
  cri_opt <- c()
  # name of each criterion, for BC of different penalties it looks as 'BC', 'BC 2', ...
  cri_name <- c()
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
    cri_opt <- c(cri_opt, bc_opt)
  }

  res <- list()
  res$cri_name <- cri_name
  res$cri_opt <- cri_opt
  res$PI <- PI[1:nWN]

  res
  
}