#' Get candidates from marginal ordering of variables
#'
#' This function uses a weighting method to assign weights to 
#' each subset in a list of candidate subsets,
#' and then obtain marginal weights for variables to obtain 
#' a new candidate set with nested subsets. 
#'
#' @param initCandidates List of vectors each indicating an initial subset
#' @param dat List containing design Matrix X and observation Vector y
#' @param weight_method String indicating which method to use for weighting
#' @return candidates List of vectors each representing a new candidate subset
#' @return w_var Vector of each variable's marginal importance  
#' @return ind_var Vector of each variable's index in the original system 1...p
#' 
# @export

getCandidates <- function(initCandidates, 
                          dat, 
                          weight_method = 'ARM'){

  X1 <- dat$X
  y1 <- dat$y

  if (weight_method == 'unchange'){
    
    candidates <- initCandidates
  
  }else{
    n1 <- nrow(X1)
    p <- ncol(X1)

    #random subsample instead of 1:n1 is essential for (non-iid) real data 
    n1a <- floor(n1/2)
    ind_n1a <- sample(1:n1, n1a, replace = FALSE) 
    X1a <- X1[ind_n1a,]
    y1a <- y1[ind_n1a]
    X1b <- X1[-ind_n1a,]
    y1b <- y1[-ind_n1a]
    nCan <- length(initCandidates) 
    wlog <- rep(NA, nCan)

    for (k in 1:nCan){
      if ( length(initCandidates[[k]]) >= (nrow(X1a) - 1) ){
        # give small weights to those cannot be fit
        wlog[k] <- -1e4
      }else{ 
        ls <- stats::lsfit(X1a[,initCandidates[[k]]], y1a, intercept = TRUE)
        d <- length(initCandidates[[k]])
        sig2 <- mean((ls$residuals)^2) * n1a / (n1a - d)
        D <- sum((y1b - ls$coef[1] - matrix(X1b[,initCandidates[[k]]], nrow=n1-n1a) 
              %*% matrix(ls$coef[2:(1+d)], ncol=1))^2)
        C <- d * log(exp(1) * p / d) + 2 * log(d+2)
        wlog[k] <- -C - n1a/2 * log(sig2) - D/(2*sig2)
      }
    }
    # obtain normalized weight of each subset 
    w <- exp(wlog - max(wlog))
    w <- w/sum(w) 
    
    # obtain unique variable index (under the original index systme of initCandidates)
    varInd <- c()
    for (k in 1:nCan){
      varInd <- c(varInd, initCandidates[[k]])
    }
    varInd <- unique(varInd)
    nvar <- length(varInd)

    w_var <- rep(0, nvar)  #weight of each variable, extracted from models 
    for (k in 1:nCan){
      for (j in 1:length(initCandidates[[k]])){
        v <- initCandidates[[k]][j]
        w_var[which(varInd == v)] <- w_var[which(varInd == v)] + w[k]
      }
    }

    # obtain normalized marginal weight of each variable
    w_var <- w_var / sum(w_var) 
    orderVar <- sort(w_var, index.return = TRUE, decreasing = TRUE)$ix
    candidates <- vector('list', nvar)
    for (k in 1:nvar){
      candidates[[k]] <- varInd[orderVar[1:k]]
    }
  }

  list(candidates = candidates, w_var = w_var, ind_var = varInd)
}
