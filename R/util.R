#=============== subroutines for high dim BC ====================
# setwd('/Users/Jie/Dropbox/MyResearch2016/Bridge Criterion/R for high dim/bc')
# install.packages('matrixStats')
# install.packages('MASS')
# install.packages('ncvreg')
library('matrixStats')
library('MASS')
library('ncvreg')


#get the candidate set
mergeSupp <- function(supp1, supp2) {
  unique( c(supp1, supp2) )
}



#compute the number of overfitted and underfitted vars, based on selected supp
getDeviation <- function(supp, supp0){
  nOver <- setdiff(supp, supp0)
  nUnder <- setdiff(supp0, supp)
  c(nOver, nUnder)
}

inverseLogit <- function(x) {
  2 * (exp(x)-1) / (exp(x)+1)
}

plotMu <- function(mu, scale, skew) { 
  graphics::plot(mu, scale * inverseLogit(mu/scale-skew), lty=3, col="blue", xlab="mu", ylab="g(mu)") 
  graphics::lines(mu, mu, lty=1, col="black")
}


#get the PI defined by Liu and Yang before 
#mse, dim are the MSE and dimension of each candidate model 
#k_gicn is the model index selected by GIC_n -- a consistent procedure
#sig2 is the estimated noise level 
#n is the sample size 
getLiuPI <- function( mse, dim, k_gicn, sig2, n ){
  if (dim[k_gicn] == min(dim)){
    PI <- n
  }else{
    ind_dim_smaller <- which(dim < dim[k_gicn])
    nextDim <- max( dim[ind_dim_smaller] )
    ind_sub_model <- which( dim == nextDim )
    n_ind_sub_model <- length(ind_sub_model)
    pis = rep(0,n_ind_sub_model)
    for (j in 1:n_ind_sub_model){
      a <- ind_sub_model[j]
      pis[j] <- (mse[a] * n + dim[a] * sig2 * log(n) - n * sig2) / 
        (mse[k_gicn] * n + dim[k_gicn] * sig2 * log(n) - n * sig2)
    }
    PI <- min(pis)
  }
  PI  
}
 
# calculate several optimum of data driven wn and average
chooseWNAdaptively <- function(y, X, candidates, sig2, Iter){
  
  WN_opts <- rep(0,Iter)
  numCand <- length(candidates)
  n <- length(y)
  mse <- rep(0,numCand)
  bc_tr <- rep(0,numCand)
  dim <- rep(0,numCand)
  for (k in 1:numCand){
    dim[k] <- length(candidates[[k]])
  }

  for (i in 1:Iter){
    
    #Method 1: random split
    #P <- randperm(n)
    #n_temp <- floor(n/2)
    #ind_tr <- P(1:n_temp)
    #ind_te <- P(n_temp+1:end)
  
    #Method 2: resample (bootstrap) 
    ind_tr <- sample(1:n, n, replace = TRUE) 
    ind_te <- 1:n
  
    y_tr <- y[ind_tr]
    X_tr <- X[ind_tr,]
    y_te <- y[ind_te]
    X_te <- X[ind_te,]
    
    n_tr <- length(ind_tr)
    n_te <- length(ind_te)
  
    WNs <- n / 1.5^(floor(log(n)/log(1.5)):1)
    e <- rep(0,length(WNs))
    for (j in 1:length(WNs)){ 
      WN <- WNs[j]
      for (k in 1:numCand){ 
        ls <- stats::lsfit(X_tr[,candidates[[k]]], y_tr, intercept = TRUE)
        mse[k] <- mean((ls$residuals)^2)  
        bc_tr[k] <- mse[k] + 2 * sig2 * sum(1./(1:dim[k])) * WN / n_tr
      }
      k_tr <- which.min(bc_tr)
      ls <- stats::lsfit(X_tr[,candidates[[k_tr]]], y_tr, intercept = TRUE)
      nd <- length(candidates[[k_tr]])
      e[j] <- sum((y_te - ls$coef[1] - matrix(X_te[,candidates[[k_tr]]], ncol=nd) %*% matrix(ls$coef[2:(1+nd)], ncol=1)  )^2) / n_te
    }
    ind_opt <- which.min(e)
    WN_opts[i] <- WNs[ind_opt]
  }

  WN <- exp(mean(log(WN_opts)))

}