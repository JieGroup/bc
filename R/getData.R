#====================== different ways to generate synthetic data ====================

#' Generate data using support of size 11
#'
#' This function generates synthetic data, corresponding to model M5 in the paper. 
#'
#' @param n Integer value of sample size
#' @param p Integer number of variables
#' @param beta_true Input value of coefficients
#' @return List of generated design matrix X, observation y, expected observation mu, 
#' and support supp
# @export

genData_11 <- function(n, p, beta_true = NULL){
  sigma <- 1.5 #5
  X <- genDesignMat(n,p,0)
  beta <- c(10, 5, 5, 2.5, 2.5, 1.25, 1.25, 0.675,0.675,0.3125,0.3125, rep(0, p-11))
  mu <- X %*% as.matrix(beta, ncol=1)
  y <- mu + stats::rnorm(n, mean = 0, sd = sigma)
  list(X=X, y=y, mu=mu, supp=c(1:11))
}

#compute the true snr of the above model, about 50
# snr <- rep(NA, 100)
# for (k in 1:100){
#   X <- genDesignMat(n,p,0)
#   beta <- c(10, 5, 5, 2.5, 2.5, 1.25, 1.25, 0.675,0.675,0.3125,0.3125, rep(0, p-11))
#   mu <- X %*% as.matrix(beta, ncol=1)
#   snr[k] <- mean(mu^2)/sigma^2
# }
#mean(snr)
#sd(snr)/10
#73.17967 (0.8017347) for n = 150





#' Generate data using support of size 15
#'
#' This function generates synthetic data, corresponding to model M6 in the paper. 
#'
#' @param n Integer value of sample size
#' @param p Integer number of variables
#' @param beta_true Input value of coefficients
#' @return List of generated design matrix X, observation y, expected observation mu, 
#' and support supp
# @export

genData_15 <- function(n, p, beta_true = NULL){
  sigma <- 1.5
  X <- matrix(NA, nrow=n, ncol=p)
  X[,1:15] <- genDesignMat(n,15,0.5)
  X[,16:p] <- genDesignMat(n,(p-15),0)
  beta <- c(rep(2.5, 5), rep(1.5, 5), rep(0.5, 5), rep(0, p-15))
  mu <- X %*% matrix(beta, ncol=1)
  y <- mu + stats::rnorm(n, mean = 0, sd = sigma)
  list(X=X, y=y, mu=mu, supp=c(1:15))
}
#snr_1 <- mean(mu^2)/sigma^2 
#snr_0 <- mean((scale*inverseLogit(mu/scale))^2)/sigma^2
#plotMu(mu, scale)
#graphics::plot(mu0, mu, lty=3, col="blue", xlab="mu", ylab="g(mu)") 
#graphics::lines(mu0, mu0, lty=1, col="black")

#compute the true snr of the above model, about 50
# snr <- rep(NA, 100)
# for (k in 1:100){
#   X <- matrix(NA, nrow=n, ncol=p)
#   X[,1:15] <- genDesignMat(n,15,0.5)
#   X[,16:p] <- genDesignMat(n,(p-15),0)
#   beta <- c(rep(2.5, 5), rep(1.5, 5), rep(0.5, 5), rep(0, p-15))
#   mu <- X %*% as.matrix(beta, ncol=1)
#   snr[k] <- mean(mu^2)/sigma^2
# }
# mean(snr)
# sd(snr)/10
#51.89004 (0.5352169) for n=150



#' Generate data using support of size 1
#'
#' This function generates synthetic data, corresponding to model M7 in the paper. 
#'
#' @param n Integer value of sample size
#' @param p Integer number of variables
#' @param beta_true Input value of coefficients
#' @return List of generated design matrix X, observation y, expected observation mu, 
#' and support supp
# @export

#general rule1: large supp0 requires heavier penalty in BC (due to finite sample size)
#general rule2: higher SNR and smaller active support favor BIC-type even in the mis-specified case 
genData_1 <- function(n, p, beta_true = NULL){
  sigma <- 1.5 
  X <- genDesignMat(n,p,0)
  beta <- c(rep(10.5,1), rep(0, p-1)) 
  mu <- X %*% as.matrix(beta, ncol=1)
  y <- mu + stats::rnorm(n, mean = 0, sd = sigma)
  list(X=X, y=y, mu=mu, supp=c(1:1))
}

#snr_1 <- mean(mu^2)/sigma^2 #4 for supp0=1:5
#snr_0 <- mean((scale*inverseLogit(mu/scale-skew))^2)/sigma^2 #1.5 for supp0=1:5
#plotMu(mu, scale, skew)
#graphics::plot(mu0, mu, lty=3, col="blue", xlab="mu", ylab="g(mu)") 
#graphics::lines(mu0, mu0, lty=1, col="black")
#compute the true snr of the above model, about 50
# snr <- rep(NA, 100)
# for (k in 1:100){
#   X <- genDesignMat(n,p,0)
#   #X[,1:15] <- genDesignMat(n,15,0.5)
#   beta <- c(rep(10.5, 1), rep(0, p-1)) #general rule: higher SNR and smaller active support favor BIC-type even in the mis-specified case 
#   mu <- X %*% as.matrix(beta, ncol=1)
#   snr[k] <- mean(mu^2)/sigma^2
# }
# mean(snr)
# sd(snr)/10
#49.11117 (0.5793626)




#' Generate data using support of size p
#'
#' This function generates synthetic data, corresponding to model M7 in the paper. 
#'
#' @param n Integer value of sample size
#' @param p Integer number of variables
#' @param beta_true Input value of coefficients
#' @return List of generated design matrix X, observation y, expected observation mu, 
#' and support supp
# @export

genData_p <- function(n, p, beta_true = NULL){
  sigma <- 5 #1.5 
  X <- genDesignMat(n,p,0)
  beta <- 10/c(1:p)^1.5
  mu <- X %*% as.matrix(beta, ncol=1)
  y <- mu + stats::rnorm(n, mean = 0, sd = sigma)
  list(X=X, y=y, mu=mu, supp=c(1:p))
}
#compute the true snr of the above model, about 50
# snr <- rep(NA, 100)
# for (k in 1:100){
#   sigma <- 5 
#   X <- genDesignMat(n,p,0)
#   beta <- 10/c(1:p)^1.5
#   mu <- X %*% as.matrix(beta, ncol=1)
#   snr[k] <- mean(mu^2)/sigma^2
# }
# mean(snr)
# sd(snr)/10
#53.35 (0.666478) for n=150

