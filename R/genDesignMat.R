#' Return a design matrix
#'
#' This function generates a design matrix of given size,
#' whose rows are independent Gaussian with decaying covariance.
#'
#' @param n Integer value of sample size
#' @param p Integer number of variables
#' @param r Rate (0 to 1) of decay so the (i,j)th covariance is r^|i-j|
#' @return Matrix of size n x p
#' @export
#' @examples
#' X <- genDesignMat(150, 200, 0.5)

genDesignMat <- function(n, p, r){
  V <- matrix(0, nrow=p, ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      V[i,j] <-  r^(abs(i-j))
    }
  }
  MASS::mvrnorm(n, mu = rep(0,p), Sigma = V)
}
