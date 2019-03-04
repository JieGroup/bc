#' Get supports, i.e. the indices of selected variables, from a solution path
#'
#' This function uses a matrix Beta that represents a solution path to calculate
#' all the supports along solution path. 
#'
#' @param Beta Matrix of size ('number of variables' + 1) x 'number of lambdas'
#' @return List of 'number of lambdas' vectors, each representing a support 
#' (excluding intercept term)
# @export

getSupp <- function(Beta) {
  nLambda <- ncol(Beta)
  supp <- vector("list", nLambda) 
  for (i in 1:nLambda){
    #original index start with interception term which we need to remove
    supp[[i]] <- setdiff( which(Beta[,i]!=0), 1 )
    supp[[i]] <- supp[[i]] - 1 
  }
  if (length(supp[[1]]) == 0){
    supp[[1]] <- {}
  }
  supp
}
