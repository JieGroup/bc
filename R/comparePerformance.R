#' Compare MCP, SCAD, LASSO, and BCVO for high dimensional regression
#'
#' This function compares the performance of MCP, SCAD, LASSO, and BCVO on
#' on a single dataset, in terms of prediction loss and over/under-fitting.
#'
#' @param dat Matrix of design that is used for training
#' @param dat_test Matrix of design that is used for comparing performance
#' @param hd_methods Vector of strings ("MCP", "SCAD", "lasso") indicating which penalized regression
#' methods are used to generate initial candidate subsets
#' @param dfmax Integer of maximum number of variables selected by each 
#' penalized regression
#' @param ratio Value (0 to 1) of proportion of data in obtaining initial 
#' candidate subsets
#' @param weight_method String indicating which method to use for weighting
#' @param criteria Vector of strings ('GIC2, GICn, Cp, AIC, BIC, BC') indicating criteria to use
#' in selecting the final subset
#' @param penaltyBC Vector of non-default penalty values in using BC, default to NULL so that
#' only the suggested value n^(1/3) will be considered
#' @param adaptiveBC Boolean indicating whether to use adaptive BC, default to FALSE so that
#' only the suggested value n^(1/3) will be considered
#' @param methodPI String ('BC', 'Drop', or 'Both') indicating the method to calculate BC, default to 'BC'
#' @return List of prediction loss, prediction correlation, number of overfitted variables,  
#' number of underfitted variables, prediction residual, and selected subsets
# @export

comparePerformance <- function(dat, 
                                dat_test, 
                                hd_methods,
                                dfmax, 
                                ratio, 
                                weight_method,
                                criteria = c('GIC2', 'AIC', 'BIC', 'BC'), 
                                penaltyBC = NULL, 
                                adaptiveBC = TRUE, 
                                methodPI = 'BC'){

  X <- dat$X 
  y <- dat$y
  X_test <- dat_test$X 
  # use mu instead of y to calculate loss
  mu_test <- dat_test$mu  
  supp0 <- dat_test$supp
  
  # performance metrics
  n <- nrow(X)
  p <- ncol(X)
  err_mcp_x <- NULL
  err_scad_x <- NULL 
  err_lasso_x <- NULL
  err_bc_x <- NULL
  
  # fit three popular high dim methods and get the mean prediction error on a test data  
  # note that we use p instead of dfmax here, 
  # for authentic comparison, we should not put upper bound dfmax
  nHDMethod <- 3
  res_mcp <- regPenalized(dat, "MCP", p, X_test)
  res_scad <- regPenalized(dat, "SCAD", p, X_test)
  res_lasso <- regPenalized(dat, "lasso", p, X_test)
  
  if (length(X_test) == 0){ 
    stop('X_test has size 0 in compare_performance!')
  }

  err_mcp_x <- mu_test-res_mcp$pre_opt
  err_scad_x <- mu_test-res_scad$pre_opt
  err_lasso_x <- mu_test-res_lasso$pre_opt
  loss <- c( mean((err_mcp_x)^2), mean((err_scad_x)^2), mean((err_lasso_x)^2) )
  corr <- c( stats::cor(mu_test, res_mcp$pre_opt), 
         stats::cor(mu_test, res_scad$pre_opt), 
         stats::cor(mu_test, res_lasso$pre_opt) )
  nOver <- c( length(setdiff(res_mcp$supp_opt, supp0)), 
              length(setdiff(res_scad$supp_opt, supp0)), 
              length(setdiff(res_lasso$supp_opt, supp0)) )
  nUnder <- c( length(setdiff(supp0, res_mcp$supp_opt)), 
               length(setdiff(supp0, res_scad$supp_opt)), 
               length(setdiff(supp0, res_lasso$supp_opt)) )

  res <- regBCVO(dat, hd_methods, dfmax, ratio, weight_method,
                  criteria, penaltyBC, adaptiveBC, methodPI, dat_test)

  loss <- c(loss, res$loss)
  corr <- c(corr, res$corr)
  nOver <- c(nOver, res$nOver)
  nUnder <- c(nUnder, res$nUnder)
  
  #the final loss for each method 
  list(loss = loss, nOver = nOver, nUnder = nUnder, 
       corr = corr, PI = res$PI,
       mcp_supp = res_mcp$supp_opt, scad_supp = res_scad$supp_opt, 
       lasso_supp = res_lasso$supp_opt, 
       cri_name = res$cri_name, cri_selection = res$cri_selection,
       err_mcp_x = err_mcp_x, err_scad_x = err_scad_x,
       err_lasso_x = err_lasso_x, err_bc_x = res$err_bc_x
      )
}
