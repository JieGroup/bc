#' Run synthetic experiments for comparing the performance of high dimensional regression
#'
#' This function generates synthetic data and run BCVO, LASSO, or others, 
#' and summarize its performance from several replications. 
#'
#' @param genData Function to generate synthetic training and testing data
#' @param n Integer value of sample size
#' @param p Integer number of variables
#' @param nRep Integer number (>1) of replications in the experiment
#' @param ratio Fraction (0 to 1) of the data used to generate solution paths 
#' @param weight_method String ('simple' or 'ARM') of method to use in weighting subsets 
#' @param fixedBeta Boolean indicator of whether to use fixed design in test
#' @return Printout of means and standard errors of prediction loss, 
#' number of overfitted variables, and number of underfitted variables

runRepExperiment <- function(genData, n, p, nRep, ratio, weight_method, fixedBeta){ 
  
  dfmax <- floor(sqrt(n))
  loss <- c()
  nOver <- c()
  nUnder <- c()
  PI <- c()
  err_mcp <- matrix(NA, nrow = 1000, nRep)
  err_scad <- matrix(NA, nrow = 1000, nRep)
  err_lasso <- matrix(NA, nrow = 1000, nRep)
  err_bc <- matrix(NA, nrow = 1000, nRep)
  
  # MAJOR difference with the ordinary function runRepExperiment(): fix test X
  # if (fixedBeta){
  #   dat_test <- genData(1000, p) 
  # }else{
  #   dat_test <- genData(1000, p, dat$beta_true)
  # }
  
  for (i in 1:nRep){
    #start.time <- Sys.time()
    
    if (fixedBeta){
      dat <- genData(n, p)
      dat_test <- genData(1000, p)
    }else{
      dat <- genData(n, p)
      dat_test <- genData(1000, p, dat$beta_true)
    }
    
    res <- comparePerformance(dat, 
                              dat_test, 
                              hd_methods = c("MCP", "SCAD", "lasso"),
                              dfmax, 
                              ratio, 
                              weight_method,
                              criteria = c('GIC2', 'GICn', 'BC'), 
                              penaltyBC = NULL, 
                              adaptiveBC = TRUE, 
                              methodPI = 'Both')

    loss <- rbind(loss, res$loss)
    nOver <- rbind(nOver, res$nOver)
    nUnder <- rbind(nUnder, res$nUnder)
    PI <- rbind(PI, res$PI)
    
    err_mcp[,i] <- res$err_mcp_x
    err_scad[,i] <- res$err_scad_x
    err_lasso[,i] <- res$err_lasso_x
    err_bc[,i] <- res$err_bc_x
    
    #end.time <- Sys.time()
    #print(paste0('time elapsed: ', end.time - start.time))

  }
  
  #print the result 
  print('The used criteria in order')
  print(c("MCP", "SCAD", "lasso", res$cri_name))
  
  print('loss')
  print(colMeans(loss)) 
  print(matrixStats::colSds(loss)/sqrt(nRep))
  
  print('nOver')
  print(colMeans(nOver)) 
  print(matrixStats::colSds(nOver)/sqrt(nRep))
  
  print('nUnder')
  print(colMeans(nUnder)) 
  print(matrixStats::colSds(nUnder)/sqrt(nRep))
  
  print('PI')
  print(colMeans(PI)) 
  print(matrixStats::colSds(PI)/sqrt(nRep))
  
  #analyze the bias-variance tradeoff of mcp, scad, lasso, bc (in order)
  bias_mcp <- mean((rowMeans(err_mcp))^2)
  bias_scad <- mean((rowMeans(err_scad))^2)
  bias_lasso <- mean((rowMeans(err_lasso))^2)
  bias_bc <- mean((rowMeans(err_bc))^2)
  
  var_mcp <- 0
  var_scad <- 0
  var_lasso <- 0
  var_bc <- 0
  
  for(j in 1:1000){
    var_mcp <- var_mcp + stats::var(err_mcp[j,])/1000
    var_scad <- var_scad + stats::var(err_scad[j,])/1000
    var_lasso <- var_lasso + stats::var(err_lasso[j,])/1000
    var_bc <- var_bc + stats::var(err_bc[j,])/1000
  }
  
  biases <- c(bias_mcp, bias_scad, bias_lasso, bias_bc)
  vares <- c(var_mcp, var_scad, var_lasso, var_bc)
  print('Bias of mcp, scad, lasso, bc')
  print(round(biases, 1))
  
  print('Var of mcp, scad, lasso, bc')
  print(round(vares, 1))
  
  print('Var percentage of mcp, scad, lasso, bc')
  print(vares/(vares+biases) * 100)
  
}


 
 