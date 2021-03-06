% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regBCVO.R
\name{regBCVO}
\alias{regBCVO}
\title{Run BCVO regression for high dimensional data}
\usage{
regBCVO(dat, hd_methods = c("MCP", "SCAD", "lasso"), dfmax = NULL,
  ratio = 0.7, weight_method = "ARM", criteria = c("BC"),
  penaltyBC = NULL, adaptiveBC = FALSE, methodPI = "BC",
  dat_test = NULL)
}
\arguments{
\item{dat}{List containing design matrix X and observation y}

\item{hd_methods}{Vector of strings ("MCP", "SCAD", "lasso") indicating which
penalized regression methods are used to generate initial candidate subsets}

\item{dfmax}{Integer of maximum number of variables selected by each 
penalized regression, default to floor(sqrt(n)) where n is the sample size}

\item{ratio}{Value (0 to 1) of proportion of data in obtaining initial 
candidate subsets}

\item{weight_method}{String indicating which method to use for weighting}

\item{criteria}{Vector of strings ('GIC2, GICn, Cp, AIC, BIC, BC') indicating criteria to use
in selecting the final subset}

\item{penaltyBC}{Vector of non-default penalty values in using BC, default to NULL so that
only the suggested value n^(1/3) will be considered.  If other values are given, the default 
will be overwritten.}

\item{adaptiveBC}{Boolean indicating whether to use adaptive BC, default to FALSE so that
only the values in penaltyBC will be considered. If TRUE then penaltyBC is appened 
with another data-driven value}

\item{methodPI}{String ('BC', 'Drop'. or 'Both') indicating the method to calculate BC, default to 'BC'}

\item{dat_test}{Matrix of design that is used for comparing performance, default to NULL}
}
\value{
var_imp PI Vector of variable importance

cri_name Vector of criteria used

cri_selection List of vectors of the selected subsets

cri_coefs List of vectors of the effective coefficients

PI Vector of PI corresponding to each BC penalty in 'penaltyBC'

err_bc_x Vector of residuals

loss Vector of prediction loss on test data, default to NA

corr Vector of correlation between prediction and truth on test data, default to NA

nOver Vector of the numbers of overfitted variables on test data, default to NA

nUnder Vector of the numbers of underfitted variables on test data, default to NA
}
\description{
This function first uses a weighting method to assign weights to 
each subset in a list of candidate subsets, which are in turn 
generated by penalized regression methods, then run BC on
the remaining data to obtain the significant subset, 
and finally refit all the data with least squares.
}
\examples{
n <- 150
p <- 200
X <- genDesignMat(n,p,0.5)
beta <- rep(0, p)
beta[c(99,199)] = c(10, 5)
mu <- X \%*\% as.matrix(beta, ncol=1)
y <- mu + stats::rnorm(n)
dat <- list(X = X, y= y)
res <- regBCVO(dat)
}
