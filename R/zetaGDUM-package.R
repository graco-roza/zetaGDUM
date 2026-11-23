#' zetaGDUM: zeta diversity with abundance and multi-site uniqueness models
#'
#' Functions to (i) compute zeta diversity and species weights at any order for
#' presence-absence, abundance and dominance data, (ii) locate the breakpoint of
#' the zeta retention rate, (iii) partition set-wise similarity into turnover and
#' nestedness components, and (iv) fit a multi-site extension of generalised
#' dissimilarity uniqueness models (GDUM) with paired bootstrap importance.
#'
#' @useDynLib zetaGDUM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rnorm rpois runif rgamma sd var cov cor median quantile plogis qlogis rlnorm
#' @importFrom utils combn
#' @importFrom parallel mclapply
"_PACKAGE"
