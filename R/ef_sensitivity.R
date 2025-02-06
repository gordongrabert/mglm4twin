#' @title Sensitivity matrix
#' @author Wagner Hugo Bonat and Eduardo Elias Ribeiro Jr
#'
#' @description Compute the sensitivity matrix associated with the
#'     Pearson estimating function.
#'
#' @param product A list of matrix.
#' @param W weights.
#' @return The sensitivity matrix associated with the Pearson estimating
#'     function.
#' @keywords internal
#' @details This function implements the equation 7 of Bonat and
#'     Jorgensen (2016).
#' @useDynLib mglm4twin
#' @importFrom Rcpp sourceCpp

ef_sensitivity <- function(product, W) {
  #sourceCpp("src/mc_sensitivity_op.cpp")
  Sensitivity <- ef_sensitivity_op(products = product, W = W)
  Sensitivity <- forceSymmetric(Sensitivity, uplo = "L")
  return(Sensitivity)
}

#ef_sensitivity <- function(product) {
#    #sourceCpp("src/mc_sensitivity_op.cpp")
#    Sensitivity <- ef_sensitivity_op(products = product)
#    Sensitivity <- forceSymmetric(Sensitivity, uplo = FALSE)
#    return(Sensitivity)
#}

#ef_sensitivity <- function(product) {
#  n_par <- length(product)
#  Sensitivity <- matrix(0, n_par, n_par)
#  Sensitivity_temp <- matrix(0, n_par, n_par)
#  Sensitivity1 <- matrix(0, n_par, n_par)
#  for (i in 1:n_par) {
#    for (j in 1:n_par) {
#      Sensitivity[i, j] <- -sum(Matrix::t(product[[i]]) * product[[j]])
#    }
#  }
#  return(Sensitivity)
#}
