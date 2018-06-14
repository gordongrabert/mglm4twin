#' @title Core of the Pearson estimating function.
#' @author Wagner Hugo Bonat
#'
#' @description Core of the Pearson estimating function.
#'
#' @param product A matrix.
#' @param inv_C A matrix.
#' @param C A matrix.
#' @param res A vector of weighted residuals.
#' @return A vector.
#' @keywords internal
#' @details It is an internal function.

ef_core_pearson <- function(product, inv_C, res) {
  output <- t(res) %*% product %*%
    (inv_C %*% res) - sum(diag(product))
  return(as.numeric(output))
}
