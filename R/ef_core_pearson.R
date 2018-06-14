#' @title Core of the Pearson estimating function.
#' @author Wagner Hugo Bonat
#'
#' @description Core of the Pearson estimating function.
#'
#' @param product A matrix.
#' @param inv_C A matrix.
#' @param C A matrix.
#' @param W A matrix.
#' @param r_w A vector of weighted residuals.
#' @return A vector.
#' @keywords internal
#' @details It is an internal function.

ef_core_pearson <- function(product, inv_C, C, W, r_w) {
    Wlambda <- product%*%inv_C
    output <- Matrix::t(r_w) %*% Wlambda %*% r_w -
      sum(Matrix::diag(Wlambda %*% W %*% C %*% W))
    return(as.numeric(output))
}
