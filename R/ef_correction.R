#' @title Pearson correction term
#' @author Wagner Hugo Bonat
#'
#' @description Compute the correction term associated with the Pearson
#'     estimating function.
#'
#' @param D_C A list of matrices.
#' @param inv_J_beta A matrix. In general it is computed based on the
#'     output of the \code{[mglm4twin]{ef_quasi_score}}.
#' @param D A matrix. In general it is the output of the
#'     \link{mt_link_function}.
#' @param inv_C A matrix. In general the output of the
#'     \link{mt_build_sigma}.
#' @return A vector with the correction terms to be used on the Pearson
#'     estimating function.
#' @keywords internal
#' @details It is an internal function useful inside the fitting
#'     algorithm.

ef_correction <- function(D_C, inv_J_beta, D, inv_C) {
    term1 <- lapply(D_C, ef_sandwich, bord1 = Matrix::t(D) %*% inv_C,
                    bord2 = inv_C %*% D)
    output <- lapply(term1,
                     function(x, inv_J_beta) sum(x * inv_J_beta),
                     inv_J_beta = inv_J_beta)
    return(unlist(output))
}
