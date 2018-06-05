#' @title Quasi-score function
#' @author Wagner Hugo Bonat
#'
#' @description Compute the quasi-score function, its sensitivy and
#'     variability matrix.
#'
#' @param D A matrix. In general the output from
#'     \code{\link[mcglm]{mc_link_function}}.
#' @param inv_C A matrix. In general the output from
#'     \code{\link[mcglm]{mc_build_C}}.
#' @param y_vec A vector.
#' @param mu_vec A vector.
#' @param weights Additional weights for the estimating function.
#' @keywords internal
#' @return The quasi-score vector, the Sensivity and variability
#'     matrices.
#' @export

ef_quasi_score <- function(D, inv_C, y_vec, mu_vec, weights = NULL) {
    res <- y_vec - mu_vec
    t_D <- Matrix::t(D)
    if(!is.null(weights)) {
      W <- Diagonal(n = length(mu_vec), weights)
      score <- t_D %*% (inv_C %*% W %*% res)
      sensitivity <- -t_D %*% inv_C %*% W %*% D
    }
    if(is.null(weights)) {
      score <- t_D %*% (inv_C %*% res)
      sensitivity <- -t_D %*% inv_C %*% D
    }
    output <- list(Score = score, Sensitivity = sensitivity,
                   Variability = -sensitivity)
    return(output)
}
