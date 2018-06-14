#' @title Quasi-score function
#' @author Wagner Hugo Bonat
#'
#' @description Compute the quasi-score function, its sensitivy and
#'     variability matrix.
#'
#' @param D A matrix. In general the output from
#'     \code{\link[mcglm]{mc_link_function}}.
#' @param inv_C A matrix. In general the output from
#'     \code{\link[mglm4twin]{mc_build_sigma}}.
#' @param C A matrix. In general the output from
#'     \code{\link[mglm4twin]{mc_build_sigma}}.
#' @param y_vec A vector.
#' @param mu_vec A vector.
#' @param W A matrix od weights for the quasi-score estimating function.
#' @keywords internal
#' @return The quasi-score vector, the Sensivity and variability
#'     matrices.
#' @export

ef_quasi_score <- function(D, inv_C, C, y_vec, mu_vec, W) {
    res <- y_vec - mu_vec
    t_D <- Matrix::t(D)
    t_D_inv_C_W <-t_D %*% inv_C %*% W
    score <- t_D_inv_C_W %*% res
    sensitivity <- -t_D_inv_C_W %*% D
    variability <- t_D_inv_C_W %*% C %*% t(t_D_inv_C_W)
    output <- list(Score = score, Sensitivity = sensitivity,
                   Variability = variability)
    return(output)
}
