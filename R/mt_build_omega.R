#' @title Build omega matrix
#' @name mt_build_omega
#' @author Wagner Hugo Bonat
#'
#' @description This function builds \eqn{\Omega} matrix.
#'
#' @param tau A vector
#' @param Z A list of matrices in general obtained from mt_twin function.
#' @keywords internal
#' @return A list with the \eqn{\Omega} matrix and its derivatives 
#' with respect to \eqn{\tau}.

mt_build_omega <- function(tau, Z) {
        Omega <- mt_matrix_linear_predictor(tau = tau, Z = Z)
        output <- list(Omega = Omega, D_Omega = Z)
    return(output)
}
