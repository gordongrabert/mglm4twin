#' @title Derivatives of V^{1/2} with respect to beta.
#' @author Wagner Hugo Bonat,  \email{wbonat@@ufpr.br}
#'
#' @description Compute the derivatives of \eqn{V^{1/2}} matrix with
#'     respect to the regression parameters beta.
#'
#' @param D A matrix.
#' @param D_V_sqrt_mu A matrix.
#' @return A list of matrices, containg the derivatives of \eqn{V^{1/2}}
#'     with respect to the regression parameters.

mt_derivative_V_sqrt_beta <- function(D, D_V_sqrt_mu) {
    n_beta <- dim(D)[2]
    n_obs <- dim(D)[1]
    output <- list()
    for (i in 1:n_beta) {
      output[[i]] <- Diagonal(n_obs, D_V_sqrt_mu * D[, i])
    }
    return(output)
}

