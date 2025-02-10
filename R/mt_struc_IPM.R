#' @title Build structural matrix for independent pathway model
#' @name mt_struc_ipm
#' @author Gordon Grabert, \email{wbonat@@ufpr.br}
#'
#' @description This function builds structural matrices to assist
#' multivariate analysis of twin data. It is an internal function used
#' in the mc_twin function.
#'
#' @param n_resp Number of response variables.
#' @keywords internal
#' @return A list of matrices.


mt_struc_ipm <- function(n_resp) {
  n_resp_c <- n_resp + 1

  # Generate specific matrices
  specific <- lapply(1:n_resp, function(i) {
    mat <- Matrix(0, nrow = n_resp, ncol = n_resp_c, sparse = TRUE)
    mat[i, i + 1] <- 1
    mat
  })

  # Generate common matrices
  common <- lapply(1:n_resp, function(i) {
    mat <- Matrix(0, nrow = n_resp, ncol = n_resp_c, sparse = TRUE)
    mat[i, 1] <- 1
    mat
  })

  return(c(specific, common))
}

