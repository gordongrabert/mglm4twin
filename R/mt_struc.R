#' @title Build structural matrix
#' @name mt_struc
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description This function builds structural matrices to assist
#' multivariate analysis of twin data. It is an internal function used
#' in the mc_twin function.
#'
#' @param n_resp Number of response variables.
#' @keywords internal
#' @return A list of matrices.

mt_struc <- function(n_resp){
  position <- combn(n_resp, 2)
  off_diag <- list()
  n_par <- n_resp*(n_resp-1)/2
  diag <- list()
  for(i in 1:n_resp) {
    diag_temp <- Matrix(0, ncol = n_resp, nrow = n_resp)
    diag_temp[i,i] <- 1
    diag[[i]] <- diag_temp
  }
  for(i in 1:n_par){
    off_diag_temp <- Matrix(0, ncol = n_resp, nrow = n_resp)
    off_diag_temp[position[1,i],position[2,i]] <- off_diag_temp[position[2,i],position[1,i]] <- 1
    off_diag[[i]] <- off_diag_temp
  }
  return(c(diag, off_diag))
}
