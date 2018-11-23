#' @title Build twin matrix linear predictor components
#' @name mt_twin
#' @author Wagner Hugo Bonat, \email{wbonat@ufpr.br}
#'
#' @description This function provides the components of the matrix
#' linear predictor suitable for fitting multivariate twin models.
#'
#' @param N_DZ Number of DZ twin pairs.
#' @param N_MZ Number of MZ twin pairs.
#' @param n_resp Number of response variables.
#' @param model a string specifying the name of the twin model.
#' Options are: \code{"E"}, \code{"AE"}, \code{"CE"},
#'     \code{"ACE"} and \code{"ADE"} models.
#' @param formula formula specifying regression model for each dispersion component.
#' @param data Data set.
#' @keywords internal
#' @export
#' @return A list of matrices.

mt_twin <- function(N_DZ, N_MZ, n_resp, model, formula = NULL, data = NULL) {
  # Non-diagonal elements
  a_dz <- 0.5
  d_dz <- 0.25
  ####################################################################
  ## Univariate matrices #############################################
  ####################################################################
  # MZ twin
  MZ_struc <- list()
  MZ_struc$A <- Matrix(c(1, 1, 1, 1), 2, 2)
  MZ_struc$C <- MZ_struc$A
  MZ_struc$D <- MZ_struc$A
  MZ_struc$E <- Diagonal(2, 1)
  # DZ twin
  DZ_struc <- list()
  DZ_struc$A <- Matrix(c(1, a_dz, a_dz, 1), 2, 2)
  DZ_struc$C <- MZ_struc$A
  DZ_struc$D <- Matrix(c(1, d_dz, d_dz, 1), 2, 2)
  DZ_struc$E <- Diagonal(2, 1)
  ####################################################################
  # Extending to the number of observed twins#########################
  ####################################################################
  I_DZ <- Diagonal(N_DZ, 1)
  I_MZ <- Diagonal(N_MZ, 1)
  DZ = lapply(DZ_struc, function(x, I_DZ)kronecker(I_DZ, x), I_DZ = I_DZ)
  MZ = lapply(MZ_struc, function(x, I_MZ)kronecker(I_MZ, x), I_MZ = I_MZ)
  Z_all <- Map(bdiag, DZ, MZ)
  ####################################################################
  ## Extending to multivariate responses #############################
  ####################################################################
  if(n_resp > 1) {
    Z_struc <- mt_struc(n_resp = n_resp)
    twin_A <- lapply(Z_struc, function(x, A)
      kronecker(x, A), A = Z_all$A)
    twin_C <- lapply(Z_struc, function(x, C)
      kronecker(x, C), C = Z_all$C)
    twin_D <- lapply(Z_struc, function(x, D)
      kronecker(x, D), D = Z_all$D)
    twin_E <- lapply(Z_struc, function(x, E)
      kronecker(x, E), E = Z_all$E)
  }
  ####################################################################
  ## Selecting the diferent twin models ##############################
  ####################################################################
  if(n_resp > 1) {
    if(model == "E") {
      output <- twin_E
    }
    if(model == "AE") {
      output <- c(twin_E, twin_A)
    }
    if(model == "CE") {
      output <- c(twin_E, twin_C)
    }
    if(model == "ACE") {
      output <- c(twin_E, twin_A, twin_C)
    }
    if(model == "ADE") {
      output <- c(twin_E, twin_A, twin_D)
    }
  }
  if(n_resp == 1) {
    if(model == "E") {
      output <- list(Z_all$E)
    }
    if(model == "AE") {
      output <- c(Z_all$E, Z_all$A)
    }
    if(model == "CE") {
      output <- c(Z_all$E, Z_all$C)
    }
    if(model == "ACE") {
      output <- c(Z_all$E, Z_all$A, Z_all$C)
    }
    if(model == "ADE") {
      output <- c(Z_all$E, Z_all$A, Z_all$D)
    }
  }
  if(!is.null(formula)) {
    if(length(output) != length(formula)) {
      print("Error: Number of formula does not match number of dispersion components")
    }
    if(length(output) == length(formula)) {
      X_list <- lapply(formula, model.matrix, data = data)
      new_output <- list()
      list_final <- list()
      for(i in 1:length(output)) {
        list_temp <- list()
        for(j in 1:ncol(X_list[[i]])) {
          list_temp[[j]] <- X_list[[i]][,j]*output[[i]]
        }
        list_final[[i]] <- list_temp
      }
      output <- do.call(c,list_final)
    }
}
  return(output)
}


