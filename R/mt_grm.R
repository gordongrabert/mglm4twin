#' @title Build genetic relationship matrix linear predictor components
#' @name mt_grm
#' @author Gordon Grabert, \email{wbonat@ufpr.br}
#'
#' @description This function provides the components of the matrix
#' linear predictor suitable for fitting multivariate genome-based models
#' by using a precomputed genetic relationship matrix.
#'
#' @param n Number of observations.
#' @param grm Genetic relationship matrix.
#' @param n_resp Number of response variables.
#' @param model a string specifying the name of the grm model.
#' Options are: \code{"E"}, \code{"AE"}
#' @param formula formula specifying regression model for each dispersion component.
#' @param data Data set.
#' @keywords internal
#' @export
#' @return A list of matrices.


mt_grm <- function(n, grm, n_resp, model, formula = NULL, data = NULL) {

  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################

  grm_sparse <- Matrix(grm, sparse = F)
  A <- as(grm_sparse, "dsCMatrix")
  E <- diag(n)
  # n = n
  # model = "AE"
  # formula = list("formE1" = ~ 1, "formE2" = ~ 1, "formE12" = ~ 1,
  #                "formA1" = ~ 1, "formA2" = ~ 1, "formA12" = ~ 1)
  # data = pheno

  ####################################################################
  ## Extending to multivariate responses #############################
  ####################################################################
  if(n_resp > 1) {
    Z_struc <-  mglm4twin:::mt_struc(n_resp = n_resp)
    ind_A <- lapply(Z_struc, function(x, A)
      kronecker(x, A), A = A)
    ind_E <- lapply(Z_struc, function(x, E)
      kronecker(x, E), E = E)
  }
  ####################################################################
  ## Selecting the diferent twin models ##############################
  ####################################################################
  if(n_resp > 1) {
    if(model == "AE") {
      output <- c(ind_A, ind_E) # TODO changed ordering of A and E matrix compared to mt_twin
    }
  }
  if(n_resp == 1) {
    if(model == "E") {
      output <- list(ind_E)
    }
    if(model == "AE") {
      output <- c(ind_A, ind_E)
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
