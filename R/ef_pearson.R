#' @title Pearson estimating function
#' @author Wagner Hugo Bonat
#'
#' @description Compute the Pearson estimating function its sensitivity
#'     and variability matrices.
#'
#' @param y_vec A vector.
#' @param mu_vec A vector.
#' @param Cfeatures A list of matrices.
#' @param inv_J_beta A matrix.
#' @param D A matrix.
#' @param correct Logical.
#' @param compute_variability Logical.
#' @param compute_sensitivity Logical.
#' @return A list with three components: (i) a vector of quasi-score
#'     values, (ii) the sensitivity and (iii) variability matrices
#'     associated with the Pearson estimating function.
#' @keywords internal
#' @details Compute the Pearson estimating function its sensitivity and
#'     variability matrices.  For more details see Bonat and Jorgensen
#'     (2016) equations 6, 7 and 8.

ef_pearson <- function(y_vec, mu_vec, Cfeatures, inv_J_beta = NULL,
                       D = NULL, correct = FALSE,
                       compute_variability = FALSE) {
  product <- lapply(Cfeatures$D_Sigma, ef_multiply,
                    bord2 = Cfeatures$inv_Sigma)
  res <- y_vec - mu_vec
  pearson_score <- unlist(lapply(product, ef_core_pearson,
                                 inv_C = Cfeatures$inv_Sigma, res = res))
  sensitivity <- ef_sensitivity(product)
  output <- list(Score = pearson_score, Sensitivity = sensitivity,
                 Extra = product)
  if (correct == TRUE) {
    correction <- ef_correction(D_C = Cfeatures$D_Sigma,
                                inv_J_beta = inv_J_beta, D = D,
                                inv_C = Cfeatures$inv_Sigma)
    output <- list(Score = pearson_score - correction,
                   Sensitivity = sensitivity, Extra = product)
  }
  if (compute_variability == TRUE) {
    variability <- ef_variability(sensitivity = sensitivity,
                                  product = product,
                                  inv_C = Cfeatures$inv_Sigma,
                                  C = Cfeatures$Sigma, res = res)
    output$Variability <- variability
  }
  return(output)
}
