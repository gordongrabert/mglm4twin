#' @title Auxiliary function for computing standard errors
#' @name aux_summary
#' @author Wagner Hugo Bonat, \email{wbonat@ufpr.br}
#'
#' @description This function provides alternatives ways to compute standard errors.
#'
#' @param object An object of mglm4twin class.
#' @param formula formula specifying regression model for each dispersion component.
#' @param type Type of standard errors: otimist, robust1 or robust2.
#' @param data Data set.
#' @param ... additional arguments affecting the summary produced.
#' @keywords internal
#' @export
#' @return Internal function returns standard errors.

aux_summary <- function(object, formula, type = "otimist",  data, ...) {
  estimates = c(object$Regression, object$Covariance)
  n_beta <- length(object$Regression)
  n_total <- length(estimates)
  std_error <- sqrt(diag(object$vcov))
  if(type == "robust1") {
    std_error <- sqrt(diag(object$vcov))
  }
  if(type == "otimist") {
    std_error2 <- sqrt(diag(-2*object$joint_inv_sensitivity[c(n_beta+1):n_total,c(n_beta+1):n_total]))
    std_error[c(n_beta+1):n_total] <- std_error2
  }
  if(type == "robust2") {
    MT = mt_emp_variability(object, ...)
    std_error <- sqrt(diag(object$joint_inv_sensitivity%*%MT%*% t(object$joint_inv_sensitivity)))
  }
  Z <- estimates/std_error
  nomes_beta <- do.call(c,lapply(object$mu_list, function(x)colnames(x$D)))
  nomes_cov <- lapply(formula, function(x, data) colnames(model.matrix(x, data = data)), data = data)
  nomes_cov <- do.call(c, nomes_cov)
  nomes <- names(formula)
  IDX <- list()
  dd <- data.frame("Age" = 1)
  for(i in 1:length(nomes)) {
    nn <- colnames(model.matrix(formula[[i]], data = data))
    IDX[[i]] = paste(nomes[i], nn)
  }
  output <- data.frame("Parameter" = do.call(c, IDX),
                       "Estimates" = estimates[c(n_beta+1):n_total],
                       "Std. Error" = std_error[c(n_beta+1):n_total],
                       "Z value" = Z[c(n_beta+1):n_total])
  return(output)
}
