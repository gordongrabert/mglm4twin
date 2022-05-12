#' @title Wald test for dispersion components
#' @name mt_anova_mglm
#' @author Wagner Hugo Bonat, \email{wbonat@ufpr.br}
#'
#' @description This function provides Wald test for dispersion components.
#'
#' @param object An object of mglm4twin class.
#' @param formula formula specifying regression model for each dispersion component.
#' @param data Data set.
#' @param ... additional arguments affecting the summary produced.
#' @keywords internal
#' @export
#' @return Internal function returns Wald statistics.

mt_anova_mglm <- function(object, formula, data, ...) {
  list_Z <- lapply(formula, model.matrix, data = object$data)
  n_idx <- list()
  n_resp = length(list_Z)
  for(i in 1:n_resp) {
    n_idx[[i]] <- rep(i, dim(list_Z[[i]])[2])
  }
  n_beta <- length(object$Regression)
  n_par <- length(c(object$Regression, object$Covariance))
  vv <- object$vcov[c(n_beta+1):n_par,c(n_beta+1):n_par]
  idx.vec <- do.call(c, n_idx)
  temp.vcov <- list()
  temp.cov <- list()
  for (i in 1:n_resp) {
    idx.id <- idx.vec == i
    temp.vcov[[i]] <- vv[idx.id, idx.id]
    temp.cov[[i]] <- object$Covariance[idx.id]
  }
  saida <- list()
  names <- names(formula)
  X2.resp <- list()
  for (j in 1:n_resp) {
    temp <- as.numeric(
      t(temp.cov[[j]]) %*%
        solve(as.matrix(temp.vcov[[j]])) %*%
        temp.cov[[j]])
    nbeta.test <- length(temp.cov[[j]])
    X2.resp[[j]] <-
      data.frame(Covariate = names[j],
                 Chi.Square = round(temp, 4), Df = nbeta.test,
                 p.value = round(pchisq(temp, nbeta.test,
                                        lower.tail = FALSE), 8))
  }
  saida <- do.call(rbind, X2.resp)
  cat("Wald test for dispersion components\n")
  cat("Call: ")
  print(formula)
  cat("\n")
  print(saida)
  cat("\n")
  return(invisible(saida))
}
