#' @title Empirical variability matrix
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Compute an empirical version of the variability matrix
#' used to estimate robust standard errors.
#'
#' @param object an object or a list of objects representing a model
#' of \code{mglm4twin} class.
#' @param id Index indicating independent unit samples.
#' @return Returns a matrix.
#'
#' @source Bonat, W. H. (2018). Multiple Response Variables Regression
#' Models in R: The mcglm Package. Journal of Statistical Software, 84(4):1--30.
#' @importFrom stats cov
#' @export

mt_emp_variability <- function(object, id) {
  Cfeatures <- mt_build_sigma(mu = object$mu, tau = object$Covariance,
                              power = do.call(c, object$list_initial$power),
                              Z =  object$matrix_pred,
                              variance = do.call(c,object$variance),
                              Ntrial = object$Ntrial,
                              power_fixed = object$power_fixed,
                              inverse = TRUE,
                              compute_derivative_beta = FALSE)
  y_vec <- as.vector(object$observed)
  mu_vec <- object$fitted
  temp = data.frame("Y" = y_vec, "mu" = mu_vec, "ID" = as.factor(id),
             "idx" = 1:length(y_vec))
  temp_group <- split(temp, temp$ID)
  nrow <- length(temp_group)
  n_par <- length(object$Regression) + length(object$Covariance)
  mat <- matrix(NA, ncol = n_par, nrow = nrow)
  D <- bdiag(lapply(object$mu, function(x) x$D))
  for(i in 1:nrow) {
    Cfeatures_group <- list()
    idx = temp_group[[i]]$idx
    Cfeatures_group$Sigma <- Cfeatures$Sigma[idx,idx]
    Cfeatures_group$inv_Sigma <- Cfeatures$inv_Sigma[idx,idx]
    Cfeatures_group$D_Sigma <- lapply(Cfeatures$D_Sigma,
                                      function(x, idx) x[idx,idx],
                                      idx = idx)
    qscore <- ef_quasi_score(D = D[idx,], Cfeatures_group$inv_Sigma,
                             temp_group[[i]]$Y, mu_vec = temp_group[[i]]$mu)$Score

    score_pearson <- ef_pearson(y_vec = temp_group[[i]]$Y,
                          mu_vec = temp_group[[i]]$mu,
                          Cfeatures = Cfeatures_group)$Score
    mat[i,] <- c(as.numeric(qscore), score_pearson)
  }
  return(cov(mat)*nrow)
}
