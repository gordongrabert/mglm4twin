#' @title Compute correlation and its standard errors
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Compute correlations and standard errors using the delta method.
#'
#' @param Estimates Table of estimates, standard errors and parameter names.
#' @param vcov Matrix of variance and covariance.
#' @param component Logical indicating for which component the correlation
#' will be computed.
#' @param n_resp Numeric. Number of response variables.
#' @keywords internal
#' @details It is an internal function useful in general for summary
#' function associated with Twin models.

mt_compute_rho <- function(Estimates, vcov, component, n_resp) {
  if(n_resp == 1) {
    tab_rho <- "Univariate model no correlation available"
  }
  if(n_resp > 1) {
  idx <- combn(n_resp,2)
  n_cov <- dim(idx)[2]
  rho <- list()
  rho_std <- list()
  for(i in 1:n_cov) {
    E1_temp <- paste0(component, idx[1,i])
    E2_temp <- paste0(component, idx[2,i])
    E12_temp <- paste0(component, paste0(idx[1,i],idx[2,i]))
    Point <- Estimates[which(Estimates$Parameters %in% c(E1_temp, E2_temp, E12_temp)),]$Estimates
    VCOV_temp <- vcov[c(E1_temp, E2_temp, E12_temp),c(E1_temp, E2_temp, E12_temp)]
    rho[[i]] <- Point[3]/sqrt(Point[1]*Point[2])
    rho_std[[i]] <- mt_delta_method(fx = ~ x3/sqrt(x1*x2), point = Point,
                                    cov = VCOV_temp)
  }
  rho_names <- paste0("rho", paste0(idx[1,],idx[2,]))
  tab_rho <- data.frame("Estimates" = unlist(rho),
                          "std.error" = unlist(rho_std))
  tab_rho$"z value" <- tab_rho[, 1]/tab_rho[, 2]
  tab_rho$"Pr(>|z|)"  <- 2*pnorm(-abs(tab_rho[, 1]/tab_rho[, 2]))
  rownames(tab_rho) <- rho_names
  }
  return(tab_rho)
}

