#' @title Compute genetic measures and their standard errors
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Compute genetic measures and their standard errors
#' using the delta method.
#'
#' @param Estimates Table of estimates, standard errors and parameter names.
#' @param vcov Matrix of variance and covariance.
#' @param model String.
#' @param n_resp Numeric. Number of response variables.
#' @keywords internal
#' @details It is an internal function useful in general for summary
#' function associated with Twin models.

mt_compute_gen <- function(Estimates, vcov, model, n_resp) {
  AA <- matrix(NA, ncol = n_resp, nrow = n_resp)
  AA_std <- matrix(NA, ncol = n_resp, nrow = n_resp)
  EE <- matrix(NA, ncol = n_resp, nrow = n_resp)
  EE_std <- matrix(NA, ncol = n_resp, nrow = n_resp)
  CC <- matrix(NA, ncol = n_resp, nrow = n_resp)
  CC_std <- matrix(NA, ncol = n_resp, nrow = n_resp)
  DD <- matrix(NA, ncol = n_resp, nrow = n_resp)
  DD_std <- matrix(NA, ncol = n_resp, nrow = n_resp)
  for(i in 1:n_resp) {
    for(j in i:n_resp) {
      #print(c(i,j))
      if(model == "ACE") {
        if(i == j) {
          Param <- paste0(c("E","A", "C"), i)
        }
        if(i != j) {
          Param <- paste0(c("E","A","C"), i, j)
        }
        Point <- Estimates[which(Estimates$Parameters %in% Param),]$Estimates
        VCOV_temp <- vcov[Param,Param]
        EE[i,j] <- Point[1]/(Point[1] + Point[2] + Point[3])
        AA[i,j] <- Point[2]/(Point[1] + Point[2] + Point[3])
        CC[i,j] <- Point[3]/(Point[1] + Point[2] + Point[3])
        EE_std[i,j] <- mt_delta_method(fx = ~ x1/(x1+x2+x3),
                                       point = Point, cov = VCOV_temp)
        AA_std[i,j] <- mt_delta_method(fx = ~ x2/(x1+x2+x3),
                                       point = Point, cov = VCOV_temp)
        CC_std[i,j] <- mt_delta_method(fx = ~ x3/(x1+x2+x3),
                                       point = Point, cov = VCOV_temp)
        output <- list(EE, AA, CC, EE_std, AA_std, CC_std)
      }
      if(model == "ADE") {
        if(i == j) {
          Param <- paste0(c("E","A", "D"), i)
        }
        if(i != j) {
          Param <- paste0(c("E","A","D"), i, j)
        }
        Point <- Estimates[which(Estimates$Parameters %in% Param),]$Estimates
        VCOV_temp <- vcov[Param,Param]
        EE[i,j] <- Point[1]/(Point[1] + Point[2] + Point[3])
        AA[i,j] <- Point[2]/(Point[1] + Point[2] + Point[3])
        DD[i,j] <- Point[3]/(Point[1] + Point[2] + Point[3])
        EE_std[i,j] <- mt_delta_method(fx = ~ x1/(x1+x2+x3),
                                         point = Point, cov = VCOV_temp)
        AA_std[i,j] <- mt_delta_method(fx = ~ x2/(x1+x2+x3),
                                         point = Point, cov = VCOV_temp)
        DD_std[i,j] <- mt_delta_method(fx = ~ x3/(x1+x2+x3),
                                         point = Point, cov = VCOV_temp)
        output <- list(EE, AA, DD, EE_std, AA_std, DD_std)
      }
      if(model == "AE") {
        if(i == j) {
          Param <- paste0(c("E","A"), i)
        }
        if(i != j) {
          Param <- paste0(c("E","A"), i, j)
        }
        Point <- Estimates[which(Estimates$Parameters %in% Param),]$Estimates
        VCOV_temp <- vcov[Param,Param]
        EE[i,j] <- Point[1]/(Point[1] + Point[2])
        EE_std[i,j] <- mt_delta_method(fx = ~ x1/(x1+x2),
                                       point = Point, cov = VCOV_temp)
        AA[i,j] <- Point[2]/(Point[1] + Point[2])
        AA_std[i,j] <- mt_delta_method(fx = ~ x2/(x1+x2),
                                       point = Point, cov = VCOV_temp)
        output <- list(EE, AA, EE_std, AA_std)
      }
      if(model == "CE") {
        if(i == j) {
          Param <- paste0(c("E","C"), i)
        }
        if(i != j) {
          Param <- paste0(c("E","C"), i, j)
        }
        Point <- Estimates[which(Estimates$Parameters %in% Param),]$Estimates
        VCOV_temp <- vcov[Param,Param]
        EE[i,j] <- Point[1]/(Point[1] + Point[2])
        EE_std[i,j] <- mt_delta_method(fx = ~ x1/(x1+x2),
                                       point = Point, cov = VCOV_temp)
        CC[i,j] <- Point[2]/(Point[1] + Point[2])
        CC_std[i,j] <- mt_delta_method(fx = ~ x2/(x1+x2),
                                       point = Point, cov = VCOV_temp)
        output <- list(EE, CC, EE_std, CC_std)
      }
    }
  }
  return(output)
}
