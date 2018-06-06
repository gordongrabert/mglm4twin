#' @title Print
#' @name print.mglm4twin
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description The default print method for an object of \code{mglm4twin} class.
#'
#' @param x fitted model objects of class \code{mglm4twin} as produced
#' by \code{mglm4twin()}.
#' @param ... further arguments passed to or from other methods.
#'
#' @seealso \code{summary}.
#' @rdname print.mglm4twin
#' @method print mglm4twin
#' @export

print.mglm4twin <- function(x, ...) {
  object <- x
  n_resp <- length(object$beta_names)
  regression <- mt_updateBeta(list_initial = list(),
                              betas = object$Regression,
                              information = object$Information,
                              n_resp = n_resp)
  for (i in 1:n_resp) {
    cat("Call: ")
    print(object$linear_pred[[i]])
    cat("\n")
    cat("Link function:", object$link[[i]])
    cat("\n")
    cat("Variance function:", object$variance[[i]])
    cat("\n")
    names(regression[[1]][[i]]) <- object$beta_names[[i]]
    cat("Regression:\n")
    print(regression[[1]][[i]])
    cat("\n")
  }
    cat("Dispersion:\n")
    tau_temp <- object$Covariance
    names(tau_temp) <- rep("", length(tau_temp))
    print(tau_temp)
    cat("\n")
}

#' @title Model Coefficients
#' @name coef.mglm4twin
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Extract model coefficients for objects
#'              of \code{mglm4twin} class.
#'
#' @param object an object of \code{mglm4twin} class.
#' @param std.error logical. If \code{TRUE} returns the standard errors
#'     for the estimates. Default is \code{FALSE}.
#' @param model Type of the fitted model. Options are \code{E}, \code{AE},
#'     \code{CE}, \code{ACE} and \code{ADE}.
#' @param response For which response regression coefficients are extracted.
#'     Default = NULL. It means all estimates are returned including
#'     dispersion estimates for all responses.
#' @param ... additional arguments affecting the summary produced. Note
#'     that there is no extra options for \code{mglm4twin} object class.
#'
#' @return A \code{data.frame} with parameters names, estimates,
#'     response variable number and parameters type.
#'
#' @method coef mglm4twin
#' @export

coef.mglm4twin <- function(object, std.error = FALSE, model,
                           response = NULL, ...) {
  Estimates <- c(object$Regression, object$Covariance)
  cod_beta <- list()
  beta_num <- list()
  power_code <- list()
  for(i in 1:n_resp) {
    cod_beta[[i]] <- rep(i, object$Information$n_betas[[i]])
    beta_num[[i]] <- rep(0:c(object$Information$n_betas[[i]]-1))
    if(object$Information$n_power[[i]] != 0) {
      power_code[[i]] <- rep(i, object$Information$n_power[[i]])
    }
  }
  cod_beta <- unlist(cod_beta)
  beta_nul <- unlist(beta_num)
  beta_nome <- paste0("beta", paste0(cod_beta, beta_nul))
  if(length(power_code) != 0) {
    power_nome <- paste0("power", unlist(power_code))
  }
  # Genetic components
  idx <- combn(n_resp, 2)
  E_main <- paste0("E", 1:n_resp)
  A_main <- paste0("A", 1:n_resp)
  C_main <- paste0("C", 1:n_resp)
  D_main <- paste0("D", 1:n_resp)
  E_cross <- list()
  A_cross <- list()
  C_cross <- list()
  D_cross <- list()
  for(i in 1:dim(idx)[2]) {
    E_cross[[i]] <- paste0("E", paste0(idx[1,i],idx[2,i]))
    A_cross[[i]] <- paste0("A", paste0(idx[1,i],idx[2,i]))
    C_cross[[i]] <- paste0("C", paste0(idx[1,i],idx[2,i]))
    D_cross[[i]] <- paste0("D", paste0(idx[1,i],idx[2,i]))
  }
  E_cross <- unlist(E_cross)
  A_cross <- unlist(A_cross)
  C_cross <- unlist(C_cross)
  D_cross <- unlist(D_cross)
  if(model == "E") {
    cod_cov <- c(E_main, E_cross)
  }
  if(model == "AE") {
    cod_cov <- c(E_main, E_cross, A_main, A_cross)
  }
  if(model == "CE") {
    cod_cov <- c(E_main, E_cross, C_main, C_cross)
  }
  if(model == "ACE") {
    cod_cov <- c(E_main, E_cross, A_main, A_cross, C_main, C_cross)
  }
  if(model == "ADE") {
    cod_cov <- c(E_main, E_cross, A_main, A_cross, D_main, D_cross)
  }
  if(length(power_code) != 0) {
    Parameters <- c(beta_nome, power_nome, cod_cov)
  }
  if(length(power_code) == 0) {
    Parameters <- c(beta_nome, cod_cov)
  }
  coef_temp <- data.frame(Parameters = Parameters,
                          Estimates = Estimates)
  if (std.error == TRUE) {
      coef_temp$std.error = sqrt(diag(object$vcov))
  }
   if(!is.null(response)) {
     coef_temp$Response <- c(cod_beta, rep(NA, object$Information$n_cov))
     coef_temp <- coef_temp[which(coef_temp$Response %in% response), ]
   }
  return(coef_temp)
}

#' @title Summarizing
#' @name summary.mglm4twin
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description The default summary method for an object of \code{mglm4twin} class.
#'
#' @param object an object of \code{mglm4twin} class.
#' @param model String indicating twin model.
#' @biometric Logical. If TRUE biometric measures such as heritability,
#' common environment, genetic correlation etc are printed. (Default = FALSE).
#' @param ... additional arguments affecting the summary produced. Note
#'     the there is no extra options for mglm4twin object class.
#'
#' @return Print a \code{mglm4twin} object.
#' @seealso \code{print}.
#' @method summary mglm4twin
#' @export

summary.mglm4twin <- function(object, model, biometric = FALSE,  ...) {
  out <- coef(object, std.error = TRUE, model = model)
  VCOV <- vcov(object, model = model)
  n_resp <- length(object$beta_names)
  output <- list()
  for(i in 1:n_resp) {
    tab_beta <- coef(object, std.error = TRUE, model = model,
                     response = i, type = "beta")[, 2:3]
    tab_beta$"z value" <- tab_beta[, 1]/tab_beta[, 2]
    tab_beta$"Pr(>|z|)"  <- 2*pnorm(-abs(tab_beta[, 1]/tab_beta[, 2]))
    rownames(tab_beta) <- object$beta_names[[i]]
    output[i][[1]]$Regression <- tab_beta
  }
  out$"z value" <- out[, 2]/out[, 3]
  out$"Pr(>|z|)"  <- 2*pnorm(-abs(out[, 2]/out[, 3]))
  for (i in 1:n_resp) {
    cat("Call: ")
    print(object$linear_pred[[i]])
    cat("\n")
    cat("Link function:", object$link[[i]])
    cat("\n")
    cat("Variance function:", object$variance[[i]])
    cat("\n")
    cat("Regression:\n")
    print(output[i][[1]]$Regression)
    cat("\n")
    if(object$Information$n_power[[i]] != 0) {
      cat("Power:\n")
      tt = paste0("power",i)
      rownames(out)[as.numeric(rownames(out[which(out$Parameters == tt),-1]))] <- tt
      print(out[which(out$Parameters == tt),-1])
      cat("\n")
    }
  }
  n_beta <- sum(unlist(object$Information$n_betas))
  n_power <- sum(unlist(object$Information$n_power))
  n_cov <- object$Information$n_taus
  n_comp <- n_resp + n_resp*(n_resp-1)/2
  n_total <- dim(out)[1]
  out_cov <- out[c(n_beta+n_power+1):n_total,]
  rownames(out_cov) <- out_cov$Parameters
  out_cov$"z value" <- out_cov[, 2]/out_cov[, 3]
  out_cov$"Pr(>|z|)"  <- 2*pnorm(-abs(out_cov[, 2]/out_cov[, 3]))
  cat("Dispersion:\n")
  if(model == "E") {
    cat("Environment component:\n")
    print(out_cov[1:c(n_comp),-1])
    cat("\n")
  }
  if(model == "AE") {
    cat("Genetic component:\n")
    print(out_cov[c(n_comp+1):c(2*n_comp),-1])
    cat("\n")
    cat("Environment component:\n")
    print(out_cov[1:c(n_comp),-1])
    cat("\n")
  }
  if(model == "CE") {
    cat("Common environment component:\n")
    print(out_cov[c(n_comp+1):c(2*n_comp),-1])
    cat("\n")
    cat("Environment component:\n")
    print(out_cov[1:c(n_comp),-1])
    cat("\n")
  }
  if(model == "ACE") {
    cat("Genetic component:\n")
    print(out_cov[c(n_comp+1):c(2*n_comp),-1])
    cat("\n")
    cat("Common environment component:\n")
    print(out_cov[c(2*n_comp+1):c(3*n_comp),-1])
    cat("\n")
    cat("Environment component:\n")
    print(out_cov[1:c(n_comp),-1])
    cat("\n")
  }
  if(model == "ADE") {
    cat("Genetic component:\n")
    print(out_cov[c(n_comp+1):c(2*n_comp),-1])
    cat("\n")
    cat("D component:\n")
    print(out_cov[c(2*n_comp+1):c(3*n_comp),-1])
    cat("\n")
    cat("Environment component:\n")
    print(out_cov[1:c(n_comp),-1])
    cat("\n")
  }
  if(biometric == TRUE) {
    if(model == "E") {
      tab_rho_e <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "E", n_resp = n_resp)
      cat("Environmental correlation:\n")
      print(tab_rho_e)
      cat("\n")
    }
    if(model == "AE") {
      tab_rho_a <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "A", n_resp = n_resp)
      tab_rho_e <- mt_compute_rho(Estimates = out_cov, vcov = VCOV,
                                  component = "E", n_resp = n_resp)


    }

  }


  names(object$con$correct) <- ""
  iteration_cov <- length(na.exclude(object$IterationCovariance[, 1]))
  names(iteration_cov) <- ""
  names(object$con$method) <- ""
  cat("Algorithm:", object$con$method)
  cat("\n")
  cat("Correction:", object$con$correct)
  cat("\n")
  cat("Number iterations:", iteration_cov)
  cat("\n")
  return(invisible(output))
}

#' @title Variance-Covariance Matrix
#' @name vcov.mglm4twin
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Returns the variance-covariance matrix for an object
#' of \code{mglm4twin} class.
#'
#' @param object an object of \code{mglm4twin} class.
#' @param model String indicating twin model.
#' @param ... additional arguments affecting the summary produced. Note
#'     that there is no extra options for mcglm object class.
#'
#' @return A variance-covariance matrix.
#'
#' @method vcov mglm4twin
#' @export

vcov.mglm4twin <- function(object, model, ...) {
  cod <- coef(object, model = model)$Parameters
  colnames(object$vcov) <- cod
  rownames(object$vcov) <- cod
  return(object$vcov)
}
