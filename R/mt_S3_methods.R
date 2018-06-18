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
#' @importFrom stats pnorm deriv
#' @method coef mglm4twin
#' @export

coef.mglm4twin <- function(object, std.error = FALSE, model,
                           response = NULL, ...) {
  n_resp <- length(object$beta_names)
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
  E_main <- paste0("E", 1:n_resp)
  A_main <- paste0("A", 1:n_resp)
  C_main <- paste0("C", 1:n_resp)
  D_main <- paste0("D", 1:n_resp)
  if(n_resp > 1) {
    idx <- combn(n_resp, 2)
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
  }
  if(model == "E") {
    if(n_resp > 1) {
      cod_cov <- c(E_main, E_cross)
    }
    if(n_resp == 1) {
      cod_cov <- E_main
    }
  }
  if(model == "AE") {
    if(n_resp > 1) {
      cod_cov <- c(E_main, E_cross, A_main, A_cross)
    }
    if(n_resp == 1) {
      cod_cov <- c(E_main, A_main)
    }
  }
  if(model == "CE") {
    if(n_resp > 1) {
      cod_cov <- c(E_main, E_cross, C_main, C_cross)
    }
    if(n_resp == 1) {
      cod_cov <- c(E_main, C_main)
    }
  }
  if(model == "ACE") {
    if(n_resp > 1) {
      cod_cov <- c(E_main, E_cross, A_main, A_cross, C_main, C_cross)
    }
    if(n_resp == 1) {
      cod_cov <- c(E_main, A_main, C_main)
    }
  }
  if(model == "ADE") {
    if(n_resp > 1) {
      cod_cov <- c(E_main, E_cross, A_main, A_cross, D_main, D_cross)
    }
    if(n_resp == 1) {
      cod_cov <- c(E_main, A_main, D_main)
    }
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
#' @param biometric Logical. If TRUE biometric measures such as heritability,
#' common environment, genetic correlation etc are printed. (Default = FALSE).
#' @param ... additional arguments affecting the summary produced. Note
#'     the there is no extra options for mglm4twin object class.
#'
#' @return Print a \code{mglm4twin} object.
#' @seealso \code{print}.
#' @method summary mglm4twin
#' @export

summary.mglm4twin <- function(object, model, biometric = FALSE,  ...) {
  n_resp <- length(object$beta_names)
  out <- coef(object, std.error = TRUE, model = model)
  VCOV <- vcov(object, model = model)
  if(n_resp > 1) {
    cross_terms <- combn(n_resp, 2)
  }
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
  output$Dispersion <- out_cov
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
      if(n_resp > 1) {
        output$E_cross <- tab_rho_e
      }
      cat("Environmental correlation:\n")
      print(tab_rho_e)
      cat("\n")
    }
    if(model == "AE") {
      tab_rho_a <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "A", n_resp = n_resp)
      tab_rho_e <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "E", n_resp = n_resp)
      gen <- mt_compute_gen(Estimates = out, vcov = VCOV,
                            model = model, n_resp = n_resp)
      # Environment main effects
      Env_main <- data.frame("Estimates" = diag(gen[[1]]),
                             "std.error" = diag(gen[[3]]))
      Env_main$"z value" <- Env_main[, 1]/Env_main[, 2]
      Env_main$"Pr(>|z|)"  <- 2*pnorm(-abs(Env_main[, 1]/Env_main[, 2]))
      rownames(Env_main) <- paste0("e", 1:n_resp)
      output$Env_main <- Env_main
      # Environment cross effects
      Env_cross <- "Bivariate environmentality no available."
      if(n_resp > 1) {
      Env_cross <- data.frame("Estimates" = gen[[1]][upper.tri(gen[[1]])],
                              "std.error" = gen[[3]][upper.tri(gen[[3]])])
      Env_cross$"z value" <- Env_cross[, 1]/Env_cross[, 2]
      Env_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(Env_cross[, 1]/Env_cross[, 2]))
      rownames(Env_cross) <- paste0("e", cross_terms[1,], cross_terms[2,])
      output$Env_cross <- Env_cross
      }
      # Genetic main effects
      A_main <- data.frame("Estimates" = diag(gen[[2]]),
                             "std.error" = diag(gen[[4]]))
      A_main$"z value" <- A_main[, 1]/A_main[, 2]
      A_main$"Pr(>|z|)"  <- 2*pnorm(-abs(A_main[, 1]/A_main[, 2]))
      rownames(A_main) <- paste0("h", 1:n_resp)
      output$A_main <- A_main
      # Genetic cross effects
      A_cross <- "Bivariate heritability no available."
      if(n_resp > 1) {
      A_cross <- data.frame("Estimates" = gen[[2]][upper.tri(gen[[2]])],
                              "std.error" = gen[[4]][upper.tri(gen[[4]])])
      A_cross$"z value" <- A_cross[, 1]/A_cross[, 2]
      A_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(A_cross[, 1]/A_cross[, 2]))
      rownames(A_cross) <- paste0("h", cross_terms[1,], cross_terms[2,])
      output$A_cross <- A_cross
      }
      cat("Measures associated with genetic structure:\n")
      cat("Genetic correlation:\n")
      print(tab_rho_a)
      cat("\n")
      cat("Heritability:\n")
      print(A_main)
      cat("\n")
      cat("Bivariate heritability:\n")
      print(A_cross)
      cat("\n")
      cat("Measures associated with environment structure:\n")
      cat("Environment correlation:\n")
      print(tab_rho_e)
      cat("\n")
      cat("Environmentality:\n")
      print(Env_main)
      cat("\n")
      cat("Bivariate environmentality:\n")
      print(Env_cross)
      cat("\n")
    }
    if(model == "CE") {
      tab_rho_c <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "C", n_resp = n_resp)
      tab_rho_e <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "E", n_resp = n_resp)
      gen <- mt_compute_gen(Estimates = out, vcov = VCOV,
                            model = model, n_resp = n_resp)
      # Environment main effects
      E_main <- data.frame("Estimates" = diag(gen[[1]]),
                             "std.error" = diag(gen[[3]]))
      E_main$"z value" <- E_main[, 1]/E_main[, 2]
      E_main$"Pr(>|z|)"  <- 2*pnorm(-abs(E_main[, 1]/E_main[, 2]))
      rownames(E_main) <- paste0("e", 1:n_resp)
      output$E_main <- E_main
      # Environment cross effects
      E_cross <- "Bivariate environmentality no available."
      if(n_resp > 1) {
      E_cross <- data.frame("Estimates" = gen[[1]][upper.tri(gen[[1]])],
                            "std.error" = gen[[3]][upper.tri(gen[[3]])])
      E_cross$"z value" <- E_cross[, 1]/E_cross[, 2]
      E_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(E_cross[, 1]/E_cross[, 2]))
      rownames(E_cross) <- paste0("e", cross_terms[1,], cross_terms[2,])
      output$E_cross <- E_cross
      }
      # Common environment main effects
      C_main <- data.frame("Estimates" = diag(gen[[2]]),
                           "std.error" = diag(gen[[4]]))
      C_main$"z value" <- C_main[, 1]/C_main[, 2]
      C_main$"Pr(>|z|)"  <- 2*pnorm(-abs(C_main[, 1]/C_main[, 2]))
      rownames(C_main) <- paste0("c", 1:n_resp)
      output$C_main <- C_main
      # Common environment cross effects
      C_cross <- "Bivariate common environmentality no available."
      if(n_resp > 1) {
      C_cross <- data.frame("Estimates" = gen[[2]][upper.tri(gen[[2]])],
                            "std.error" = gen[[4]][upper.tri(gen[[4]])])
      C_cross$"z value" <- C_cross[, 1]/C_cross[, 2]
      C_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(C_cross[, 1]/C_cross[, 2]))
      rownames(C_cross) <- paste0("c", cross_terms[1,], cross_terms[2,])
      output$C_cross <- C_cross
      }
      cat("Measures associated with common environment structure:\n")
      cat("Common environment correlation:\n")
      print(tab_rho_c)
      cat("\n")
      cat("Common environment:\n")
      print(C_main)
      cat("\n")
      cat("Bivariate common environmentality:\n")
      print(C_cross)
      cat("\n")
      cat("Measures associated with environment structure:\n")
      cat("Environment correlation:\n")
      print(tab_rho_e)
      cat("\n")
      cat("Environmentality:\n")
      print(E_main)
      cat("\n")
      cat("Bivariate environmentality:\n")
      print(E_cross)
      cat("\n")
    }
    if(model == "ACE") {
      tab_rho_a <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "A", n_resp = n_resp)
      tab_rho_c <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "C", n_resp = n_resp)
      tab_rho_e <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "E", n_resp = n_resp)
      gen <- mt_compute_gen(Estimates = out, vcov = VCOV,
                            model = model, n_resp = n_resp)
      # Environment main effects
      E_main <- data.frame("Estimates" = diag(gen[[1]]),
                            "std.error" = diag(gen[[4]]))
      E_main$"z value" <- E_main[, 1]/E_main[, 2]
      E_main$"Pr(>|z|)"  <- 2*pnorm(-abs(E_main[, 1]/E_main[, 2]))
      rownames(E_main) <- paste0("e", 1:n_resp)
      output$E_main <- E_main
      # Environment cross effects
      E_cross <- "Bivariate environmetality no available."
      if(n_resp > 1) {
      E_cross <- data.frame("Estimates" = gen[[1]][upper.tri(gen[[1]])],
                            "std.error" = gen[[4]][upper.tri(gen[[4]])])
      E_cross$"z value" <- E_cross[, 1]/E_cross[, 2]
      E_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(E_cross[, 1]/E_cross[, 2]))
      rownames(E_cross) <- paste0("e", cross_terms[1,], cross_terms[2,])
      output$E_cross <- E_cross
      }
      # Genetic main effects
      A_main <- data.frame("Estimates" = diag(gen[[2]]),
                           "std.error" = diag(gen[[5]]))
      A_main$"z value" <- A_main[, 1]/A_main[, 2]
      A_main$"Pr(>|z|)"  <- 2*pnorm(-abs(A_main[, 1]/A_main[, 2]))
      rownames(A_main) <- paste0("h", 1:n_resp)
      output$A_main <- A_main
      # Genetic cross effects
      A_cross <- "Bivariate heritability no available."
      if(n_resp > 1) {
      A_cross <- data.frame("Estimates" = gen[[2]][upper.tri(gen[[2]])],
                            "std.error" = gen[[5]][upper.tri(gen[[5]])])
      A_cross$"z value" <- A_cross[, 1]/A_cross[, 2]
      A_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(A_cross[, 1]/A_cross[, 2]))
      rownames(A_cross) <- paste0("h", cross_terms[1,], cross_terms[2,])
      output$A_cross <- A_cross
      }
      # Commmon environment main effects
      C_main <- data.frame("Estimates" = diag(gen[[3]]),
                           "std.error" = diag(gen[[6]]))
      C_main$"z value" <- C_main[, 1]/C_main[, 2]
      C_main$"Pr(>|z|)"  <- 2*pnorm(-abs(C_main[, 1]/C_main[, 2]))
      rownames(C_main) <- paste0("c", 1:n_resp)
      output$C_main <- C_main
      # Common environment cross effects
      C_cross <- "Bivariate common environmentality no available."
      if(n_resp > 1) {
      C_cross <- data.frame("Estimates" = gen[[3]][upper.tri(gen[[3]])],
                            "std.error" = gen[[6]][upper.tri(gen[[6]])])
      C_cross$"z value" <- C_cross[, 1]/C_cross[, 2]
      C_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(C_cross[, 1]/C_cross[, 2]))
      rownames(C_cross) <- paste0("c", cross_terms[1,], cross_terms[2,])
      output$C_cross <- C_cross
      }
      cat("Measures associated with genetic structure:\n")
      cat("Genetic correlation:\n")
      print(tab_rho_a)
      cat("\n")
      cat("Heritability:\n")
      print(A_main)
      cat("\n")
      cat("Bivariate heritability:\n")
      print(A_cross)
      cat("\n")
      cat("Measures associated with common environment structure:\n")
      cat("Common environment correlation:\n")
      print(tab_rho_c)
      cat("\n")
      cat("Common environment:\n")
      print(C_main)
      cat("\n")
      cat("Bivariate common environment:\n")
      print(C_cross)
      cat("\n")
      cat("Measures associated with environment structure:\n")
      cat("Environment correlation:\n")
      print(tab_rho_e)
      cat("\n")
      cat("Environmentality:\n")
      print(E_main)
      cat("\n")
      cat("Bivariate environmentality:\n")
      print(E_cross)
      cat("\n")
    }
    if(model == "ADE") {
      tab_rho_a <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "A", n_resp = n_resp)
      tab_rho_d <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "D", n_resp = n_resp)
      tab_rho_e <- mt_compute_rho(Estimates = out, vcov = VCOV,
                                  component = "E", n_resp = n_resp)
      gen <- mt_compute_gen(Estimates = out, vcov = VCOV,
                            model = model, n_resp = n_resp)
      # Environment main effects
      E_main <- data.frame("Estimates" = diag(gen[[1]]),
                           "std.error" = diag(gen[[4]]))
      E_main$"z value" <- E_main[, 1]/E_main[, 2]
      E_main$"Pr(>|z|)"  <- 2*pnorm(-abs(E_main[, 1]/E_main[, 2]))
      rownames(E_main) <- paste0("e", 1:n_resp)
      output$E_main <- E_main
      # Environment cross effects
      E_cross <- "Bivariate environmentality no available."
      if(n_resp > 1) {
      E_cross <- data.frame("Estimates" = gen[[1]][upper.tri(gen[[1]])],
                            "std.error" = gen[[4]][upper.tri(gen[[4]])])
      E_cross$"z value" <- E_cross[, 1]/E_cross[, 2]
      E_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(E_cross[, 1]/E_cross[, 2]))
      rownames(E_cross) <- paste0("e", cross_terms[1,], cross_terms[2,])
      output$E_cross <- E_cross
      }
      # Genetic main effects
      A_main <- data.frame("Estimates" = diag(gen[[2]]),
                           "std.error" = diag(gen[[5]]))
      A_main$"z value" <- A_main[, 1]/A_main[, 2]
      A_main$"Pr(>|z|)"  <- 2*pnorm(-abs(A_main[, 1]/A_main[, 2]))
      rownames(A_main) <- paste0("h", 1:n_resp)
      output$A_main <- A_main
      # Genetic cross effects
      A_cross <- "Bivariate heritability no available."
      if(n_resp > 1) {
      A_cross <- data.frame("Estimates" = gen[[2]][upper.tri(gen[[2]])],
                            "std.error" = gen[[5]][upper.tri(gen[[5]])])
      A_cross$"z value" <- A_cross[, 1]/A_cross[, 2]
      A_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(A_cross[, 1]/A_cross[, 2]))
      rownames(A_cross) <- paste0("h", cross_terms[1,], cross_terms[2,])
      output$A_cross <- A_cross
      }
      # D main effects
      D_main <- data.frame("Estimates" = diag(gen[[3]]),
                           "std.error" = diag(gen[[6]]))
      D_main$"z value" <- D_main[, 1]/D_main[, 2]
      D_main$"Pr(>|z|)"  <- 2*pnorm(-abs(D_main[, 1]/D_main[, 2]))
      rownames(D_main) <- paste0("d", 1:n_resp)
      output$D_main <- D_main
      # D cross effects
      D_cross <- "Bivariate D components no available."
      if(n_resp > 1) {
      D_cross <- data.frame("Estimates" = gen[[3]][upper.tri(gen[[3]])],
                            "std.error" = gen[[6]][upper.tri(gen[[6]])])
      D_cross$"z value" <- D_cross[, 1]/D_cross[, 2]
      D_cross$"Pr(>|z|)"  <- 2*pnorm(-abs(D_cross[, 1]/D_cross[, 2]))
      rownames(D_cross) <- paste0("d", cross_terms[1,], cross_terms[2,])
      output$D_cross <- D_cross
      }
      cat("Measures associated with genetic structure:\n")
      cat("Genetic correlation:\n")
      print(tab_rho_a)
      cat("\n")
      cat("Heritability:\n")
      print(A_main)
      cat("\n")
      cat("Bivariate heritability:\n")
      print(A_cross)
      cat("\n")
      cat("Measures associated with D structure:\n")
      cat("D correlation:\n")
      print(tab_rho_d)
      cat("\n")
      cat("D:\n")
      print(D_main)
      cat("\n")
      cat("Bivariate D:\n")
      print(D_cross)
      cat("\n")
      cat("Measures associated with environment structure:\n")
      cat("Environment correlation:\n")
      print(tab_rho_e)
      cat("\n")
      cat("Environmentality:\n")
      print(E_main)
      cat("\n")
      cat("Bivariate environmentality:\n")
      print(E_cross)
      cat("\n")
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
