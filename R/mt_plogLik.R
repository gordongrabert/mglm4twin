#' @title Gaussian Pseudo-loglikelihood
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Extract the Gaussian pseudo-loglikelihood (plogLik)
#' value for objects of \code{mglm4twin} class.
#'
#' @param object an object or a list of objects representing a model
#' of \code{mglm4twin} class.
#' @param verbose logical. Print or not the plogLik value.
#' @return Returns the value of the Gaussian pseudo-loglikelihood.
#' @export

plogLik <- function(object, verbose = TRUE) {
    my.gauss <- function(b, det.Sigma, inv.Sigma){
        n <- length(b)
        dens <- (-n/2)*log(2*pi) - 0.5*det.Sigma -
            0.5*t(b)%*%inv.Sigma%*%b
        return(dens)
    }
    if(class(object) == "mglm4twin") {
    b <- c(as.matrix(object$observed)) - c(as.matrix(object$fitted))
    det.Sigma <- determinant(object$C)$modulus
    ll <- round(as.numeric(my.gauss(b = b, det.Sigma = det.Sigma,
                              inv.Sigma = object$inv_C)),2)
    df <- length(object$Regression) + length(object$Covariance)
    if (verbose) cat("Pseudo log Lik.", ll, sprintf("(df=%d)",df))
    return(invisible(list("plogLik" = ll, "df" = df)))
    }
    if(class(object) == "list") {
      Y <- do.call(c,lapply(object, function(x)as.numeric(x$observed)))
      mu <-do.call(c,lapply(object, function(x)as.numeric(x$fitted)))
      b <- Y - mu
      C.list <- lapply(object, function(x)x$C)
      C <- bdiag(C.list)
      inv_C.list <- lapply(object, function(x)x$inv_C)
      inv_C <- bdiag(inv_C.list)
      det.Sigma <- determinant(C)$modulus
      ll <- round(as.numeric(my.gauss(b = b, det.Sigma = det.Sigma,
                                      inv.Sigma = inv_C)),2)
      #df <- length(object$Regression) + length(object$Covariance)
      df <- sum(do.call(c,lapply(object, function(x)as.numeric(length(x$Regression) + length(x$Covariance)))))

      if (verbose) cat("Pseudo log Lik.", ll, sprintf("(df=%d)",df))
      return(invisible(list("plogLik" = ll, "df" = df)))
    }
}

