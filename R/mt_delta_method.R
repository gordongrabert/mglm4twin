#' @title Delta Method
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Compute standard errors for functions of model parameters
#' using the delta method.
#'
#' @param fx A function (string like formula) of model parameters.
#' @param point Vector of parameter estimates.
#' @param cov Variance-covariance matrix.
#' @param std Logical. If TRUE returns the standard error, otherwise
#' return the new variance-covariance.
#' @keywords internal
#' @details It is an internal function useful in general for summary
#' function associated with Twin models.
#' @export

mt_delta_method <- function (fx, point, cov, std = TRUE) {
    cov <- as.matrix(cov)
    n <- length(point)
    if (!is.list(fx))
        fx <- list(fx)
    if ((dim(cov)[1] != n) || (dim(cov)[2] != n))
        stop(paste("Covariances should be a ", n, " by ", n,
                   " matrix"))
    syms <- paste("x", 1:n, sep = "")
    for (i in 1:n) assign(syms[i], point[i])
    gdashmu <- t(sapply(fx, function(form) {
        as.numeric(attr(eval(deriv(form, syms)), "gradient"))
    }))
    new.covar <- gdashmu %*% cov %*% t(gdashmu)
    if (std) {
        new.se <- sqrt(diag(new.covar))
        new.se
    }
    else new.covar
}
