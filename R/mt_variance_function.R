#' @title Variance Functions
#' @author Wagner Hugo Bonat, \email{wbonat@@ufpr.br}
#'
#' @description Compute the variance function and its derivatives with
#' respect to regression, dispersion and power parameters.
#'
#' @param mu a numeric vector. In general the output from
#'     \code{\link{mt_link_function}}.
#' @param power a numeric value (\code{tweedie} and \code{binomialP}) or
#'     a vector (\code{binomialPQ}) of the power parameters.
#' @param Ntrial number of trials, useful only when dealing with
#'     binomial response variables.
#' @param variance a string specifying the name (\code{constant, tweedie, binomialP
#'     or binomialPQ}) of the variance function.
#' @param derivative_power logical if compute (TRUE) or not (FALSE) the
#'     derivatives with respect to the power parameter.
#' @param derivative_mu logical if compute (TRUE) or not (FALSE) the
#'     derivative with respect to the mu parameter.
#' @return A list with from one to four elements depending on the
#'     arguments.
#'
#' @seealso \code{\link{mt_link_function}}.
#'
#' @details The function \code{mt_variance_function} computes three
#'     features related with the variance function. Depending on the
#'     logical arguments, the function returns \eqn{V^{1/2}} and its
#'     derivatives with respect to the parameters power and mu,
#'     respectivelly.  The output is a named list, completely
#'     informative about what the function has been computed.  For
#'     example, if \code{derivative_power = TRUE} and \code{derivative_mu = TRUE}.
#'     The output will be a list, with three elements: V_sqrt, D_V_sqrt_power and
#'     D_V_sqrt_mu.
#'
#' @usage mt_variance_function(mu, power, Ntrial, variance,
#'                            derivative_power, derivative_mu)
#'
#' @source Bonat, W. H. and Jorgensen, B. (2016) Multivariate
#'     covariance generalized linear models.
#'     Journal of Royal Statistical Society - Series C 65:649--675.
#'
#' @export
#'
#' @examples
#' x1 <- seq(-1, 1, l = 5)
#' X <- model.matrix(~x1)
#' mu <- mt_link_function(beta = c(1, 0.5), X = X, offset = NULL,
#'                        link = "logit")
#' mt_variance_function(mu = mu$mu, power = c(2, 1), Ntrial = 1,
#'                      variance = "binomialPQ",
#'                      derivative_power = TRUE, derivative_mu = TRUE)
#'

## Generic variance function -------------------------------------------
mt_variance_function <- function(mu, power, Ntrial, variance,
                                 derivative_power,
                                 derivative_mu) {
    assert_that(is.logical(derivative_power))
    assert_that(is.logical(derivative_mu))
    switch(variance,
           constant = {
             output <- mt_constant(mu = mu,
                                   derivative_power = derivative_power,
                                   derivative_mu = derivative_mu)
           },
           tweedie = {
               output <- mt_tweedie(mu = mu, power = power,
                                    derivative_power = derivative_power,
                                    derivative_mu = derivative_mu)
           },
           binomialP = {
               output <- mt_binomialP(mu = mu, power = power,
                                      Ntrial = Ntrial,
                                      derivative_power = derivative_power,
                                      derivative_mu = derivative_mu)
           },
           binomialPQ = {
               output <- mt_binomialPQ(mu = mu, power = power,
                                       Ntrial = Ntrial,
                                       derivative_power = derivative_power,
                                       derivative_mu = derivative_mu)
           },
           stop(gettextf("%s variance function not recognised",
                         sQuote(variance)), domain = NA))
    return(output)
}

#' @rdname mt_variance_function
#' @usage mt_tweedie(mu, power, Ntrial, derivative_power, derivative_mu)
## Tweedie variance function ---------------------------------------------
mt_tweedie <- function(mu, power, Ntrial,
                       derivative_power, derivative_mu) {
    ## The observed value can be zero, but not the expected value.
    assert_that(all(mu > 0))
    assert_that(is.number(power))
    mu.power <- mu^power
    sqrt.mu.power <- sqrt(mu.power)
    n <- length(mu)
    output <- list()
    output$V_sqrt = Diagonal(n = n, sqrt.mu.power)
    if(derivative_power == TRUE) {
      output$D_V_sqrt_power = Diagonal(n = n, (mu.power * log(mu))/(2 * sqrt.mu.power))
    }
    if(derivative_mu == TRUE) {
      output$D_V_sqrt_mu = (mu^(power -  1) * power)/(2 * sqrt.mu.power)
    }
    return(output)
}

#' @rdname mt_variance_function
#' @usage mt_binomialP(mu, power, Ntrial,
#'                     derivative_power, derivative_mu)
## BinomialP variance function
## -----------------------------------------
mt_binomialP <- function(mu, power, Ntrial,
                         derivative_power, derivative_mu) {
    ## The observed value can be 0 and 1, but not the expected value
    assert_that(all(mu > 0))
    assert_that(all(mu < 1))
    assert_that(is.number(power))
    assert_that(all(Ntrial > 0))
    constant <- (1/Ntrial)
    mu.power <- mu^power
    mu.power1 <- (1 - mu)^power
    mu1mu <- constant * (mu.power * mu.power1)
    sqrt.mu1mu <- sqrt(mu1mu)
    n <- length(mu)
    output <- list()
    output$V_sqrt = Diagonal(n = n, sqrt.mu1mu)
    if(derivative_power == TRUE) {
      output$D_V_sqrt_power = Diagonal(n = n, (log(1 - mu) * mu1mu +
                                                 log(mu) * mu1mu)/(2 * sqrt.mu1mu))
    }
    if(derivative_mu == TRUE) {
      output$D_V_sqrt_mu = (constant *
                       (mu.power1 * (mu^(power - 1)) * power) -
                       constant * (((1 - mu)^(power - 1)) *
                                     mu.power * power))/(2 * sqrt.mu1mu)
    }
    return(output)
}
#' @rdname mt_variance_function
#' @usage mt_binomialPQ(mu, power, Ntrial,
#'                      derivative_power, derivative_mu)
## BinomialPQ variance function ----------------------------------------
mt_binomialPQ <- function(mu, power, Ntrial,
                          derivative_power, derivative_mu) {
    ## The observed value can be 0 and 1, but not the expected value
    assert_that(all(mu > 0))
    assert_that(all(mu < 1))
    assert_that(length(power) == 2)
    assert_that(all(Ntrial > 0))
    constant <- (1/Ntrial)
    p <- power[1]
    q <- power[2]
    mu.p <- mu^p
    mu1.q <- (1 - mu)^q
    mu.p.mu.q <- mu.p * mu1.q
    mu1mu <- mu.p.mu.q * constant
    sqrt.mu1mu <- sqrt(mu1mu)
    n <- length(mu)
    denominator <- 2 * sqrt(mu1mu)
    output <- list()
    output$V_sqrt = Diagonal(n = n, sqrt.mu1mu)
    if(derivative_power == TRUE) {
      output$D_V_sqrt_p = Diagonal(n = n,+(mu.p.mu.q * log(mu))/denominator)
      output$D_V_sqrt_q = Diagonal(n = n,+(mu.p.mu.q * log(1 - mu))/denominator)
    }
    if(derivative_mu == TRUE) {
      output$D_V_inv_sqrt_mu = -(constant * (mu1.q * (mu^(p - 1)) * p) -
                            constant * (((1 - mu)^(q - 1)) * mu.p * q))/(2 * (mu1mu^1.5))
    }
    return(output)
}

#' @rdname mt_variance_function
#' @usage mt_constant(mu, power, Ntrial, derivative_power, derivative_mu)
## Constant variance function ----------------------------------------
mt_constant <- function(mu, power, Ntrial,
                        derivative_power, derivative_mu) {
  ## The observed value can be zero, but not the expected value.
  n <- length(mu)
  output <- list()
  output$V_sqrt = Diagonal(n = n, 1)
  if(derivative_power == TRUE) {
    output$D_V_sqrt_power = Diagonal(n = n, 0)
  }
  if(derivative_mu == TRUE) {
    output$D_V_sqrt_mu = rep(0, n)
  }
  return(output)
}
