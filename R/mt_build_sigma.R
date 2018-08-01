#' @title Build variance-covariance matrix
#' @author Wagner Hugo Bonat,  \email{wbonat@@ufpr.br}
#'
#' @description This function builds a variance-covariance matrix, based
#'     on the variance function and Omega matrix.
#'
#'@param mu List with expected values and derivatives. In general the output from
#'     \code{\link{mt_link_function}}.
#'@param tau A numeric vector.
#'@param power A list of numeric or numeric vector. It should be one number for
#'     all variance functions except binomialPQ, in that case the
#'     argument specifies both p and q. Each slot of the list should specify the
#'     power parameter for each response variable.
#'@param Z A list of matrices.
#'@param Ntrial Number of trials for binomial data.
#'@param variance String specifing the variance function: constant,
#'     tweedie, poisson_tweedie, binomialP or binomialPQ.
#'@param power_fixed Logical if the power parameter is fixed at initial
#'     value (TRUE). In the case power_fixed = FALSE the power parameter
#'     will be estimated.
#'@param compute_derivative_beta Logical. Compute or not the derivative
#'     with respect to regression parameters.
#'@keywords internal
#'@return A list with \eqn{\Sigma},
#'     \eqn{\Sigma^{-1}} and the derivative of \eqn{\Sigma} with respect
#'     to the power, tau and beta parameters.
#'@import Matrix
#'@seealso \code{\link{mt_link_function}},
#'     \code{\link{mt_variance_function}}, \code{\link{mt_build_omega}}.

mt_build_sigma <- function(mu, tau, power, Z, Ntrial,
                           variance, power_fixed, inverse = FALSE,
                           compute_derivative_beta = TRUE) {
  mu_list <- lapply(mu, function(x) x$mu)
  D_list <- lapply(mu, function(x) x$D)
  output <- list()
  n_resp <- length(mu_list)
  n <- length(mu_list[[1]])
  idx = rep(unlist(variance), each = n)
  Omega <- mt_build_omega(tau = tau, Z = Z)
  variance_work <- variance
  power_fixed_w <- lapply(power_fixed, function(x)!x)
  variance_work[which(variance_work == "poisson_tweedie")] <- "tweedie"
  V.features <- Map(mt_variance_function, mu = mu_list, power = power,
                    Ntrial = Ntrial, #weights = weights,
                    variance = variance_work,
                    derivative_power = power_fixed_w,
                    derivative_mu = compute_derivative_beta)
  V_sqrt <- bdiag(lapply(V.features, function(x)x$V_sqrt))
  V_O <- Omega$Omega%*%V_sqrt
  Sigma <- V_sqrt%*%V_O
  mu_vec <- do.call(c, mu_list)
  Matrix::diag(Sigma) <- Matrix::diag(Sigma) + mu_vec*as.numeric(idx == "poisson_tweedie")
  D_Sigma_tau <- lapply(Omega$D_Omega,
                        function(x,V_sqrt)V_sqrt%*%x%*%V_sqrt,
                        V_sqrt = V_sqrt)
  D_V_sqrt_power <- lapply(V.features,function(x)x$D_V_sqrt_power)
  D_V_sqrt_power <- D_V_sqrt_power[unlist(power_fixed_w)]
  if(n_resp > 1) {
    ID <- as.factor(seq(1:n_resp))
    X <- model.matrix(~ ID -1)
    DD <- apply(X, 1, Diagonal, n = n_resp)
    DD <- lapply(DD, function(x)Matrix(x, sparse = TRUE))
    DD2 <- DD[unlist(power_fixed_w)]
  }
  if(n_resp == 1) {
    DD <- list(Matrix(1, ncol = 1, nrow = 1))
    DD2 <- Matrix(1, ncol = 1, nrow = 1)
    DD2 <- DD2[unlist(power_fixed_w)]
  }
  D_V_sqrt_power <- Map(kronecker, DD2, D_V_sqrt_power)
  D_Sigma_power <- lapply(D_V_sqrt_power,
                          function(x, V_O)
                            x%*%V_O + Matrix::t(V_O)%*%x,
                          V_O = V_O)
  if(length(D_Sigma_power) != 0) {
    output$D_Sigma <- D_Sigma_power
  }
  output$D_Sigma <- c(output$D_Sigma, D_Sigma_tau)
  output$Sigma <- Matrix::forceSymmetric(Sigma)
  if(inverse == TRUE) {
    output$inv_Sigma <- Matrix::chol2inv(Matrix::chol(output$Sigma))
  }
  if(compute_derivative_beta == TRUE) {
    D_V_sqrt_mu <- lapply(V.features, function(x)x$D_V_sqrt_mu)
    D_V_sqrt_beta <- Map(mt_derivative_V_sqrt_beta, D = D_list, D_V_sqrt_mu)
    output_temp <- list()
    for(i in 1:n_resp) {
      D_V_sqrt_beta_i <- lapply(D_V_sqrt_beta[[i]],
                                function(x,DDi)kronecker(DDi, x), DDi = DD[[i]])
      output_temp[[i]] <- lapply(D_V_sqrt_beta_i,
                                     function(x, V_O)
                                       x%*%V_O + Matrix::t(V_O)%*%x,
                                     V_O = V_O)
      if(variance[i] == "poisson_tweedie") {
        DD_elem <- apply(D_list[[i]], 2,
                         function(x, DDi)kronecker(DDi, x),
                         DDi = Matrix::diag(DD[[i]]))
        n_obs <- dim(V_sqrt)[1]
        n_beta <- dim(DD_elem)[2]
        for(j in 1:n_beta) {
          diag(output_temp[[i]][[j]]) <- DD_elem[,j] + Matrix::diag(output_temp[[i]][[j]])
        }
      }
    }
    output$D_Sigma_beta <- unlist(output_temp)
  }
  return(output)
}

