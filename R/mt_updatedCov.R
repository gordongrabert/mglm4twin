#' @title Updated covariance parameters
#' @author Wagner Hugo Bonat
#'
#' @description This function update a list of covariance parameters. It
#'     will be useful only inside the fitting algorithm.
#'
#' @param list_initial A list of initial values.
#' @param covariance A vector with actual covariance parameters values.
#' @param information A list with information about the number of
#'     parameters in the model. In general the output from
#'     \link{mt_getInformation}.
#' @param list_power_fixed A list of logicals indicating if the power
#'     parameter should be estimated or not.
#' @param n_resp A numeric specyfing the number of response variables.
#' @keywords internal
#' @return A list with updated values of the covariance parameters.

mt_updateCov <- function(list_initial, covariance, list_power_fixed,
                         information, n_resp) {
    power_cod <- list()
    tau_cod <- rep("tau", information$n_taus)
    for (i in 1:n_resp) {
        power_cod[[i]] <- rep(paste("power", i, sep = ""),
                              information$n_power[[i]])
    }
    dd <- data.frame("Parameter" = covariance,
                     "Cod" = c(unlist(power_cod), tau_cod))
    for(i in 1:n_resp) {
      if(list_power_fixed[i] == FALSE)
      list_initial$power[[i]] <- dd[which(dd$Cod == paste("power", i, sep = "")),]$Parameter
    }
    list_initial$tau <- dd[which(dd$Cod == "tau"),]$Parameter
    return(list_initial)
}
