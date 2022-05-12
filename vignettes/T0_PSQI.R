## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  install_git("wbonat/mglm4twin")

## ---- eval=TRUE, error=FALSE, message=FALSE, warning=FALSE--------------------
library(mglm4twin)
packageVersion("mglm4twin")

## -----------------------------------------------------------------------------
data(t0psqi)
head(t0psqi)

## -----------------------------------------------------------------------------
form_T0 <- T0 ~ Age + Gender + Group + Type*Twin_pair
form_PSQI <- PSQI ~ Age + Gender + Group + Type*Twin_pair

## -----------------------------------------------------------------------------
E <- mt_twin(N_DZ = 135, N_MZ = 116, n_resp = 2, model = "E")
AE <- mt_twin(N_DZ = 135, N_MZ = 116, n_resp = 2, model = "AE")
CE <- mt_twin(N_DZ = 135, N_MZ = 116, n_resp = 2, model = "CE")
ACE <- mt_twin(N_DZ = 135, N_MZ = 116, n_resp = 2, model = "ACE")

## -----------------------------------------------------------------------------
## Gaussian model (link = 'identity' and variance = 'constant')
Gauss_E <- mglm4twin(linear_pred = c(form_T0,  form_PSQI), 
                     matrix_pred = E, data = t0psqi)
Gauss_AE <- mglm4twin(linear_pred = c(form_T0,  form_PSQI), 
                      matrix_pred = AE, data = t0psqi)
Gauss_CE <- mglm4twin(linear_pred = c(form_T0,  form_PSQI), 
                      matrix_pred = CE, data = t0psqi)
Gauss_ACE <- mglm4twin(linear_pred = c(form_T0,  form_PSQI), 
                       matrix_pred = ACE, data = t0psqi)

## Mixed types
fit_E <- mglm4twin(linear_pred = c(form_T0,  form_PSQI), 
                   matrix_pred = E, 
                   link = c("log","logit"), 
                   variance = c("tweedie","binomialP"), 
                   control_algorithm = list(tuning = 0.25, max_iter = 100),
                   power_fixed = c(F,F), data = t0psqi)
fit_AE <- mglm4twin(linear_pred = c(form_T0,  form_PSQI), 
                    matrix_pred = AE, 
                    link = c("log","logit"), 
                    variance = c("tweedie","binomialP"), 
                    control_algorithm = list(tuning = 0.25, max_iter = 100),
                    power_fixed = c(F,F), data = t0psqi)
fit_CE <- mglm4twin(linear_pred = c(form_T0,  form_PSQI), 
                    matrix_pred = CE, 
                    link = c("log","logit"), 
                    variance = c("tweedie","binomialP"), 
                    control_algorithm = list(tuning = 0.25, max_iter = 100),
                    power_fixed = c(F,F), data = t0psqi)
fit_ACE <- mglm4twin(linear_pred = c(form_T0,  form_PSQI), 
                     matrix_pred = ACE, 
                     link = c("log","logit"), 
                     variance = c("tweedie","binomialP"), 
                     control_algorithm = list(tuning = 0.25, max_iter = 100),
                     power_fixed = c(F,F), data = t0psqi)

## ---- echo = FALSE, results = "hide"------------------------------------------
Table6 <- cbind(rbind(gof(Gauss_E), gof(Gauss_AE), gof(Gauss_CE), 
                      gof(Gauss_ACE))[,-c(4)],
                rbind(gof(fit_E), gof(fit_AE), gof(fit_CE), 
                      gof(fit_ACE))[,-c(4)])


## -----------------------------------------------------------------------------
Table6

## -----------------------------------------------------------------------------
summary(fit_AE, model = "AE", biometric = TRUE)

