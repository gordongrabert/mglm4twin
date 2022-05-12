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
data(anthro)
head(anthro)

## Standardize
anthro$age <- (anthro$age - mean(anthro$age))/sd(anthro$age)
anthro$weight <- (anthro$weight - mean(anthro$weight))/sd(anthro$weight)
anthro$height <- (anthro$height - mean(anthro$height))/sd(anthro$height)


## -----------------------------------------------------------------------------
form_Wt <- weight ~ age + Group*Twin_pair
form_Ht <- height ~ age + Group*Twin_pair

## -----------------------------------------------------------------------------
## ACE model
biv0 <- list("formE1" = ~ age, "formE2" = ~ age, "formE12" = ~ age,
             "formA1" = ~ age, "formA2" = ~ age, "formA12" = ~ age,
             "formC1" = ~ age, "formC2" = ~ age, "formC12" = ~ age)
Z_biv0 <- mt_twin(N_DZ = 327, N_MZ = 534, n_resp = 2, model = "ACE",
                    formula = biv0, data = anthro)

## Special case of ACE model
biv4 <- list("formE1" = ~ age, "formE2" = ~ age, "formE12" = ~ 1,
             "formA1" = ~ age, "formA2" = ~ 1, "formA12" = ~ 1)
Z_biv4 <- mt_twin(N_DZ = 327, N_MZ = 534, n_resp = 2, model = "AE",
                  formula = biv4, data = anthro)

## -----------------------------------------------------------------------------
## Initial values
control_initial <- list()
control_initial$regression <- list("R1" = c(0.13, 0.10, -0.20, -0.02, 0.037),
                                   "R2" = c(0.23, 0.01, -0.27, -0.11, 0.11))
control_initial$power <- list(c(0), c(0))
control_initial$tau <- c(0.15, 0, 0.12, rep(0,15))

fit_0 <- mglm4twin(linear_pred = c(form_Wt, form_Ht), matrix_pred = Z_biv0,
                          control_initial = control_initial,
                          control_algorithm = list(tuning = 0.5),
                          power_fixed = c(TRUE, TRUE),
                          data = anthro)

control_initial$tau <- c(0.15, 0, 0.12, rep(0,6))
fit_4 <- mglm4twin(linear_pred = c(form_Wt, form_Ht), matrix_pred = Z_biv4,
                   control_initial = control_initial,
                   control_algorithm = list(tuning = 0.5),
                   power_fixed = c(TRUE, TRUE),
                   data = anthro)

