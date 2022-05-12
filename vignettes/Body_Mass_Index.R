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
data(bmi)
head(bmi)

## -----------------------------------------------------------------------------
linear_pred <- bmi ~ Group*Twin_pair

## -----------------------------------------------------------------------------
ACE = mt_twin(N_DZ = 5576/2, N_MZ = 2966/2, n_resp = 1, model = "ACE")
AE = mt_twin(N_DZ = 5576/2, N_MZ = 2966/2, n_resp = 1, model = "AE")

## -----------------------------------------------------------------------------
bmi$age <- (bmi$age - mean(bmi$age))/sd(bmi$age)
list_form <- list(formE = ~ age + gender, formA = ~ age + gender,
                  formC = ~ age + gender)

## -----------------------------------------------------------------------------
ACE_reg = mt_twin(N_DZ = 5576/2, N_MZ = 2966/2, n_resp = 1,
                  model = "ACE", formula = list_form, data = bmi)

## -----------------------------------------------------------------------------
list_form2 <- list(formE = ~ age + gender, formA = ~ age + gender)
AE_reg = mt_twin(N_DZ = 5576/2, N_MZ = 2966/2, n_resp = 1,
                 model = "AE", formula = list_form2, data = bmi)

## -----------------------------------------------------------------------------
link = "identity"
variance = "constant"

## -----------------------------------------------------------------------------
## Standard ACE model
fit_ACE <- mglm4twin(linear_pred = c(linear_pred), 
                     matrix_pred = ACE, data = bmi)

## Standard AE model
fit_AE <- mglm4twin(linear_pred = c(linear_pred), 
                    matrix_pred = AE, data = bmi)

## ACE regression on the dispersion
fit_ACE_reg <- mglm4twin(linear_pred = c(linear_pred), 
                         matrix_pred = ACE_reg, data = bmi)

## AE regression on the dispersion
fit_AE_reg <- mglm4twin(linear_pred = c(linear_pred), 
                        matrix_pred = AE_reg, data = bmi)

