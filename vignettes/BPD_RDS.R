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
data(bpdrds)
head(bpdrds)

## -----------------------------------------------------------------------------
form_BPD <- BPD ~ BW + GA + gender + Group*Twin_pair
form_RDS <- RDS ~ BW + GA + gender + Group*Twin_pair

## -----------------------------------------------------------------------------
## Univariate models
uni_E <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 1, model = "E")
uni_AE <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 1, model = "AE")
uni_CE <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 1, model = "CE")
uni_ACE <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 1, model = "ACE")

## Bivariate models
biv_E <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 2, model = "E")
biv_AE <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 2, model = "AE")
biv_CE <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 2, model = "CE")
biv_ACE <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 2, model = "ACE")

## -----------------------------------------------------------------------------
link = c("logit", "logit")
variance = c("binomialP", "binomialP")

## -----------------------------------------------------------------------------
## Univariate fit

# Univariate E model
fitE_BPD <- mglm4twin(linear_pred = c(form_BPD), matrix_pred = uni_E,
                      link = c("logit"), variance = c("binomialP"),
                      data = bpdrds)
fitE_RDS <- mglm4twin(linear_pred = c(form_RDS), matrix_pred = uni_E,
                      link = c("logit"), variance = c("binomialP"),
                      data = bpdrds)

# Univariate AE model
fitAE_BPD <- mglm4twin(linear_pred = c(form_BPD), matrix_pred = uni_AE,
                      link = c("logit"), variance = c("binomialP"),
                      data = bpdrds)
fitAE_RDS <- mglm4twin(linear_pred = c(form_RDS), matrix_pred = uni_AE,
                       link = c("logit"), variance = c("binomialP"),
                       data = bpdrds)

# Univariate CE model
fitCE_BPD <- mglm4twin(linear_pred = c(form_BPD), matrix_pred = uni_CE,
                       link = c("logit"), variance = c("binomialP"),
                       data = bpdrds)
fitCE_RDS <- mglm4twin(linear_pred = c(form_RDS), matrix_pred = uni_CE,
                       link = c("logit"), variance = c("binomialP"),
                       data = bpdrds)

# Univariate ACE model
fitACE_BPD <- mglm4twin(linear_pred = c(form_BPD), matrix_pred = uni_ACE,
                       link = c("logit"), variance = c("binomialP"),
                       data = bpdrds)
fitACE_RDS <- mglm4twin(linear_pred = c(form_RDS), matrix_pred = uni_ACE,
                        link = c("logit"), variance = c("binomialP"),
                        data = bpdrds)

## Bivariate fit

# E model
fitE_biv <- mglm4twin(linear_pred = c(form_BPD, form_RDS), matrix_pred = biv_E, 
                      link = c("logit","logit"), 
                      variance = c("binomialP","binomialP"), 
                      data = bpdrds)
# AE model
fitAE_biv <- mglm4twin(linear_pred = c(form_BPD, form_RDS), matrix_pred = biv_AE, 
                       link = c("logit","logit"), 
                       variance = c("binomialP","binomialP"), 
                       data = bpdrds)
# CE model
fitCE_biv <- mglm4twin(linear_pred = c(form_BPD, form_RDS), matrix_pred = biv_CE, 
                       link = c("logit","logit"), 
                       variance = c("binomialP","binomialP"), 
                       data = bpdrds)
# ACE model
fitACE_biv <- mglm4twin(linear_pred = c(form_BPD, form_RDS), matrix_pred = biv_ACE, 
                        link = c("logit", "logit"),
                        variance = c("binomialP","binomialP"),
                        data = bpdrds)

## -----------------------------------------------------------------------------
## Univariate models
uni <- round(rbind("E" = gof(list(fitE_BPD, fitE_RDS)), 
                   "AE" = gof(list(fitAE_BPD, fitAE_RDS)) ,
                   "CE" = gof(list(fitCE_BPD, fitCE_RDS)),
                   "ACE" = gof(list(fitACE_BPD, fitACE_RDS))), 2)

## Bivariate models
multi <- round(rbind("E" = gof(fitE_biv), "AE" = gof(fitAE_biv), 
                     "CE" = gof(fitCE_biv), "ACE" = gof(fitACE_biv)), 2)


## ---- echo = FALSE, results = "hide"------------------------------------------
Table3 <- cbind(t(uni)[c(1,3,5,2),c(4,3,2,1)], t(multi)[c(1,3,5,2),c(4,3,2,1)])

## -----------------------------------------------------------------------------
Table3

