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
data(psydis)
head(psydis)

## -----------------------------------------------------------------------------
linear_pred <- y ~ 1

## -----------------------------------------------------------------------------
ex1_E <- mt_twin(N_DZ = 440, N_MZ = 590, n_resp = 1, model = "E")
ex1_AE <- mt_twin(N_DZ = 440, N_MZ = 590, n_resp = 1, model = "AE")
ex1_CE <- mt_twin(N_DZ = 440, N_MZ = 590, n_resp = 1, model = "CE")
ex1_ACE <- mt_twin(N_DZ = 440, N_MZ = 590, n_resp = 1, model = "ACE")

## -----------------------------------------------------------------------------
link = "logit"
variance = "binomialP"

## -----------------------------------------------------------------------------
fit_ex1_E <- mglm4twin(linear_pred = c(linear_pred), matrix_pred = ex1_E, 
                       link = c(link), variance = c(variance), data = psydis)
fit_ex1_AE <- mglm4twin(linear_pred = c(linear_pred), matrix_pred = ex1_AE, 
                        link = c(link), variance = c(variance), data = psydis)
fit_ex1_CE <- mglm4twin(linear_pred = c(linear_pred), matrix_pred = ex1_CE, 
                        link = c(link), variance = c(variance), data = psydis)
fit_ex1_ACE <- mglm4twin(linear_pred = c(linear_pred), matrix_pred = ex1_ACE, 
                         link = c(link), variance = c(variance), data = psydis)

## -----------------------------------------------------------------------------
summary(fit_ex1_ACE, model = "ACE") ## Standard output

## -----------------------------------------------------------------------------
summary(fit_ex1_ACE, model = "ACE", biometric = TRUE) ## Standard output

## ---- echo = FALSE, results = "hide"------------------------------------------
Estimates <- cbind(c(coef(fit_ex1_E, model = "E")$Estimates[2], rep(NA, 2)),
c(coef(fit_ex1_AE, model = "AE")$Estimates[2:3], rep(NA, 1)),
c(coef(fit_ex1_CE, model = "CE")$Estimates[2], 0, coef(fit_ex1_CE, model = "CE")$Estimates[3]),
coef(fit_ex1_ACE, model = "ACE")$Estimates[2:4])

Std_error <- cbind(c(coef(fit_ex1_E, model = "E", 
                          std.error = TRUE)$std.error[2], rep(NA, 2)),
                   c(coef(fit_ex1_AE, model = "AE", 
                          std.error = TRUE)$std.error[2:3], rep(NA, 1)),
                   c(coef(fit_ex1_CE, model = "CE", 
                          std.error = TRUE)$std.error[2], 0, 
                     coef(fit_ex1_CE, model = "CE", 
                          std.error = TRUE)$std.error[3]),
                   coef(fit_ex1_ACE, model = "ACE", 
                        std.error = TRUE)$std.error[2:4])
ll <- c(as.numeric(plogLik(fit_ex1_E, verbose = FALSE)), 
        as.numeric(plogLik(fit_ex1_AE, verbose = FALSE)),
        as.numeric(plogLik(fit_ex1_CE, verbose = FALSE)),
        as.numeric(plogLik(fit_ex1_ACE, verbose = FALSE)))
Table1 <- cbind(Estimates[,1], Std_error[,1], Estimates[,2], Std_error[,2],
                Estimates[,3], Std_error[,3], Estimates[,4], Std_error[,4])
h2 <- summary(fit_ex1_AE, biometric = TRUE, model = "AE")
h2_ACE <- summary(fit_ex1_ACE, biometric = TRUE, model = "ACE")
H <- c(NA, NA, c(h2$A_main$Estimates, h2$A_main$std.error), NA, NA,
       c(h2_ACE$A_main$Estimates, h2_ACE$A_main$std.error))
Table1 <- round(rbind(Table1, H, ll), 2)
rownames(Table1) <- c("E", "A", "C", "h2", "ll")
colnames(Table1) <- c("Estimate", "Std_Error", "Estimate", "Std_Error", 
                      "Estimate", "Std_Error", "Estimate", "Std_Error")

## -----------------------------------------------------------------------------
Table1

