#' @title Psychiatric disorders
#' @name psydis
#'
#' @description Psychiatric disorders in 1030 (440 DZ and 590 MZ) Caucasian
#'     female twin-pairs sampled from the Virginia Twin Registry.
#'     Lifetime psychiatric illness is a binary trait and was diagnosed using an
#'     adapted version of the Structured Clinical Interview for DSM-II-R Diagnosis.
#'
#' \itemize{
#'
#' \item \code{y} - Binary trait (disease presence YES - 1; NO - 0).
#'
#' \item \code{Group} - Twin zygosity (DZ - dizygotic; MZ - monozygotic).
#'
#' \item \code{Twin} - Code of twin pair.
#'
#' \item \code{Twin_pair} - Code of twin within the pair (1 and 2).
#'
#' }
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(psydis)
#'
#' @format a \code{data.frame} with 2060 records and 4 variables.
#'
#' @source Neale, M. C. and Maes, H. H. (2004) . Methodology for Genetic Studies
#'     of Twins and Families. Tech. rep., Virginia Common wealth University,
#'     Department of Psychiatry.
#'
#' @source Rabe-Hesketh, S., Skrondal, A. and Gjessing, H. K. (2008)
#'     Biometrical modeling of twin and family data using standard mixed
#'     model software. Biometrics, 64, 280–288.
#'
#' @source Bonat, W. H. and v. B. Hjelmborg, J. (2020) Multivariate Generalized
#'     Linear Models for Twin and Family data. to appear.
#'
#' @examples
#' require(mglm4twin)
#' data(psydis, package="mglm4twin")
#' ex1_form <- y ~ 1
#' ex1_AE <- mt_twin(N_DZ = 440, N_MZ = 590, n_resp = 1, model = "AE")
#' ex1_AE <- mglm4twin(c(ex1_form), matrix_pred = ex1_AE,
#'                     link = c("logit"), variance = c("binomialP"),
#'                    data = psydis)
#' summary(ex1_AE, model = "AE")
#' summary(ex1_AE, model = "AE", biometric = TRUE)
NULL

#' @title Body mass index
#' @name bmi
#'
#' @description It is a fairly common data set from the `mets` package.
#'     The dataset consists of 11188 observations, however, in the `mglm4twin`
#'     package we considered only paired twin-pairs. Thus, we opted to circulate
#'     the data in this new form to avoid mistakes. The resulting dataset consists
#'     of 4271(2788 DZ and 1483 MZ) twin-pairs.
#'
#' \itemize{
#'
#' \item \code{bmi} - Continuous trait (body mass index).
#'
#' \item \code{age} - Twin age.
#'
#' \item \code{gender} - Twin gender (male and female).
#'
#' \item \code{Group} - Twin zygosity (DZ - dizygotic; MZ - monozygotic).
#'
#' \item \code{Twin} - Code of twin pair.
#'
#' \item \code{Twin_pair} - Code of twin within the pair (1 and 2).
#'
#' }
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(bmi)
#'
#' @format a \code{data.frame} with 8542 records and 6 variables.
#'
#' @source Holst, K. K. and Scheike, T. H. and Hjelmborg, J. B. (2016).
#'     The Liability Threshold Model for Censored twin Data. Computational
#'     Statistics and Data Analysis 93, pp. 324-335. doi: 10.1016/j.csda.2015.01.014
#'
#' @source Bonat, W. H. and Hjelmborg, J. v. B. (2020) Multivariate Generalized
#'     Linear Models for Twin and Family data. to appear.
#'
#' @examples
#' require(mglm4twin)
#' data(bmi, package="mglm4twin")
#' form = bmi ~ Group*Twin_pair
#' ACE = mt_twin(N_DZ = 5576/2, N_MZ = 2966/2, n_resp = 1, model = "ACE")
#' fit_ACE <- mglm4twin(linear_pred = c(form), matrix_pred = ACE, data = bmi)
#  summary(fit_ACE, model = "ACE")
NULL

#' @title Bronchopulmonary dysplasia and respiratory distress syndrome on preterm infants
#' @name bpdrds
#'
#' @description We use the dataset analysed by Feng et al. (2009) regarding
#'     bronchopulmonary dysplasia (BPD) and respiratory distress syndrome (RDS)
#'     on preterm infants. Both diseases are lung related and expected to have a
#'     genetic component. The dataset consists of 200 twin-pairs being 137 DZ and 63 MZ.
#'     Additionally, we considered the covariates: birth weight (BW),
#'     gestation age (GA) and gender (female and male).
#'
#' \itemize{
#'
#' \item \code{Twin} - Code of twin pair.
#'
#' \item \code{gender} - Twin age gender (male and female).
#'
#' \item \code{GA} - Gestation age.
#'
#' \item \code{BW} - Birth weight.
#'
#' \item \code{RDS} - Respiratory distress syndrome (binary).
#'
#' \item \code{BPD} - Bronchopulmonary dysplasia (binary).
#'
#' \item \code{Group} - Twin zygosity (DZ - dizygotic; MZ - monozygotic).
#'
#' \item \code{Twin_pair} - Code of twin within the pair (1 and 2).
#'
#' }
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(bmi)
#'
#' @format a \code{data.frame} with 400 records and 8 variables.
#'
#' @source Feng, R., Zhou, G., Zhang, M. and Zhang, H. (2009)
#'     Analysis of twin data using sas. Biometrics, 65, 584–589.
#'
#' @source Bonat, W. H. and Hjelmborg, J. v. B. (2020) Multivariate Generalized
#'     Linear Models for Twin and Family data. to appear.
#'
#' @examples
#' require(mglm4twin)
#' data(bpdrds, package="mglm4twin")
#' form_BPD <- BPD ~ BW + GA + gender + Group*Twin_pair
#' form_RDS <- RDS ~ BW + GA + gender + Group*Twin_pair
#' AE <- mt_twin(N_DZ = 137, N_MZ = 63, n_resp = 2, model = "AE")
#' fitAE <- mglm4twin(linear_pred = c(form_BPD, form_RDS), matrix_pred = AE,
#'                    link = c("logit","logit"),
#'                    variance = c("binomialP","binomialP"), data = bpdrds)
NULL

#' @title Anthropometric measures (weight and height)
#' @name anthro
#'
#' @description Anthropometric measures (weight and height) on 861 (327 DZ and 534 MZ) twin-pairs.
#' Furthermore, we explore the flexibility of our proposed model class and model the dispersion.
#' The data set is available as an example in the OpenMx package (Neale et al., 2016).
#' We customize the data set for our needs, so make it available organized for
#' use in the mlm4twin package.
#'
#' \itemize{
#'
#' \item \code{weight} - Twin weight.
#'
#' \item \code{height} - Twin height.
#'
#' \item \code{age} - Twin age.
#'
#' \item \code{Group} - Twin zygosity (DZ - dizygotic; MZ - monozygotic).
#'
#'\item \code{Twin} - Twin code.
#'
#' \item \code{Twin_pair} - Code of twin within the pair (1 and 2).
#'
#' }
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(anthro)
#'
#' @format a \code{data.frame} with 1722 records and 6 variables.
#'
#' @source Neale, M. C., Hunter, M. D., Pritikin, J. N., Zahery, M., Brick, T. R.,
#'    Kirkpatrick, R. M., Estabrook, R., Bates, T. C., Maes, H. H. and Boker, S. M. (2016)
#'    OpenMx 2.0: Extended structural equation and statistical modeling.
#'    Psychometrika, 81, 535–549.
#'
#' @source Bonat, W. H. and Hjelmborg, J. v. B. (2020) Multivariate Generalized
#'     Linear Models for Twin and Family data. to appear.
#'
#' @examples
#' require(mglm4twin)
#' data(anthro, package="mglm4twin")
#' anthro$age <- (anthro$age - mean(anthro$age))/sd(anthro$age)
#' anthro$weight <- (anthro$weight - mean(anthro$weight))/sd(anthro$weight)
#' anthro$height <- (anthro$height - mean(anthro$height))/sd(anthro$height)
#' form_Wt <- weight ~ age + Group*Twin_pair
#' form_Ht <- height ~ age + Group*Twin_pair
#' biv0 <- list("formE1" = ~ age, "formE2" = ~ age, "formE12" = ~ age,
#'              "formA1" = ~ age, "formA2" = ~ age, "formA12" = ~ age,
#'              "formC1" = ~ age, "formC2" = ~ age, "formC12" = ~ age)
#' Z_biv0 <- mt_twin(N_DZ = 327, N_MZ = 534, n_resp = 2, model = "ACE",
#'                  formula = biv0, data = anthro)
#' control_initial <- list()
#' control_initial$regression <- list("R1" = c(0.13, 0.10, -0.20, -0.02, 0.037),
#'                                    "R2" = c(0.23, 0.01, -0.27, -0.11, 0.11))
#' control_initial$power <- list(c(0), c(0))
#' control_initial$tau <- c(0.15, 0, 0.12, rep(0,15))
#' fit_0 <- mglm4twin(linear_pred = c(form_Wt, form_Ht), matrix_pred = Z_biv0,
#'                    control_initial = control_initial,
#'                    control_algorithm = list(tuning = 0.5),
#'                    power_fixed = c(TRUE, TRUE), data = anthro)
NULL

#' @title Sleep’s quality
#' @name t0psqi
#'
#' @description Data set concerning sleep’s quality in a sample of
#'     250 (135 DZ and 116 MZ) Danish twin pairs. The traits are cortisone
#'     levels when waking up (T0) and PSQI (Pittsburgh Sleep Quality Index).
#'     It is a simulated data set based on the parameter estimates obtained
#'     fitting the model to a motivating real data set. The code for the simulation
#'     is available in the folder data-raw.
#'
#' \itemize{
#'
#' \item \code{Twin_pair} - Code of twin within the pair (1 and 2).
#'
#' \item \code{Twin_id} - Twin code.
#'
#' \item \code{Age} - Twin age.
#'
#' \item \code{Type} - Twin zygosity (DZ - dizygotic; MZ - monozygotic).
#'
#'\item \code{Gender} - Gender (Male and Female).
#'
#' \item \code{Group} - Treatment group, it is categorical covariate for composing
#'                      the linear predictor.
#'
#' \item \code{T0} - Cortisone levels when waking up (continuous trait).
#'
#' \item \code{PSQI} - Pittsburgh Sleep Quality Index (bounded trait) divided
#'                     by 21 (scale maximum).
#' }
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(t0psqi)
#'
#' @format a \code{data.frame} with 502 records and 8 variables.
#'
#' @source Bonat, W. H. and Hjelmborg, J. v. B. (2020) Multivariate Generalized
#'     Linear Models for Twin and Family data. to appear.
#'
#' @examples
#' require(mglm4twin)
#' form_T0 <- T0 ~ Age + Gender + Group + Type*Twin_pair
#' form_PSQI <- PSQI ~ Age + Gender + Group + Type*Twin_pair
#' AE <- mt_twin(N_DZ = 135, N_MZ = 116, n_resp = 2, model = "AE")
#' fit_AE <- mglm4twin(linear_pred = c(form_T0,  form_PSQI),
#'                    matrix_pred = AE,
#'                    link = c("log","logit"),
#'                    variance = c("tweedie","binomialP"),
#'                    control_algorithm = list(tuning = 0.25, max_iter = 100),
#'                    power_fixed = c(FALSE,FALSE), data = t0psqi)
NULL
