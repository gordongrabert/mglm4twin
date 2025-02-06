# Load libraries
library(MASS)
library(copula)
library(ggplot2)
library(ggExtra)
library(Matrix)
library(AGHmatrix)

#### Use SNP data from AGHmatrix #####

data(snp.pine)

# Step 1: Define Sample Size & Covariance Matrices
n <- 926 # Number of individuals

#Computing the additive relationship matrix based on VanRaden 2008
GRM <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9,
               maf=0.05, method="VanRaden")

GRM <- nearPD(GRM)$mat


# Genetic Covariance Matrix (Sigma_G)
Sigma_G <- matrix(c(0.7, 0.5,  # Trait 1 variance & covariance
                    0.5, 0.3), # Trait 2 variance & covariance
                  nrow = 2, byrow = TRUE)

# Environmental Covariance Matrix (Sigma_E)
Sigma_E <- matrix(c(0.3, 0.2,  # Trait 1 variance & covariance
                    0.2, 0.7), # Trait 2 variance & covariance
                  nrow = 2, byrow = TRUE)

# Step 2: Simulate from Gaussian Copula
# Construct full covariance matrix for copula

# Step 2: Construct the Full Covariance Matrix
Sigma_total <- kronecker(GRM, Sigma_G) + kronecker(diag(n), Sigma_E)

Sigma_total <- nearPD(Sigma_total)$mat

# Convert to a correlation matrix
Sigma_total <- cov2cor(Sigma_total)

# Step 3: Simulate from Gaussian Copula
copula_model <- normalCopula(param = P2p(Sigma_total), dim = nrow(Sigma_total), dispstr = "un")
U <- rCopula(1, copula_model)
# Reshape the vector into a matrix with 2 columns (will create pairs)
pheno_uni <- matrix(U, ncol = 2, byrow = TRUE)

# Step 4: Transform to Normal Marginals
Z <- qnorm(pheno_uni)  # Convert uniform [0,1] to standard normal
P <- qpois(pheno_uni, lambda = c(4,9))
PB <- data.frame(X1 = qpois(pheno_uni[,1], lambda = c(2,5)),
                 X2 = qbeta(pheno_uni[,2], shape1 = 2, shape2 = 3))
B <- qbeta(pheno_uni, shape1 = 5, shape2 = 5)

plot(Z, main = "Bivariate Gaussian", col = rgb(0, 0, 1, 0.5), pch = 16)
plot(P, main = "Bivariate Poisson", col = rgb(0, 0, 1, 0.5), pch = 16)
plot(B, main = "Bivariate Beta", col = rgb(0, 0, 1, 0.5), pch = 16)

# Create bivariate Gaussian plot with marginal distributions
p1 <- ggplot(data = data.frame(Z), aes(x = X1, y = X2)) +
  geom_point(alpha = 0.5, color = rgb(0, 0, 1, 0.5)) +
  ggtitle("Bivariate Gaussian") +
  xlab("Phenotype 1") + ylab("Phenotype 2") +
  theme_minimal(base_size = 30)
p1

# Add marginal distributions for Gaussian data
p1_marginal <- ggMarginal(p1, type = "density", fill = "blue")
print(p1_marginal)

# Create bivariate Poisson plot with marginal distributions
p2 <- ggplot(data = data.frame(P), aes(x = X1, y = X2)) +
  geom_point(alpha = 0.5, color = "red") +
  ggtitle("Bivariate Poisson") +
  theme_minimal()

# Add marginal distributions for Poisson data
p2_marginal <- ggMarginal(p2, type = "density", fill = "red")
print(p2_marginal)

# Create bivariate Beta plot with marginal distributions
p3 <- ggplot(data = data.frame(B), aes(x = X1, y = X2)) +
  geom_point(alpha = 0.5, color = rgb(0, 0, 1, 0.5)) +
  ggtitle("Bivariate Beta") +
  theme_minimal()

# Add marginal distributions for Poisson data
p3_marginal <- ggMarginal(p3, type = "density", fill = "orange")
print(p3_marginal)


# Create bivariate Beta plot with marginal distributions
p4 <- ggplot(data = data.frame(PB), aes(x = X1, y = X2)) +
  geom_point(alpha = 0.5, color = "red") +
  xlab("Poisson (Gene Expression)") + ylab("Beta (DNA Methylation)") +
  ggtitle("Bivariate Non-Gaussian") +
  theme_minimal(base_size = 30)

# Add marginal distributions for Poisson data
p4_marginal <- ggMarginal(p4, type = "density", fill = "red")
print(p4_marginal)


### McGLM ####

pheno = as.data.frame(Z)
V1 = pheno$V1
V2 = pheno$V2

I <- diag(n)
I <- as(I, "dtCMatrix")
linear_pred_1 <- V1 ~ 1
linear_pred_2 <- V2 ~ 1

#grm <- make_positive_definite_eigen(grm)
grm_sparse <- Matrix(GRM, sparse = F)
grm_dsC <- as(grm_sparse, "dsCMatrix")


A <- list(grm_dsC, I = I)

res <- mglm4twin(linear_pred = c(linear_pred_1),
                 matrix_pred = c(A),
                 data = as.data.frame(pheno))

A.est <- res$Covariance[1]
I.est <- res$Covariance[2]

print(c(A.est, I.est))


mt_twin <- function(N_DZ, N_MZ, n_resp, model, formula = NULL, data = NULL) {
  # Non-diagonal elements

  N_DZ = 4
  N_MZ = 4
  n_resp = 2

  a_dz <- 0.5
  d_dz <- 0.25
  ####################################################################
  ## Univariate matrices #############################################
  ####################################################################
  # MZ twin
  MZ_struc <- list()
  MZ_struc$A <- Matrix(c(1, 1, 1, 1), 2, 2)
  MZ_struc$C <- MZ_struc$A
  MZ_struc$D <- MZ_struc$A
  MZ_struc$E <- Diagonal(2, 1)
  # DZ twin
  DZ_struc <- list()
  DZ_struc$A <- Matrix(c(1, a_dz, a_dz, 1), 2, 2)
  DZ_struc$C <- MZ_struc$A
  DZ_struc$D <- Matrix(c(1, d_dz, d_dz, 1), 2, 2)
  DZ_struc$E <- Diagonal(2, 1)
  ####################################################################
  # Extending to the number of observed twins#########################
  ####################################################################
  I_DZ <- Diagonal(N_DZ, 1)
  I_MZ <- Diagonal(N_MZ, 1)
  DZ = lapply(DZ_struc, function(x, I_DZ)kronecker(I_DZ, x), I_DZ = I_DZ)
  MZ = lapply(MZ_struc, function(x, I_MZ)kronecker(I_MZ, x), I_MZ = I_MZ)
  Z_all <- Map(bdiag, DZ, MZ)
  ####################################################################
  ## Extending to multivariate responses #############################
  ####################################################################
  if(n_resp > 1) {
    Z_struc <- mt_struc(n_resp = n_resp)
    twin_A <- lapply(Z_struc, function(x, A)
      kronecker(x, A), A = Z_all$A)
    twin_C <- lapply(Z_struc, function(x, C)
      kronecker(x, C), C = Z_all$C)
    twin_D <- lapply(Z_struc, function(x, D)
      kronecker(x, D), D = Z_all$D)
    twin_E <- lapply(Z_struc, function(x, E)
      kronecker(x, E), E = Z_all$E)
  }
  ####################################################################
  ## Selecting the diferent twin models ##############################
  ####################################################################
  if(n_resp > 1) {
    if(model == "E") {
      output <- twin_E
    }
    if(model == "AE") {
      output <- c(twin_E, twin_A)
    }
    if(model == "CE") {
      output <- c(twin_E, twin_C)
    }
    if(model == "ACE") {
      output <- c(twin_E, twin_A, twin_C)
    }
    if(model == "ADE") {
      output <- c(twin_E, twin_A, twin_D)
    }
  }
  if(n_resp == 1) {
    if(model == "E") {
      output <- list(Z_all$E)
    }
    if(model == "AE") {
      output <- c(Z_all$E, Z_all$A)
    }
    if(model == "CE") {
      output <- c(Z_all$E, Z_all$C)
    }
    if(model == "ACE") {
      output <- c(Z_all$E, Z_all$A, Z_all$C)
    }
    if(model == "ADE") {
      output <- c(Z_all$E, Z_all$A, Z_all$D)
    }
  }
  if(!is.null(formula)) {
    if(length(output) != length(formula)) {
      print("Error: Number of formula does not match number of dispersion components")
    }
    if(length(output) == length(formula)) {
      X_list <- lapply(formula, model.matrix, data = data)
      new_output <- list()
      list_final <- list()
      for(i in 1:length(output)) {
        list_temp <- list()
        for(j in 1:ncol(X_list[[i]])) {
          list_temp[[j]] <- X_list[[i]][,j]*output[[i]]
        }
        list_final[[i]] <- list_temp
      }
      output <- do.call(c,list_final)
    }
  }
  return(output)
}


mt_GRM <- function(n, GRM, n_resp, model, formula = NULL, data = NULL) {
  ####################################################################
  ## Univariate matrices #############################################
  ####################################################################
  n = n
  grm_sparse <- Matrix(GRM, sparse = F)
  A <- as(grm_sparse, "dsCMatrix")
  E <- diag(n)
  model = "AE"
  formula = list("formE1" = ~ 1, "formE2" = ~ 1, "formE12" = ~ 1,
                 "formA1" = ~ 1, "formA2" = ~ 1, "formA12" = ~ 1)
  data = pheno

  ####################################################################
  # Extending to the number of observed twins#########################
  ####################################################################

  # I <- diag(n)
  #
  # ind = lapply(struc, function(x, I)kronecker(I, x), I = I)
  #
  # Z_all <- Map(bdiag, ind)

  ####################################################################
  ## Extending to multivariate responses #############################
  ####################################################################
  if(n_resp > 1) {
    Z_struc <-  mglm4twin:::mt_struc(n_resp = n_resp)
    ind_A <- lapply(Z_struc, function(x, A)
      kronecker(x, A), A = A)
    ind_E <- lapply(Z_struc, function(x, E)
      kronecker(x, E), E = E)
  }
  ####################################################################
  ## Selecting the diferent twin models ##############################
  ####################################################################
  if(n_resp > 1) {
    if(model == "AE") {
      output <- c(ind_E, ind_A)
    }

  }
  if(n_resp == 1) {
    if(model == "E") {
      output <- list(Z_all$E)
    }
    if(model == "AE") {
      output <- c(Z_all$E, Z_all$A)
    }
    if(model == "CE") {
      output <- c(Z_all$E, Z_all$C)
    }
    if(model == "ACE") {
      output <- c(Z_all$E, Z_all$A, Z_all$C)
    }
    if(model == "ADE") {
      output <- c(Z_all$E, Z_all$A, Z_all$D)
    }
  }
  if(!is.null(formula)) {
    if(length(output) != length(formula)) {
      print("Error: Number of formula does not match number of dispersion components")
    }
    if(length(output) == length(formula)) {
      X_list <- lapply(formula, model.matrix, data = data)
      new_output <- list()
      list_final <- list()
      for(i in 1:length(output)) {
        list_temp <- list()
        for(j in 1:ncol(X_list[[i]])) {
          list_temp[[j]] <- X_list[[i]][,j]*output[[i]]
        }
        list_final[[i]] <- list_temp
      }
      output <- do.call(c,list_final)
    }
  }
  return(output)
}


form <- list("formE1" = ~ 1, "formE2" = ~ 1, "formE12" = ~ 1,
             "formA1" = ~ 1, "formA2" = ~ 1, "formA12" = ~ 1)

mat <- mt_GRM(n = n, GRM = GRM, n_resp = 2, model = "AE")


res <- mglm4twin(linear_pred = c(linear_pred_1, linear_pred_2),
                 matrix_pred = c(mat),
                 data = as.data.frame(pheno))






link = rep("identity", 2)
variance = rep("constant", 2)

## Initial values
control_initial <- list()
control_initial$tau <- c(0.15, 0, 0.12, rep(0,3))

res <- mglm4twin(linear_pred = c(V1 ~ 1, V2 ~ 1),
                 matrix_pred = mat,
                 control_initial = control_initial,
                 power_fixed = c(TRUE, TRUE),
                 control_algorithm = list(tuning = 0.5),
                 data = as.data.frame(pheno[1:100,]))



### Poisson ###
pheno = as.data.frame(Z)
V1 = pheno$V1
V2 = pheno$V2

I <- diag(n)
I <- as(I, "dtCMatrix")
linear_pred_1 <- c(V1 ~ 1)
linear_pred_2 <- c(V2 ~ 1)

#grm <- make_positive_definite_eigen(grm)
grm_sparse <- Matrix(GRM, sparse = F)
grm_dsC <- as(grm_sparse, "dsCMatrix")

A <- list(grm_dsC, I = I)
res <- mglm4twin(linear_pred = c(linear_pred_1), matrix_pred = c(A), data = as.data.frame(pheno))

A.est <- res$Covariance[1]
I.est <- res$Covariance[2]

print(c(A.est, I.est))

### Beta ###
pheno = as.data.frame(B)
V1 = pheno$V1
V2 = pheno$V2

I <- diag(n)
I <- as(I, "dtCMatrix")
linear_pred_1 <- c(V1 ~ 1)
linear_pred_2 <- c(V2 ~ 1)

#grm <- make_positive_definite_eigen(grm)
grm_sparse <- Matrix(GRM, sparse = F)
grm_dsC <- as(grm_sparse, "dsCMatrix")


A1 <- list(grm_dsC, I = I)

link = rep("logit", 1)
variance = rep("binomialP", 1)

res <- mglm4twin(linear_pred = c(linear_pred_1), matrix_pred = c(A1), data = as.data.frame(pheno), link = link, variance = variance)

A.est <- res$Covariance[1]
I.est <- res$Covariance[2]

print(c(A.est, I.est))

data(anthro)
head(anthro)

## Standardize
anthro$age <- (anthro$age - mean(anthro$age))/sd(anthro$age)
anthro$weight <- (anthro$weight - mean(anthro$weight))/sd(anthro$weight)
anthro$height <- (anthro$height - mean(anthro$height))/sd(anthro$height)

form_Wt <- weight ~ age + Group*Twin_pair
form_Ht <- height ~ age + Group*Twin_pair

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


