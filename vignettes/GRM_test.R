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

# GRM <- nearPD(GRM)$mat


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

#Sigma_total <- nearPD(Sigma_total)$mat

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

#### McGLM ####

pheno = as.data.frame(Z)
V1 = pheno$V1
V2 = pheno$V2

I <- diag(n)
I <- as(I, "dtCMatrix")
linear_pred_1 <- V1 ~ 1
linear_pred_2 <- V2 ~ 1

grm_sparse <- Matrix(GRM, sparse = F)
grm_dsC <- as(grm_sparse, "dsCMatrix")


A <- list(grm_dsC, I = I)

res <- mglm4twin(linear_pred = c(linear_pred_1),
                 matrix_pred = c(A),
                 data = as.data.frame(pheno))

A.est <- res$Covariance[1]
I.est <- res$Covariance[2]


print(c(A.est, I.est))

### Multivariate GRM ###

mt_GRM <- function(n, GRM, n_resp, model, formula = NULL, data = NULL) {
  ####################################################################
  ## Univariate matrices #############################################
  ####################################################################
  # n = n
   grm_sparse <- Matrix(GRM, sparse = F)
   A <- as(grm_sparse, "dsCMatrix")
   E <- diag(n)
  # model = "AE"
  # formula = list("formE1" = ~ 1, "formE2" = ~ 1, "formE12" = ~ 1,
  #                "formA1" = ~ 1, "formA2" = ~ 1, "formA12" = ~ 1)
  # data = pheno

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
      output <- c(ind_A, ind_E)
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


n <- 926 # Number of individuals

mat <- mt_GRM(n = 500, GRM = GRM[1:500, 1:500], n_resp = 2, model = "AE", data = pheno[1:500,])


res <- mglm4twin(linear_pred = c(linear_pred_1, linear_pred_2),
                 matrix_pred = c(mat),
                 data = as.data.frame(pheno[1:500,]))


A1 <- res$Covariance[1]
A2 <- res$Covariance[2]
A3 <- res$Covariance[3]
E1 <- res$Covariance[4]
E2 <- res$Covariance[5]
E3 <- res$Covariance[6]


cov_table <- tibble::tibble(
  Component = c("A1", "A2", "A3", "E1", "E2", "E3"),
  Covariance = c(A1, A2, A3, E1, E2, E3)
)

cov_table

# Construct covariance matrix
cov_matrix_A <- matrix(c(
  A1, A3,
  A3, A2
), nrow = 2, byrow = TRUE)

cov_matrix_A

# Construct covariance matrix
cov_matrix_E <- matrix(c(
  E1, E3,
  E3, E2
), nrow = 2, byrow = TRUE)

cov_matrix_E

##### Estimate runtime for phenotype c = 1 #####


# Define sequence of n values
n_values <- seq(100, 900, by = 200)  # Adjust upper limit if needed
times_1 <- numeric(length(n_values))  # Empty vector to store times

# Loop over different values of n
for (i in seq_along(n_values)) {
  n <- n_values[i]  # Set current n

  start_time <- Sys.time()  # Start timing

  # Run the model for the current n
  res <- mglm4twin(
    linear_pred = c(linear_pred_1),
    matrix_pred = c(grm_dsC[1:n, 1:n], I[1:n, 1:n]),
    data = as.data.frame(pheno[1:n,])
  )

  end_time <- Sys.time()  # End timing

  times_1[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))  # Store elapsed time
}

# Combine results into a data frame
time_results_1 <- data.frame(n = n_values, Time_in_seconds = times_1)

# Print results
print(time_results_1)

# Optional: Plot computation time vs n
library(ggplot2)
ggplot(time_results_1, aes(x = n, y = Time_in_seconds)) +
  geom_line() + geom_point() +
  labs(title = "Computation Time vs n, Phenotypes = 1, Type = Gaussian",
       x = "Sample Size (n)",
       y = "Time (seconds)") +
  theme_minimal()



##### Estimate runtime for phenotype c = 2 ####
# Define sequence of n values
n_values <- seq(100, 900, by = 200)  # Adjust upper limit if needed
times <- numeric(length(n_values))  # Empty vector to store times

# Loop over different values of n
for (i in seq_along(n_values)) {
  n <- n_values[i]  # Set current n

  start_time <- Sys.time()  # Start timing

  # Run the model for the current n
  res <- mglm4twin(
    linear_pred = c(linear_pred_1, linear_pred_2),
    matrix_pred = c(mt_GRM(n = n, GRM = GRM[1:n, 1:n], n_resp = 2, model = "AE", data = pheno[1:n,])),
    data = as.data.frame(pheno[1:n,])
  )

  end_time <- Sys.time()  # End timing

  times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))  # Store elapsed time
}

# Combine results into a data frame
time_results <- data.frame(n = n_values, Time_in_seconds = times)
regress <- lm(Time_in_seconds ~ I(n^2) , data = time_results)
summary(regress)

time_calc <- function(n){-6.560e+01 + 1.027e-03*n^2}
time_calc(1000)

# Print results
print(time_results)

# Optional: Plot computation time vs n
library(ggplot2)
ggplot(time_results, aes(x = n, y = Time_in_seconds)) +
  geom_line() + geom_point() +
  labs(title = "Computation Time vs n, Phenotypes = 2, Type = Gaussian",
       x = "Sample Size (n)",
       y = "Time (seconds)") +
  theme_minimal()



res <- mglm4twin(linear_pred = c(linear_pred_1),
                 matrix_pred = c(A),
                 data = as.data.frame(pheno))
