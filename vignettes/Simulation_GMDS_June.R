## Loading extra packages
require(Matrix)
require(SimCorMultRes)
require(mglm4twin) ### devtools::load_all()
require(ggplot2)
require(ggExtra)
require(AGHmatrix)
require(MASS)
require(mvnfast)
require(multivarious)
require(Rdimtools)
require(mildsvm)
require(rsvd)
require(rpca)
library(dplyr)
library(jsonlite)
library(RSpectra)
require(GGally)

#### Load Simulated Data ####

### read data from epigen ###
data <- read_json("/Users/gordonplri/Documents/epigen/epigen/sim/0_1_ASW.json")

### data extraction

df <- t(data.frame(matrix(unlist(data$genotype), nrow=length(data$genotype), byrow=TRUE)))


# Randomly select n rows
set.seed(123)
n = 5000
selected_rows <- sample(nrow(df), n)
df_selected <- df[selected_rows, ]  # Rows selected here

#Computing the additive relationship matrix based on VanRaden 2008
GRM <- Gmatrix(SNPmatrix=df_selected, missingValue=-9,
               maf=0.05, method="VanRaden")

# Search next PD matrix !IMPORTANT!
GRM <- nearPD(GRM)$mat

### Optional: Eigenvalue Decomposition 
eigen <- eigen(GRM)
plot(eigen$values)
title(main = "Eigenvalues GRM")


#### Phenotype simulation ####

#### Wagner Paper: Multi-bounded data

## Setting model parameters
E <- c(0.75, 0.7, 0.65, -0.3, 0.25, -0.4)
A <- c(0.25,0.3, 0.35, -0.15, 0.20, -0.2)
tau = c(E, A)

## Groundtruth 

# Heritability 
# h² values for lower triangle (row-wise order)
h2_vals <- A / (A + E)

# Fill into 3×3 matrix
h2_matrix <- matrix(0, 3, 3)
h2_matrix[lower.tri(h2_matrix, diag = TRUE)] <- h2_vals
h2_matrix <- h2_matrix + t(h2_matrix) - diag(diag(h2_matrix))

# Environmentatility  

e2_vals <- E / (A + E)
e2_matrix <- matrix(0, 3, 3)
e2_matrix[lower.tri(e2_matrix, diag = TRUE)] <- e2_vals
e2_matrix <- e2_matrix + t(e2_matrix) - diag(diag(e2_matrix))


## GRM structure

mt_grm_1 <- function(n, grm, n_resp, model, formula = NULL, data = NULL) {
  
  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################
  
  grm_sparse <- Matrix(grm, sparse = F)
  A <- as(grm_sparse, "dsCMatrix")
  E <- diag(n)
  # n = n
  # model = "AE"
  # formula = list("formE1" = ~ 1, "formE2" = ~ 1, "formE12" = ~ 1,
  #                "formA1" = ~ 1, "formA2" = ~ 1, "formA12" = ~ 1)
  # data = pheno
  
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
      output <- c(ind_E, ind_A) # TODO changed ordering of A and E matrix compared to mt_twin
    }
  }
  if(n_resp == 1) {
    if(model == "E") {
      output <- list(ind_E)
    }
    if(model == "AE") {
      output <- c(ind_E, ind_A)
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

GRM.cor <- cov2cor(GRM)
mat <- mt_grm_1(n = n, grm = GRM.cor, n_resp = 3, model = "AE", data = NULL)
Omega <- as.matrix(mt_matrix_linear_predictor(tau = tau, Z = mat))



# Create correlation matrix

set.seed(123)

## Regression structures
# Balanced but shuffled
sex <- sample(rep(c("Male", "Female"), each = n / 2))
trt <- sample(rep(c("Control", "Treatment"), each = n / 2))

Age <- rbeta(n, shape1 = 0.3*2, shape2 = 0.7*2)*20 + 70
age_std <- (Age - mean(Age))/var(Age)
X <- model.matrix(~ sex + age_std)
beta1 <- c(1.6956, 0.0584, -0.2576)
mu1 <- exp(X%*%beta1)/(1+exp(X%*%beta1))
beta2 <- c(-1.7930, 0.0875, 0.2382)
mu2 <- exp(X%*%beta2)/(1+exp(X%*%beta2))
beta3 <- c(-0.5363, -0.05138, -0.1528)
mu3 <- exp(X%*%beta3)/(1+exp(X%*%beta3))

## Marginal distributions
dist <- c("qbeta","qbeta","qbeta")
invcdfnames <- rep(dist, each = n)

## Simulating data set
Y1 <- list()
Y2 <- list()
Y3 <- list()

set.seed(181185)
# Pre-allocate qparameters list and names

qparameters <- vector("list", 3 * n)
names_vec <- character(3 * n)
phi <- 5

# Fill qparameters and names
for (i in 1:n) {
  # m1
  qparameters[[i]] <- list(shape1 = mu1[i] * phi, shape2 = (1 - mu1[i]) * phi)
  names_vec[i] <- paste0("m1_", i)

  # m2
  qparameters[[n + i]] <- list(shape1 = mu2[i] * phi, shape2 = (1 - mu2[i]) * phi)
  names_vec[n + i] <- paste0("m2_", i)

  # m3
  qparameters[[2 * n + i]] <- list(shape1 = mu3[i] * phi, shape2 = (1 - mu3[i]) * phi)
  names_vec[2 * n + i] <- paste0("m3_", i)

  names(qparameters) <- names_vec
}


Y <- rnorta(R = 1, cor.matrix = Omega,
            distr = invcdfnames, qparameters = qparameters)

Y1 <- Y[1:n]
Y2 <- Y[(n+1):(2*n)]
Y3 <- Y[(2*n+1):(3*n)]


data <- data.frame("Y1" = Y1, "Y2" =  Y2 , "Y3" = Y3,
                   "trt" = trt, "sex" = sex,
                   "age_std" = age_std)

hist(data$Y1)
hist(data$Y2)
hist(data$Y3)


ggpairs(data, columns = 1:3, aes(color = sex, alpha = 0.5),
        upper = list(continuous = "points")) + papaja::theme_apa()


# Generate ggpairs plot
p <- ggpairs(
  data,
  columns = 1:3,
  aes(color = sex, alpha = 0.6),
  upper = list(continuous = wrap("points", size = 1.5)),
  lower = list(continuous = wrap("points", size = 1.5)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
) +
  papaja::theme_apa(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

p

#ggsave("vignettes/figures/correlation_structure_plot.png", plot = p, width = 7, height = 7, dpi = 600)
### fitting model


data_select <- data %>% select(Y1, Y2, Y3, sex , age_std)



mt_rsvd <- function(n, grm, n_resp, model, formula = NULL, data = NULL){

  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################


  grm_sparse <- Matrix(grm, sparse = FALSE)
  A <- as(grm_sparse, "dsCMatrix")


  res <- rsvd(A)
  d <- res$d
  u <- res$u
  v <- res$v



  D <- diag(d)
  P <- res$u

  # Transform response variable matrix
  projected.data <- data
  E <- diag(nrow(projected.data))


  ####################################################################
  ## Extending to multivariate responses #############################
  ####################################################################
  output <- list()

  if (n_resp > 1) {
    Z_struc <- mglm4twin:::mt_struc(n_resp = n_resp)
    ind_A <- lapply(Z_struc, function(x) kronecker(x, D))
    ind_E <- lapply(Z_struc, function(x) kronecker(x, E))
  }

  ####################################################################
  ## Selecting the different twin models #############################
  ####################################################################
  if (n_resp == 1) {
    if (model == "E") {
      output <- list(
        matrices = c(ind_E),
        data = projected.data
      )
    } else if (model == "AE") {
      output <- list(
        matrices = c(ind_E, ind_A),
        data = projected.data
      )
    }
  } else if (n_resp > 1 && model == "AE") {
    output <- list(
      matrices = c(ind_E, ind_A),
      data = projected.data
    )
  }

  ####################################################################
  ## Applying formula-based transformations ##########################
  ####################################################################
  if (!is.null(formula)) {
    if (length(output) != length(formula)) {
      stop("Error: Number of formulas does not match number of dispersion components")
    }
    X_list <- lapply(formula, model.matrix, data = data)
    list_final <- lapply(seq_along(output), function(i) {
      lapply(seq_len(ncol(X_list[[i]])), function(j) {
        X_list[[i]][, j] * output[[i]]
      })
    })
    output <- do.call(c, list_final)
  }

  return(output)
}
mt_evd <- function(n, grm, n_resp, model, formula = NULL, data = NULL){
  
  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################
  data_numeric <- data[, sapply(data, is.numeric)]

  grm_sparse <- Matrix(grm, sparse = FALSE)
  A <- as(grm_sparse, "dsCMatrix")
  
  res <- eigen(A)
  
  A <- diag(res$values)
  
  Q <- as.matrix(res$vectors)
  
  is_prob_column <- function(x) all(x > 0 & x < 1)
  
  data_prob   <- data_numeric[, sapply(data_numeric, is_prob_column)]
  data_nonprob <- data_numeric[, !sapply(data_numeric, is_prob_column)]

  # Apply logit only where valid
  data_logit <- log(data_prob / (1 - data_prob))
  # Combine again
  data_transformed <- cbind(data_logit, data_nonprob)

  
  # Project
  P <- Q %*% as.matrix(data_transformed)
  P <- as.data.frame(P)
  colnames(P) <- colnames(data_numeric)
  P$sex <- data$sex

  # Transform response variable matrix
  
  E <- diag(nrow(A))
  
  
  ####################################################################
  ## Extending to multivariate responses #############################
  ####################################################################
  output <- list()
  
  if (n_resp > 1) {
    Z_struc <- mglm4twin:::mt_struc(n_resp = n_resp)
    ind_A <- lapply(Z_struc, function(x) kronecker(x, A))
    ind_E <- lapply(Z_struc, function(x) kronecker(x, E))
  }
  
  ####################################################################
  ## Selecting the different twin models #############################
  ####################################################################
  if (n_resp == 1) {
    if (model == "E") {
      output$matrices <- c(ind_E)
      output$data <- P
      
    } else if (model == "AE") {
      output$matrices <- c(ind_E, ind_A)
      output$data <- P
      
    }
  } else if (n_resp > 1 && model == "AE") {
    output$matrices <-  c(ind_E, ind_A)
    output$data <- P
    
  }
  
  ####################################################################
  ## Applying formula-based transformations ##########################
  ####################################################################
  if (!is.null(formula)) {
    if (length(output) != length(formula)) {
      stop("Error: Number of formulas does not match number of dispersion components")
    }
    X_list <- lapply(formula, model.matrix, data = data)
    list_final <- lapply(seq_along(output), function(i) {
      lapply(seq_len(ncol(X_list[[i]])), function(j) {
        X_list[[i]][, j] * output[[i]]
      })
    })
    output <- do.call(c, list_final)
  }
  
  return(output)
}


mat <- mt_grm_1(n=n, grm = GRM, n_resp = 3, model = "AE", data = NULL)

mat <- mt_rsvd(n=n, grm = GRM, n_resp = 3, model = "AE", data = data_select)
mat <- mt_evd(n=n, grm = GRM, n_resp = 3, model = "AE", data = data_select)




# Generate ggpairs plot
p <- ggpairs(
  mat$data,
  columns = 1:3,
  aes(color = sex, alpha = 0.6),
  upper = list(continuous = wrap("points", size = 1.5)),
  lower = list(continuous = wrap("points", size = 1.5)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
) +
  papaja::theme_apa(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

p


form_Y1 <- c(Y1 ~ sex + age_std)
form_Y2 <- c(Y2 ~ sex + age_std)
form_Y3 <- c(Y3 ~ sex + age_std)

link = rep("logit", 3)
variance = rep("binomialP", 3)


res <- mglm4twin(linear_pred = c(form_Y1, form_Y2, form_Y3),
                 matrix_pred = c(mat$matrices),
                 data = mat$data)

sum <-summary(res, model = "AE", biometric = T)
sum

# Create an empty 3×3 matrix
h2_estimate <- matrix(0, 3, 3)

# Fill diagonal
diag(h2_estimate) <- sum$A_main$Estimates

# Fill lower triangle
h2_estimate[lower.tri(h2_estimate)] <-sum$A_cross$Estimates

# Symmetrize the matrix
h2_estimate <- h2_estimate + t(h2_estimate) - diag(diag(h2_estimate))

# View result
round(h2_estimate, 3)

# Comparison Matrix Norms

# Frobenius norm of the difference
norm(h2_matrix - h2_estimate, type = "F")

### Simulation Pipeline
library(microbenchmark)
library(Matrix)

# Store results
results <- data.frame(
  n = numeric(),
  runtime_sec = numeric(),
  frob_norm = numeric()
)

# Set seeds for reproducibility
set.seed(123)

# Sample sizes to loop through
n_vals <- seq(500, 3000, by = 500)
n_vals_new <- seq(550, 800, by = 50)

summaries_list <- list()

load("vignettes/simulation_results.Rdata")  # or "results_df.Rdata" and "summaries_list.Rdata"

for (n in n_vals_new) {
  cat("Running for n =", n, "\n")

  ## 1. Subset data
  selected_rows <- sample(nrow(df), n)
  df_selected <- df[selected_rows, ]
  
  ## 2. Compute GRM and make PD
  GRM <- Gmatrix(SNPmatrix = df_selected, missingValue = -9,
                 maf = 0.05, method = "VanRaden")
  GRM <- nearPD(GRM)$mat
  GRM_cor <- cov2cor(GRM)
  
  ## 3. Generate Omega
  E <- c(0.75, 0.7, 0.65, -0.3, 0.25, -0.4)
  A <- c(0.25, 0.3, 0.35, -0.15, 0.20, -0.2)
  tau <- c(E, A)
  mat <- mt_grm_1(n = n, grm = GRM_cor, n_resp = 3, model = "AE", data = NULL)
  Omega <- as.matrix(mt_matrix_linear_predictor(tau = tau, Z = mat))
  
  ## 4. Simulate covariates and responses
  sex <- sample(rep(c("Male", "Female"), each = n / 2))
  trt <- sample(rep(c("Control", "Treatment"), each = n / 2))
  Age <- rbeta(n, shape1 = 0.3*2, shape2 = 0.7*2)*20 + 70
  age_std <- (Age - mean(Age))/var(Age)
  X <- model.matrix(~ sex + age_std)
  
  beta1 <- c(1.6956, 0.0584, -0.2576)
  beta2 <- c(-1.7930, 0.0875, 0.2382)
  beta3 <- c(-0.5363, -0.05138, -0.1528)
  
  mu1 <- exp(X %*% beta1) / (1 + exp(X %*% beta1))
  mu2 <- exp(X %*% beta2) / (1 + exp(X %*% beta2))
  mu3 <- exp(X %*% beta3) / (1 + exp(X %*% beta3))
  
  phi <- 5
  qparameters <- vector("list", 3 * n)
  invcdfnames <- rep("qbeta", 3 * n)
  for (i in 1:n) {
    qparameters[[i]] <- list(shape1 = mu1[i]*phi, shape2 = (1 - mu1[i])*phi)
    qparameters[[n + i]] <- list(shape1 = mu2[i]*phi, shape2 = (1 - mu2[i])*phi)
    qparameters[[2 * n + i]] <- list(shape1 = mu3[i]*phi, shape2 = (1 - mu3[i])*phi)
  }
  
  Y <- rnorta(R = 1, cor.matrix = Omega, distr = invcdfnames, qparameters = qparameters)
  Y1 <- Y[1:n]
  Y2 <- Y[(n+1):(2*n)]
  Y3 <- Y[(2*n+1):(3*n)]
  data <- data.frame(Y1 = Y1, Y2 = Y2, Y3 = Y3, sex = sex, age_std = age_std)
  
  ## 5. Fit model and time it
  mat <- mt_grm_1(n=n, grm = GRM, n_resp = 3, model = "AE", data = NULL)
  
  runtime <- system.time({
    res <- mglm4twin(
      linear_pred = list(Y1 ~ sex + age_std, Y2 ~ sex + age_std, Y3 ~ sex + age_std),
      matrix_pred = mat,
      link = rep("logit", 3),
      variance = rep("binomialP", 3),
      data = data
    )
  })["elapsed"]
  
  sum <- summary(res, model = "AE", biometric = TRUE)
  
  ## 6. Reconstruct estimated h2 matrix
  h2_estimate <- matrix(0, 3, 3)
  diag(h2_estimate) <- sum$A_main$Estimates
  h2_estimate[lower.tri(h2_estimate)] <- sum$A_cross$Estimates
  h2_estimate <- h2_estimate + t(h2_estimate) - diag(diag(h2_estimate))
  
  ## 7. Ground truth h2 matrix
  h2_vals <- A / (A + E)
  h2_matrix <- matrix(0, 3, 3)
  h2_matrix[lower.tri(h2_matrix, diag = TRUE)] <- h2_vals
  h2_matrix <- h2_matrix + t(h2_matrix) - diag(diag(h2_matrix))
  
  ## 8. Frobenius norm
  frob <- norm(h2_matrix - h2_estimate, type = "F")
  
  ## 9. Store results
  results <- rbind(results, data.frame(n = n, runtime_sec = runtime, frob_norm = frob))
  summaries_list[[paste0("n_", n)]] <- sum
}
# View results
print(results)

# Save results to .Rdata
#save(results, summaries_list, file = "vignettes/simulation_results.Rdata")

# Define a consistent, publication-ready theme
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

# Plot 1: Computation time
p1 <- ggplot(results, aes(x = n, y = runtime_sec)) +
  geom_line(color = "#1f77b4", linewidth = 1) +
  geom_point(shape = 21, fill = "#1f77b4", size = 2) +
  labs(
    title = "Computation Time vs Sample Size (n)",
    subtitle = "Phenotypes = 3, Type = Multibound",
    x = "Sample Size (n)",
    y = "Computation Time (seconds)"
  ) +
  theme_pub

ggsave("vignettes/figures/computation_time_vs_n.png", plot = p1, width = 6, height = 4, dpi = 600)

# Plot 2: Frobenius norm
p2 <- ggplot(results, aes(x = n, y = frob_norm)) +
  geom_line(color = "#d62728", linewidth = 1) +
  geom_point(shape = 21, fill = "#d62728", size = 2) +
  labs(
    title = "Estimation Error vs Sample Size (n)",
    subtitle = expression("Frobenius Norm of (" * hat(h)^2 * " - True " * h^2 * "), Phenotypes = 3, Type = Multibound"),
    x = "Sample Size (n)",
    y = "Frobenius Norm"
  ) +
  theme_pub

ggsave("vignettes/figures/frob_norm_vs_n.png", plot = p2, width = 6, height = 4, dpi = 600)





# Store results
results <- data.frame(
  n = numeric(),
  runtime_sec = numeric(),
  frob_norm = numeric()
)

# Set seeds for reproducibility
set.seed(123)

# Sample sizes to loop through
n_vals <- seq(500, 3000, by = 500)
n_vals_new <- seq(550, 800, by = 50)

summaries_list <- list()

load("vignettes/simulation_results.Rdata")  # or "results_df.Rdata" and "summaries_list.Rdata"

for (n in n_vals_new) {
  cat("Running for n =", n, "\n")
  
  ## 1. Subset data
  selected_rows <- sample(nrow(df), n)
  df_selected <- df[selected_rows, ]
  
  ## 2. Compute GRM and make PD
  GRM <- Gmatrix(SNPmatrix = df_selected, missingValue = -9,
                 maf = 0.05, method = "VanRaden")
  GRM <- nearPD(GRM)$mat
  GRM_cor <- cov2cor(GRM)
  
  ## 3. Generate Omega
  E <- c(0.75, 0.7, 0.65, -0.3, 0.25, -0.4)
  A <- c(0.25, 0.3, 0.35, -0.15, 0.20, -0.2)
  tau <- c(E, A)
  mat <- mt_grm_1(n = n, grm = GRM_cor, n_resp = 3, model = "AE", data = NULL)
  Omega <- as.matrix(mt_matrix_linear_predictor(tau = tau, Z = mat))
  
  ## 4. Simulate covariates and responses
  sex <- sample(rep(c("Male", "Female"), each = n / 2))
  trt <- sample(rep(c("Control", "Treatment"), each = n / 2))
  Age <- rbeta(n, shape1 = 0.3*2, shape2 = 0.7*2)*20 + 70
  age_std <- (Age - mean(Age))/var(Age)
  X <- model.matrix(~ sex + age_std)
  
  beta1 <- c(1.6956, 0.0584, -0.2576)
  beta2 <- c(-1.7930, 0.0875, 0.2382)
  beta3 <- c(-0.5363, -0.05138, -0.1528)
  
  mu1 <- exp(X %*% beta1) / (1 + exp(X %*% beta1))
  mu2 <- exp(X %*% beta2) / (1 + exp(X %*% beta2))
  mu3 <- exp(X %*% beta3) / (1 + exp(X %*% beta3))
  
  phi <- 5
  qparameters <- vector("list", 3 * n)
  invcdfnames <- rep("qbeta", 3 * n)
  for (i in 1:n) {
    qparameters[[i]] <- list(shape1 = mu1[i]*phi, shape2 = (1 - mu1[i])*phi)
    qparameters[[n + i]] <- list(shape1 = mu2[i]*phi, shape2 = (1 - mu2[i])*phi)
    qparameters[[2 * n + i]] <- list(shape1 = mu3[i]*phi, shape2 = (1 - mu3[i])*phi)
  }
  
  Y <- rnorta(R = 1, cor.matrix = Omega, distr = invcdfnames, qparameters = qparameters)
  Y1 <- Y[1:n]
  Y2 <- Y[(n+1):(2*n)]
  Y3 <- Y[(2*n+1):(3*n)]
  data <- data.frame(Y1 = Y1, Y2 = Y2, Y3 = Y3, sex = sex, age_std = age_std)
  
  ## 5. Fit model and time it
  mat <- mt_grm_1(n=n, grm = GRM, n_resp = 3, model = "AE", data = NULL)
  
  runtime <- system.time({
    res <- mglm4twin(
      linear_pred = list(Y1 ~ sex + age_std, Y2 ~ sex + age_std, Y3 ~ sex + age_std),
      matrix_pred = mat,
      link = rep("logit", 3),
      variance = rep("binomialP", 3),
      data = data
    )
  })["elapsed"]
  
  sum <- summary(res, model = "AE", biometric = TRUE)
  
  ## 6. Reconstruct estimated h2 matrix
  h2_estimate <- matrix(0, 3, 3)
  diag(h2_estimate) <- sum$A_main$Estimates
  h2_estimate[lower.tri(h2_estimate)] <- sum$A_cross$Estimates
  h2_estimate <- h2_estimate + t(h2_estimate) - diag(diag(h2_estimate))
  
  ## 7. Ground truth h2 matrix
  h2_vals <- A / (A + E)
  h2_matrix <- matrix(0, 3, 3)
  h2_matrix[lower.tri(h2_matrix, diag = TRUE)] <- h2_vals
  h2_matrix <- h2_matrix + t(h2_matrix) - diag(diag(h2_matrix))
  
  ## 8. Frobenius norm
  frob <- norm(h2_matrix - h2_estimate, type = "F")
  
  ## 9. Store results
  results <- rbind(results, data.frame(n = n, runtime_sec = runtime, frob_norm = frob))
  summaries_list[[paste0("n_", n)]] <- sum
}
# View results
print(results)


#### EVD model ####


# Store results
results.evd <- data.frame(
  n = numeric(),
  runtime_sec = numeric(),
  frob_norm = numeric()
)

# Set seeds for reproducibility
set.seed(123)

# Sample sizes to loop through
n_vals <- seq(3000, 6000, by = 1000)
#n_vals_new <- seq(550, 800, by = 50)

summaries_list.evd <- list()

#load("vignettes/simulation_results.Rdata")  # or "results_df.Rdata" and "summaries_list.Rdata"

for (n in n_vals) {
  cat("Running for n =", n, "\n")
  
  ## 1. Subset data
  selected_rows <- sample(nrow(df), n)
  df_selected <- df[selected_rows, ]
  
  ## 2. Compute GRM and make PD
  GRM <- Gmatrix(SNPmatrix = df_selected, missingValue = -9,
                 maf = 0.05, method = "VanRaden")
  GRM <- nearPD(GRM)$mat
  GRM_cor <- cov2cor(GRM)
  
  ## 3. Generate Omega
  E <- c(0.75, 0.7, 0.65, -0.3, 0.25, -0.4)
  A <- c(0.25, 0.3, 0.35, -0.15, 0.20, -0.2)
  tau <- c(E, A)
  mat <- mt_grm_1(n = n, grm = GRM_cor, n_resp = 3, model = "AE", data = NULL)
  Omega <- as.matrix(mt_matrix_linear_predictor(tau = tau, Z = mat))
  
  ## 4. Simulate covariates and responses
  sex <- sample(rep(c("Male", "Female"), each = n / 2))
  trt <- sample(rep(c("Control", "Treatment"), each = n / 2))
  Age <- rbeta(n, shape1 = 0.3*2, shape2 = 0.7*2)*20 + 70
  age_std <- (Age - mean(Age))/var(Age)
  X <- model.matrix(~ sex + age_std)
  
  beta1 <- c(1.6956, 0.0584, -0.2576)
  beta2 <- c(-1.7930, 0.0875, 0.2382)
  beta3 <- c(-0.5363, -0.05138, -0.1528)
  
  mu1 <- exp(X %*% beta1) / (1 + exp(X %*% beta1))
  mu2 <- exp(X %*% beta2) / (1 + exp(X %*% beta2))
  mu3 <- exp(X %*% beta3) / (1 + exp(X %*% beta3))
  
  phi <- 5
  qparameters <- vector("list", 3 * n)
  invcdfnames <- rep("qbeta", 3 * n)
  for (i in 1:n) {
    qparameters[[i]] <- list(shape1 = mu1[i]*phi, shape2 = (1 - mu1[i])*phi)
    qparameters[[n + i]] <- list(shape1 = mu2[i]*phi, shape2 = (1 - mu2[i])*phi)
    qparameters[[2 * n + i]] <- list(shape1 = mu3[i]*phi, shape2 = (1 - mu3[i])*phi)
  }
  
  Y <- rnorta(R = 1, cor.matrix = Omega, distr = invcdfnames, qparameters = qparameters)
  Y1 <- Y[1:n]
  Y2 <- Y[(n+1):(2*n)]
  Y3 <- Y[(2*n+1):(3*n)]
  data <- data.frame(Y1 = Y1, Y2 = Y2, Y3 = Y3, sex = sex, age_std = age_std)
  
  ## 5. Fit model and time it
  mat <- mt_evd(n=n, grm = GRM, n_resp = 3, model = "AE", data = NULL)
  
  runtime <- system.time({
    res <- mglm4twin(
      linear_pred = list(Y1 ~ sex + age_std, Y2 ~ sex + age_std, Y3 ~ sex + age_std),
      matrix_pred = mat,
      data = data
    )
  })["elapsed"]
  
  sum <- summary(res, model = "AE", biometric = TRUE)
  
  ## 6. Reconstruct estimated h2 matrix
  h2_estimate <- matrix(0, 3, 3)
  diag(h2_estimate) <- sum$A_main$Estimates
  h2_estimate[lower.tri(h2_estimate)] <- sum$A_cross$Estimates
  h2_estimate <- h2_estimate + t(h2_estimate) - diag(diag(h2_estimate))
  
  ## 7. Ground truth h2 matrix
  h2_vals <- A / (A + E)
  h2_matrix <- matrix(0, 3, 3)
  h2_matrix[lower.tri(h2_matrix, diag = TRUE)] <- h2_vals
  h2_matrix <- h2_matrix + t(h2_matrix) - diag(diag(h2_matrix))
  
  ## 8. Frobenius norm
  frob <- norm(h2_matrix - h2_estimate, type = "F")
  
  ## 9. Store results
  results.evd <- rbind(results.evd, data.frame(n = n, runtime_sec = runtime, frob_norm = frob))
  summaries_list.evd[[paste0("n_", n)]] <- sum
}
# View results
print(results.evd)










##### Factor model ####



# Matrix linear predictor ----------------------------------------------



# Simulate a random 3x3 covariance matrix
A <- matrix(rnorm(3), nrow = 3)
E <- matrix(diag(3), nrow = 3)
parameters <- (3 * 2)*2
tau <- rep(rnorm(parameters))
  
  
  
Z_struc <-  mglm4twin:::mt_struc(n_resp = 3)
ind_A <- lapply(Z_struc, function(x, A)
            kronecker(x, A), A = A)
ind_E <- lapply(Z_struc, function(x, E)
  kronecker(x, E), E = E)

Z = c(ind_E, ind_A)

mt_matrix_linear_predictor_factor_model <- function(tau, Z) {
  # Calculate number of traits p from tau vector
  p <- length(tau) / 2
  if (p != floor(p)) stop("Length of tau must be even (tau_c and tau_s for each trait).")
  
  # Split tau into tau_c and tau_s
  tau_c <- tau[1:p]
  tau_s <- tau[(p + 1):(2 * p)]
  
  # Number of unique elements in symmetric p x p matrix
  n_elements <- p * (p + 1) / 2
  if (length(Z) != n_elements) stop("Z must contain p(p+1)/2 elements.") # TODO check
  
  # Generate weighted matrices
  output <- list()
  idx <- 1
  for (i in 1:p) {
    # Diagonal (i == j)
    output[[idx]] <- (tau_c[i]^2 + tau_s[i]^2) * Z[[idx]]
    idx <- idx + 1
    # Off-diagonals (i < j)
    if (i < p) {
      for (j in (i + 1):p) {
        output[[idx]] <- (tau_c[i] * tau_c[j]) * Z[[idx]]
        idx <- idx + 1
      }
    }
  }
  
  return(Reduce("+", output))
}


# Run model
mt_matrix_linear_predictor_factor_model(tau, Z)
