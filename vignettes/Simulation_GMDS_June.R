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
n = 1000
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


mat <- mt_grm_1(n=n, grm = GRM, n_resp = 3, model = "AE", data = NULL)

mat <- mt_rsvd(n=n, grm = GRM, n_resp = 3, model = "AE", data = data_select)


form_Y1 <- c(Y1 ~ sex + age_std)
form_Y2 <- c(Y2 ~ sex + age_std)
form_Y3 <- c(Y3 ~ sex + age_std)

link = rep("logit", 3)
variance = rep("binomialP", 3)


res <- mglm4twin(linear_pred = c(form_Y1, form_Y2, form_Y3),
                 matrix_pred = c(mat),
                 link = link,
                 variance = variance,
                 data = data_select)

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
n_vals <- seq(200, 300, by = 100)

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
  
  # Generate qparameters
  phi <- 5
  qparameters <- vector("list", 3 * n)
  invcdfnames <- rep("qbeta", 3 * n)
  for (i in 1:n) {
    qparameters[[i]] <- list(shape1 = mu1[i]*phi, shape2 = (1 - mu1[i])*phi)
    qparameters[[n + i]] <- list(shape1 = mu2[i]*phi, shape2 = (1 - mu2[i])*phi)
    qparameters[[2 * n + i]] <- list(shape1 = mu3[i]*phi, shape2 = (1 - mu3[i])*phi)
  }
  
  # Simulate multivariate response
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
}

# View results
print(results)


# Optional: Plot computation time vs n
library(ggplot2)
ggplot(results, aes(x = n, y = runtime_sec)) +
  geom_line() + geom_point() +
  labs(title = "Computation Time vs n, Phenotypes = 3, Type = Multibound",
       x = "Sample Size (n)",
       y = "Time (seconds)") +
  theme_minimal()

ggplot(results, aes(x = n, y = frob_norm)) +
  geom_line() + geom_point() +
  labs(title = "Frobenius Norm (True h2 - Estimated h2) vs n, Phenotypes = 3, Type = Multibound",
       x = "Sample Size (n)",
       y = "Frobenius Norm") +
  theme_minimal()


