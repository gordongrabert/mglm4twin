## Loading extra packages
require(Matrix)
require(SimCorMultRes)
require(mglm4twin) ### devtools::load_all()
require(ggplot2)
require(ggExtra)
require(AGHmatrix)
require(MASS)
require(mvnfast)


#### Use SNP data from AGHmatrix #####

data(snp.pine)

# Step 1: Define Sample Size & Covariance Matrices
n <- 926 # Number of individuals

#Computing the additive relationship matrix based on VanRaden 2008
GRM <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9,
               maf=0.05, method="VanRaden")

# Search next PD matrix
GRM <- nearPD(GRM)$mat

pine.eigen <- eigen(GRM)
plot(pine.eigen$values)
title(main = "loblolly pine GRM")


# Genetic Covariance Matrix (Sigma_G)
Sigma_G <- matrix(c(0.6, 0.3, 0.2,  # Trait 1 variance & covariance
                    0.3, 0.5, 0.25, # Trait 2 variance & covariance
                    0.2, 0.25, 0.4),# Trait 3 variance & covariance
                  nrow = 3, byrow = TRUE)


# Environmental Covariance Matrix (Sigma_E)
Sigma_E <- matrix(c(0.5, 0.2, 0.3,  # Trait 1 variance & covariance
                    0.2, 0.6, 0.15, # Trait 2 variance & covariance
                    0.3, 0.15, 0.5),# Trait 3 variance & covariance
                  nrow = 3, byrow = TRUE)


# Step 2: Construct the Full Covariance Matrix
Sigma_total <- as.matrix(kronecker(GRM, Sigma_G) + kronecker(diag(n), Sigma_E))

# Simulate Gaussian data

# mvrnorm() from MASS package might be considerably slower
#pheno = mvrnorm(n = 1, mu = rep(0, nrow(Sigma_total)), Sigma = Sigma_total)

# you may need to adjust number of cores
pheno = rmvn(n = 1, mu = rep(0, nrow(Sigma_total)), sigma = Sigma_total, ncores = 9)

# Reshape the vector into a matrix with n columns
pheno_col <- matrix(pheno, ncol = 3, byrow = TRUE)

# Create data frame
data = as.data.frame(pheno_col)

# Prepare for mglm4twin

# Fixed effects

# When there are no fixed effects

linear_pred_1 <- V1 ~ 1
linear_pred_2 <- V2 ~ 1
linear_pred_3 <- V3 ~ 1

n <- 926 # Number of individuals

# Build matrix linear predictor
mat <- mt_grm(n = n, grm = GRM, n_resp = 3, model = "AE", data = data)

res <- mglm4twin(linear_pred = c(linear_pred_1, linear_pred_2, linear_pred_3),
                 matrix_pred = c(mat),
                 data = as.data.frame(data))

# Extract estimated variances and covariances

A11 <- res$Covariance[1]
A22 <- res$Covariance[2]
A33 <- res$Covariance[3]
A21 <- res$Covariance[4]
A31 <- res$Covariance[5]
A32 <- res$Covariance[6]

E11 <- res$Covariance[7]
E22 <- res$Covariance[8]
E33 <- res$Covariance[9]
E21 <- res$Covariance[10]
E31 <- res$Covariance[11]
E32 <- res$Covariance[12]


# Construct estimated covariance matrix A
cov_matrix_A <- matrix(c(
  A11, A21, A31,
  A21, A22, A32,
  A31, A32, A33
), nrow = 3, byrow = TRUE)

cov_matrix_A

# Construct estimated covariance matrix E
cov_matrix_E <- matrix(c(
  E11, E21, E31,
  E21, E22, E32,
  E31, E32, E33
), nrow = 3, byrow = TRUE)

cov_matrix_E


### For checking runtime for different sample sizes ##

##### Estimate runtime for phenotype c = 3 ####
# Define sequence of n values
n_values <- seq(100, 900, by = 200)  # Adjust upper limit if needed
times <- numeric(length(n_values))  # Empty vector to store times

# Loop over different values of n
for (i in seq_along(n_values)) {
  n <- n_values[i]  # Set current n

  start_time <- Sys.time()  # Start timing

  # Run the model for the current n
  res <- mglm4twin(
    linear_pred = c(linear_pred_1, linear_pred_2, linear_pred_3),
    matrix_pred = c(mt_grm(n = n, grm = GRM[1:n, 1:n], n_resp = 3, model = "AE", data = data[1:n,])),
    data = as.data.frame(data[1:n,])
  )

  end_time <- Sys.time()  # End timing

  times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))  # Store elapsed time
}

# Combine results into a data frame
time_results <- data.frame(n = n_values, Time_in_seconds = times)



# Optional: Plot computation time vs n
library(ggplot2)
ggplot(time_results, aes(x = n, y = Time_in_seconds)) +
  geom_line() + geom_point() +
  labs(title = "Computation Time vs n, Phenotypes = 3, Type = Gaussian",
       x = "Sample Size (n)",
       y = "Time (seconds)") +
  theme_minimal()



#### Speeding up computation by using Eigenvalue decomposition procedure by De Vlaming et al. (2022)  ####


##### EVD of A ####

mt_evd <- function(n, grm, n_resp, resp.m = NULL, model, n_pc = 20, formula = NULL, data = NULL){

  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################

  grm_sparse <- Matrix(grm, sparse = FALSE)
  A <- as(grm_sparse, "dsCMatrix")

  # Eigen decomposition of GRM
  evd <- eigen(A, symmetric = TRUE)
  Q <- evd$vectors
  lambda <- evd$values

  # Reduce Q
  n <- length(lambda) - n_pc
  P <- Q[, (n_pc + 1):length(lambda)]
  D <- diag(lambda[(n_pc + 1):length(lambda)])

  # Transform response variable matrix
  trans.resp.m <- t(P) %*% as.matrix(resp.m)
  E <- diag(nrow(trans.resp.m))


  ####################################################################
  ## Extending to multivariate responses #############################
  ####################################################################
  output <- list()

  if (n_resp > 1) {
    Z_struc <- mglm4twin:::mt_struc(n_resp = n_resp)
    ind_A <- lapply(Z_struc, function(x) kronecker(D, x))
    ind_E <- lapply(Z_struc, function(x) kronecker(E, x))
  }

  ####################################################################
  ## Selecting the different twin models #############################
  ####################################################################
  if (n_resp > 1) {
    if (model == "AE") {
      output$matrices <- c(ind_A, ind_E)
      output$phenotype <- trans.resp.m
    }
  } else {
    if (model == "E") {
      output <- list(ind_E)
    } else if (model == "AE") {
      output$matrices <- c(ind_A, ind_E)
      output$phenotype <- trans.resp.m
    }
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


mat <- mt_evd(grm = GRM, n_resp = 3, resp.m = as.data.frame(pheno_col), n_pc = 20, model = "AE", data = data)

# Prepare for mglm4twin

# Fixed effects

linear_pred_1 <- V1 ~ 1
linear_pred_2 <- V2 ~ 1
linear_pred_3 <- V3 ~ 1

data <- as.data.frame(mat$phenotype)

res <- mglm4twin(linear_pred = c(linear_pred_1, linear_pred_2, linear_pred_3),
                 matrix_pred = c(mat$matrices),
                 data = data)


A11 <- res$Covariance[1]
A22 <- res$Covariance[2]
A33 <- res$Covariance[3]
A21 <- res$Covariance[4]
A31 <- res$Covariance[5]
A32 <- res$Covariance[6]

E11 <- res$Covariance[7]
E22 <- res$Covariance[8]
E33 <- res$Covariance[9]
E21 <- res$Covariance[10]
E31 <- res$Covariance[11]
E32 <- res$Covariance[12]



# Construct covariance matrix
cov_matrix_A <- matrix(c(
  A11, A21, A31,
  A21, A22, A32,
  A31, A32, A33
), nrow = 3, byrow = TRUE)

cov_matrix_A

# Construct covariance matrix
cov_matrix_E <- matrix(c(
  E11, E21, E31,
  E21, E22, E32,
  E31, E32, E33
), nrow = 3, byrow = TRUE)

cov_matrix_E






