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



mt_evd_2 <- function(n, grm, n_resp, resp.m = NULL, model, n_pc = 20, formula = NULL, data = NULL){

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
  P <- Q[, 1:n_pc]
  D <- diag(lambda[1:n_pc])

  # Transform response variable matrix
  trans.resp.m <- t(P) %*% as.matrix(resp.m)
  E <- diag(nrow(trans.resp.m))


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


mat <- mt_evd_2(grm = GRM, n_resp = 3, resp.m = as.data.frame(pheno_col), n_pc = 926, model = "AE", data = data)



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


###### GRM Tensor model ####

require(TensorTools)


#### Use SNP data from AGHmatrix #####

data(snp.pine)

# Step 1: Define Sample Size & Covariance Matrices
n <- 926 # Number of individuals

#Computing the additive relationship matrix based on VanRaden 2008
GRMa <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9,
               maf=0.05, method="VanRaden")

# Search next PD matrix
GRMa <- nearPD(GRMa)$mat


#Computing the dominance relationship matrix based on VanRaden 2008
GRMd <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9,
                maf=0.05, method="Su")

# Search next PD matrix
GRMd <- nearPD(GRMd)$mat

# Additive-by-Additive Interactions

GRMaa <- GRMa * GRMa

# Constructing the tensor with additive, dominance, and interaction components
GRMtensor <- array(c(as.matrix(GRMa), as.matrix(GRMd), as.matrix(GRMaa)),
                   dim = c(nrow(GRMa), ncol(GRMa), 3))

dim(GRMtensor)

# Convert to tensor format
GRMtensor <- as.Tensor(GRMtensor)
dim(GRMtensor)

# Eigen decomposition of the tensor
teigen <- tEIG(GRMtensor, "dst")

A_eigen <- teigen$D$data[, , 1]
D_eigen <- teigen$D$data[, , 2]
AA_eigen <- teigen$D$data[, , 3]


pine.eigen <- eigen(GRMaa)
plot(pine.eigen$values)
title(main = "loblolly pine GRM")


# Genetic Covariance Matrix (Sigma_A)
Sigma_A <- matrix(c(0.6, 0.3, 0.2,
                    0.3, 0.5, 0.25,
                    0.2, 0.25, 0.4),
                  nrow = 3, byrow = TRUE)

# Genetic Covariance Matrix (Sigma_D)
Sigma_D <- matrix(c(0.3, 0.2, 0.1,
                    0.2, 0.4, 0.15,
                    0.1, 0.15, 0.3),
                  nrow = 3, byrow = TRUE)

# Genetic Covariance Matrix (Sigma_AA)
Sigma_AA <- matrix(c(0.4, 0.1, 0.2,
                     0.1, 0.35, 0.1,
                     0.2, 0.1, 0.5),
                   nrow = 3, byrow = TRUE)

# Environmental Covariance Matrix (Sigma_E)
Sigma_E <- matrix(c(0.5, 0.2, 0.3,
                    0.2, 0.6, 0.15,
                    0.3, 0.15, 0.5),
                  nrow = 3, byrow = TRUE)


# Step 2: Construct the Full Covariance Matrix
Sigma_total <- as.matrix(kronecker(GRMa, Sigma_A) +
                         kronecker(GRMd, Sigma_D) +
                        # kronecker(GRMaa, Sigma_AA) +
                         kronecker(diag(n), Sigma_E))



#Simulate Gaussian data

# mvrnorm() from MASS package might be considerably slower
#pheno = mvrnorm(n = 1, mu = rep(0, nrow(Sigma_total)), Sigma = Sigma_total)

# you may need to adjust number of cores
pheno = rmvn(n = 1, mu = rep(0, nrow(Sigma_total)), sigma = Sigma_total, ncores = 9)

# Reshape the vector into a matrix with n columns
pheno_col <- matrix(pheno, ncol = 3, byrow = TRUE)

# Create data frame
data = as.data.frame(pheno_col)



mt_tensor_evd <- function(n, grms, n_resp, resp.m = NULL, model, n_pc = 20, formula = NULL, data = NULL){

  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################


  # Constructing the tensor with additive, dominance, and interaction components
  GRMtensor <- array(c(as.matrix(grms$A), as.matrix(grms$D)), #as.matrix(grms$AA)),
                     dim = c(nrow(grms$A), ncol(grms$A), 2))

  # Convert to tensor format
  GRMtensor <- as.Tensor(GRMtensor)

  # Eigen decomposition of the tensor
  teigen <- tEIG(GRMtensor, "dst")

  # Extract matrix D
  A_eigenvalue <- teigen$D$data[1:n_pc,1:n_pc,1]
  D_eigenvalue <- teigen$D$data[1:n_pc,1:n_pc,2]
  #AA_eigenvalue <- teigen$D$data[1:n_pc,1:n_pc,3]

  # Store matrices
  D_A_sparse <- Matrix(A_eigenvalue, sparse = T)
  D_DD_sparse <- Matrix(D_eigenvalue, sparse = T)
  #D_AA_sparse <- Matrix(AA_eigenvalue, sparse = T)


  # Eigenvectors

  Q <- teigen$P$data

  # Reduced
  P <- as.Tensor(Q[,1:n_pc,])

  # Transform response variable matrix

  #trans.resp.m <- as.matrix(resp.m)


  # library(torch)
  #
  # # Convert P and resp.m to torch tensors
  # P_tensor <- torch_tensor(P$data)
  # resp_m_tensor <- torch_tensor(as.matrix(resp.m))
  #
  # # Perform the mode-3 product using Einstein summation
  # trans_resp_m_tensor <- torch_einsum("ijn,im->jm", list(P_tensor, resp_m_tensor))
  #
  # # Convert back to R matrix if needed
  # trans.resp.m <- as.matrix(trans_resp_m_tensor)

  # Einstein Summation
  trans.resp.m <- matrix(0, nrow = n_pc, ncol = ncol(resp.m))

  # Perform the mode-3 product by iterating over the third dimension of P
  for (i in 1:2) {
    # Slice the i-th "layer" of the tensor P (dimension 926 x n_pc)
    P_slice <- P$data[, , i]

    # Perform matrix multiplication (P_slice %*% resp.m)
    trans.resp.m <- trans.resp.m  + t(P_slice) %*% as.matrix(resp.m)
  }


  E <- diag(nrow(trans.resp.m))


  ####################################################################
  ## Extending to multivariate responses #############################
  ####################################################################
  output <- list()

  if (n_resp > 1) {
    Z_struc <- mglm4twin:::mt_struc(n_resp = n_resp)
    ind_A <- lapply(Z_struc, function(x) kronecker(x,  D_A_sparse ))
    ind_D <- lapply(Z_struc, function(x) kronecker(x, D_DD_sparse))
    #ind_AA <- lapply(Z_struc, function(x) kronecker(x, D_AA_sparse))
    ind_E <- lapply(Z_struc, function(x) kronecker(x, E))
  }

  ####################################################################
  ## Selecting the different twin models #############################
  ####################################################################
  if (n_resp > 1) {
    if (model == "AE") {
      output$matrices <- c(ind_A,
                           ind_D,
                           #ind_AA,
                           ind_E)
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


mat <- mt_tensor_evd(grms = c(A = GRMa, D = GRMd),
                     n_resp = 3,
                     resp.m = as.data.frame(pheno_col),
                     n_pc = 900,
                     model = "AE", data = data)


# Prepare for mglm4twin

# Fixed effects

linear_pred_1 <- V1 ~ 1
linear_pred_2 <- V2 ~ 1
linear_pred_3 <- V3 ~ 1

data <- as.data.frame(mat$phenotype)

res <- mglm4twin(linear_pred = c(linear_pred_1, linear_pred_2, linear_pred_3),
                 matrix_pred = c(mat$matrices),
                 data = data)
res$Covariance









##### Matrix Sketching: Random Projections #####

output.snp.col <- Rdimtools::do.rndproj(snp.pine, ndim = 4000, type = "gaussian")
output.snp.row <- Rdimtools::do.rndproj(t(snp.pine), ndim = 100, type = "gaussian")

dim(output.snp.row$Y)
dim(snp.pine)
dim(output.snp.col$Y)


Z = t(output.snp.row$Y) %*% t(snp.pine) %*% output.snp.col$Y

K = Z %*% t(Z)



mt_randproj <- function(n, marker_matrix, n_resp, row_dim, col_dim, resp.m = NULL, model, formula = NULL, data = NULL){

  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################

  # marker_matrix = snp.pine
  # col_dim = 4000
  # row_dim = 300

  # Compute allele frequencies
  p_vec <- colMeans(marker_matrix, na.rm = TRUE) / 2

  # Create copy of marker_matrix
  marker_imputed <- marker_matrix

  # Impute NAs with 2 * p_j (the expected value under Hardy-Weinberg)
  for (j in seq_len(ncol(marker_imputed))) {
    marker_imputed[is.na(marker_imputed[, j]), j] <- 2 * p_vec[j]
  }

  # Center the imputed matrix
  Z <- sweep(marker_imputed, 2, 2 * p_vec, FUN = "-")

  #output.snp.col <- Rdimtools::do.rndproj(Z, ndim = col_dim, type = "gaussian")
  output.snp.row <- Rdimtools::do.rndproj(t(Z), ndim = row_dim, type = "gaussian")

  S1 = t(output.snp.row$projection)
  #S2 = output.snp.col$projection

  dim(S1)
  #dim(S2)
  dim(Z)

  # sketch genotype matrix

  Zs1s2 = S1 %*% Z #%*% S2
  dim(Zs1s2)

  Z1 = t(output.snp.row$projection) %*% marker_matrix
  dim(Z1)

  K = Zs1s2 %*% t(Zs1s2)
  dim(K)

  H = K/ncol(Z)

  A <- as(H, "dsCMatrix")

  # Transform response variable matrix
  trans.resp.m <- t(output.snp.row$projection) %*% as.matrix(resp.m)

  dim(resp.m)

  E <- diag(nrow(trans.resp.m))

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


mat <- mt_randproj(n = 926, marker_matrix = snp.pine, n_resp = 3, row_dim = 100, col_dim = ncol(snp.pine), resp.m = pheno_col, model = "AE", formula = NULL, data = NULL)


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

### sparse svd
### Check sketching package srht


mt_nyström <- function(n, marker_matrix, n_resp, row_dim, resp.m = NULL, model, formula = NULL, data = NULL){

  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################

  # marker_matrix = snp.pine
  # col_dim = 4000
  # row_dim = 300
  # resp.m = pheno_col


  # Compute allele frequencies
  p_vec <- colMeans(marker_matrix, na.rm = TRUE) / 2

  # Create copy of marker_matrix
  marker_imputed <- marker_matrix

  # Impute NAs with 2 * p_j (the expected value under Hardy-Weinberg)
  for (j in seq_len(ncol(marker_imputed))) {
    marker_imputed[is.na(marker_imputed[, j]), j] <- 2 * p_vec[j]
  }

  # Center the imputed matrix
  Z <- sweep(marker_imputed, 2, 2 * p_vec, FUN = "-")

  #output.snp.col <- kfm_nystrom(t(Z), m = nrow(t(Z)), r = 4000, kernel = "radial")
  output.snp.row <- kfm_nystrom(Z, m = nrow(Z), r = row_dim, kernel = "radial")


  S1 = t(output.snp.row$dv)
  S2 = Z

  dim(S1)
  dim(S2)
  dim(Z)

  # sketch genotype matrix

  Zs1s2 = t(S1) %*% Z #%*% t(Z)
  dim(Zs1s2)

  Z1 = t(S1) %*% marker_matrix
  dim(Z1)

  K = Zs1s2 %*% t(Zs1s2)
  dim(K)

  H = K/ncol(Z)

  A <- as(H, "dsCMatrix")

  # Transform response variable matrix
  trans.resp.m <- t(S1) %*% as.matrix(pheno_col)

  dim(trans.resp.m)

  E <- diag(nrow(trans.resp.m))

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



mat <- mt_nyström(n = 926, marker_matrix = snp.pine, n_resp = 3, row_dim = 700, resp.m = pheno_col, model = "AE", formula = NULL, data = NULL)



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



mt_rsvd <- function(n, grm, n_resp, resp.m = NULL, model, formula = NULL, data = NULL){

  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################



  grm_sparse <- Matrix(grm, sparse = FALSE)
  A <- as(grm_sparse, "dsCMatrix")


  res <- rsvd(A)
  d <- res$d
  u <- res$u
  v <- res$v



  # # Eigen decomposition of GRM
  # evd <- eigen(A, symmetric = TRUE)
  # Q <- evd$eigvalsQ <- evd$vectors
  # lambda <- evd$values
  #
  # # Reduce Q
  # n <- length(lambda) - n_pc
  # P <- Q[, 1:n_pc]
  D <- diag(d)


  P <- res$u



  # Transform response variable matrix
  trans.resp.m <- t(P) %*% as.matrix(resp.m)
  E <- diag(nrow(trans.resp.m))


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


mat <- mt_rsvd(grm = GRM, n_resp = 3, resp.m = as.data.frame(pheno_col),  model = "AE", data = data)


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

round(cov_matrix_A, 2)

# Construct covariance matrix
cov_matrix_E <- matrix(c(
  E11, E21, E31,
  E21, E22, E32,
  E31, E32, E33
), nrow = 3, byrow = TRUE)

round(cov_matrix_E, 2)




mt_rpca<- function(n, grm, n_resp, resp.m = NULL, model, formula = NULL, data = NULL){

  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################



  grm_sparse <- Matrix(grm, sparse = FALSE)
  A <- as(grm_sparse, "dsCMatrix")

  grm <- as.matrix(GRM@x)

  res <- rpca::rpca(grm)
  length(res$L.svd$d)
  L <- res$L.svd$d
  L <- res$L.svd$
  S <- res$S


  u <- res$L.svd$u
  vt <- res$L.svd$vt
  L <- res$L.svd$L

  res$L.svd$d



  # # Eigen decomposition of GRM
  # evd <- eigen(A, symmetric = TRUE)
  # Q <- evd$eigvalsQ <- evd$vectors
  # lambda <- evd$values
  #
  # # Reduce Q
  # n <- length(lambda) - n_pc
  # P <- Q[, 1:n_pc]
  D <- diag(d)


  P <- res$u


  # # EigenL# # Eigen decomposition of GRM
  # evd <- eigen(A, symmetric = TRUE)
  # Q <- evd$eigvalsQ <- evd$vectors
  # lambda <- evd$values
  #
  # # Reduce Q
  # n <- length(lambda) - n_pc
  # P <- Q[, 1:n_pc]
  D <- diag(d)


  P <- res$u



  # Transform response variable matrix
  trans.resp.m <- t(P) %*% as.matrix(resp.m)
  E <- diag(nrow(trans.resp.m))


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


mat <- mt_rsvd(grm = GRM, n_resp = 3, resp.m = as.data.frame(pheno_col),  model = "AE", data = data)


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

round(cov_matrix_A, 2)

# Construct covariance matrix
cov_matrix_E <- matrix(c(
  E11, E21, E31,
  E21, E22, E32,
  E31, E32, E33
), nrow = 3, byrow = TRUE)

round(cov_matrix_E, 2)

### Test CUR and Biclustering (Co-clustering)

library(biclust)

erg <- biclust(as.matrix(snp.pine), method=BCCC(), delta=1.5, alpha=1, number=10)
erg

