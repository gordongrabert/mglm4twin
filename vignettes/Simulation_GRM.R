## Loading extra packages
require(Matrix)
require(SimCorMultRes)
require(mglm4twin)
require(ggplot2)
require(ggExtra)
require(AGHmatrix)
require(MASS)


#### Use SNP data from AGHmatrix #####

data(snp.pine)

# Step 1: Define Sample Size & Covariance Matrices
n <- 926 # Number of individuals

#Computing the additive relationship matrix based on VanRaden 2008
GRM <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9,
               maf=0.05, method="VanRaden")


# Genetic Covariance Matrix (Sigma_G)
Sigma_G <- matrix(c(0.5, 0.5, 0.3,  # Trait 1 variance & covariance
                    0.5, 0.3, 0.1,
                    0.3, 0.1, 0.2), # Trait 2 variance & covariance
                  nrow = 3, byrow = TRUE)

# Environmental Covariance Matrix (Sigma_E)
Sigma_E <- matrix(c(0.3, 0.2, 0.6,  # Trait 1 variance & covariance
                    0.2, 0.4, 0.1,
                    0.6, 0.1, 0.3), # Trait 2 variance & covariance
                  nrow = 3, byrow = TRUE)

# Step 2: Simulate from Gaussian Copula
# Construct full covariance matrix for copula

# Step 2: Construct the Full Covariance Matrix
Sigma_total <- as.matrix(kronecker(GRM, Sigma_G) + kronecker(diag(n), Sigma_E))

# Search next PD matrix
Sigma_total <- nearPD(Sigma_total)$mat

# Convert to a correlation matrix
# Sigma_total <- cov2cor(Sigma_total)

# Simulate Gaussian data
pheno = mvrnorm(n = 1, mu = rep(0, nrow(Sigma_total)), Sigma = Sigma_total)

# Reshape the vector into a matrix with 2 columns (will create pairs)
pheno_col <- matrix(pheno, ncol = 3, byrow = TRUE)

# Create data frame
data = as.data.frame(pheno_col)

# Prepare for mglm4twin
V1 = data$V1
V2 = data$V2

linear_pred_1 <- V1 ~ 1
linear_pred_2 <- V2 ~ 1
linear_pred_3 <- V3 ~ 1

n <- 926 # Number of individuals

mat <- mt_grm(n = n, grm = GRM, n_resp = 3, model = "AE", data = data)

res <- mglm4twin(linear_pred = c(linear_pred_1, linear_pred_2, linear_pred_3),
                 matrix_pred = c(mat),
                 data = as.data.frame(data))

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


