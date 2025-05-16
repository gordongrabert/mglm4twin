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

#### Load Simulated Data ####

### read data from epigen ###
data <- read_json("/Users/gordonplri/Documents/epigen/epigen/sim/0_1_ASW.json")

### extract and transform ###
snp.name <- data$snps

df <- t(data.frame(matrix(unlist(data$genotype), nrow=length(data$genotype), byrow=TRUE)))

#### Phenotype simulation ####

# Randomly select 1000 rows
set.seed(123)
selected_rows <- sample(nrow(df), 1000)
df_selected <- df[selected_rows, ]  # Rows selected here

# Step 1: Define Sample Size & Covariance Matrices
n <- nrow(df_selected) # Number of individuals

#Computing the additive relationship matrix based on VanRaden 2008
GRM <- Gmatrix(SNPmatrix=df_selected, missingValue=-9,
               maf=0.05, method="VanRaden")


# Search next PD matrix
GRM <- nearPD(GRM)$mat

eigen <- eigen(GRM)
plot(eigen$values)
title(main = "Eigenvalues GRM")


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


# you may need to adjust number of cores
pheno = rmvn(n = 1, mu = rep(0, nrow(Sigma_total)), sigma = Sigma_total, ncores = 9)

# Reshape the vector into a matrix with n columns
pheno_col <- matrix(pheno, ncol = 3, byrow = TRUE)

# Compute correlation matrices
cor_pheno <- cor(pheno_col)

# Print the results
print(cor_pheno)

# Create data frame
data = as.data.frame(pheno_col)


#### Wagner Paper: Multi-bounded data


## Setting model parameters
A <- c(0.25,0.3, 0.35, -0.15, 0.20, -0.2)
E <- c(0.75, 0.7, 0.65, -0.3, 0.25, -0.4)
tau = c(A, E)


## GRM structure
mat <- mt_grm(n = n, grm = GRM, n_resp = 3, model = "AE", data = data)
Omega <- as.matrix(mt_matrix_linear_predictor(tau = tau, Z = mat))

# Create correlation matrix
Omega.cor <- cov2cor(Omega)

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

for (i in 1:n) {
  Y <- rnorta(R = 1, cor.matrix = Omega.cor,
                distr = invcdfnames, qparameters = qparameters)
    Y1[[i]] <- Y[1]
    Y2[[i]] <- Y[2]
    Y3[[i]] <- Y[3]
}




Y1 <- c(do.call(c, Y1_DZ), do.call(c, Y1_MZ))
Y2 <- c(do.call(c, Y2_DZ), do.call(c, Y2_MZ))
Y3 <- c(do.call(c, Y3_DZ), do.call(c, Y3_MZ))

data <- data.frame("Y1" = Y1, "Y2" = Y2, "Y3" = Y3, "twin_id" = rep(1:2, 414),
                   "zyg" = rep(zyg, each = 2), "sex" = rep(sex, each = 2),
                   "age_std" = rep(age_std, each = 2))




