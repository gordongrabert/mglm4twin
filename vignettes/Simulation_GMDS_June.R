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

#### Load Simulated Data ####

### read data from epigen ###
data <- read_json("/Users/gordonplri/Documents/epigen/epigen/sim/0_1_ASW.json")

### extract and transform ###
snp.name <- data$snps

df <- t(data.frame(matrix(unlist(data$genotype), nrow=length(data$genotype), byrow=TRUE)))

#### Phenotype simulation ####

# Randomly select 100 columns
set.seed(123)  # Setting seed for reproducibility
selected_cols <- sample(ncol(df), 100)
df_selected <- df[, selected_cols]


# Step 1: Define Sample Size & Covariance Matrices
n <- nrow(df) # Number of individuals

#Computing the additive relationship matrix based on VanRaden 2008
GRM <- Gmatrix(SNPmatrix=df, missingValue=-9,
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

# mvrnorm() from MASS package might be considerably slower
#pheno = mvrnorm(n = 1, mu = rep(0, nrow(Sigma_total)), Sigma = Sigma_total)

# you may need to adjust number of cores
pheno = rmvn(n = 1, mu = rep(0, nrow(Sigma_total)), sigma = Sigma_total, ncores = 9)

# Reshape the vector into a matrix with n columns
pheno_col <- matrix(pheno, ncol = 3, byrow = TRUE)

# Create data frame
data = as.data.frame(pheno_col)


