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
data("snp.pine")
### extract and transform ###
snp.name <- data$snps

df <- t(data.frame(matrix(unlist(data$genotype), nrow=length(data$genotype), byrow=TRUE)))

#### Phenotype simulation ####

# Randomly select n rows
set.seed(123)
selected_rows <- sample(nrow(df), 2000)
df_selected <- df[selected_rows, ]  # Rows selected here

# Step 1: Define Sample Size & Covariance Matrices
n <- nrow(df_selected) # Number of individuals

#Computing the additive relationship matrix based on VanRaden 2008
GRM <- Gmatrix(SNPmatrix=df_selected, missingValue=-9,
               maf=0.05, method="VanRaden")

n = 2000

# Search next PD matrix
GRM <- nearPD(GRM)$mat

eigen <- eigen(GRM)
plot(eigen$values)
title(main = "Eigenvalues GRM")



#### Wagner Paper: Multi-bounded data


## Setting model parameters
## Setting model parameters
E <- c(0.75, 0.7, 0.65, -0.3, 0.25, -0.4)
A <- c(0.25,0.3, 0.35, -0.15, 0.20, -0.2)

### For check 0.3 + 0.15 = 0.45 = 100% 
tau = c(E, A)


## GRM structure

GRM.2 <- cov2cor(GRM)

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

mat <- mt_grm_1(n = n, grm = GRM.2, n_resp = 3, model = "AE", data = NULL)
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

Y1 <- Y[1:2000]
Y2 <- Y[2001:4000]
Y3 <- Y[4001:6000]


data <- data.frame("Y1" = Y1, "Y2" =  Y2 , "Y3" = Y3,
                   "trt" = trt, "sex" = sex,
                   "age_std" = age_std)

hist(data$Y1)
hist(data$Y2)
hist(data$Y3)


my_cols <- c("#00AFBB", "#E7B800" )
pairs(data[,1:3], pch = 19,  cex = 0.5,
      lower.panel=NULL, col = my_cols[as.factor(data$trt)])


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

data.model <- as.data.frame(mat$data)


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

summary(res, model = "AE", biometric = T)

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


#### Hybrid Simulation #####


#### Use SNP data from AGHmatrix #####

data(snp.pine)

# Step 1: Define Sample Size & Covariance Matrices
n <- 926 # Number of individuals

#Computing the additive relationship matrix based on VanRaden 2008
GRM <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9,
               maf=0.05, method="VanRaden")

# Search next PD matrix
GRM <- nearPD(GRM)$mat

# Genetic Covariance Matrix (Sigma_G)
Sigma_G <- matrix(c(0.5, 0.3, 0.2,  # Trait 1 variance & covariance
                    0.3, 0.5, 0.25, # Trait 2 variance & covariance
                    0.2, 0.25, 0.4),# Trait 3 variance & covariance
                  nrow = 3, byrow = TRUE)


# Environmental Covariance Matrix (Sigma_E)
Sigma_E <- matrix(c(0.5, 0.2, 0.3,  # Trait 1 variance & covariance
                    0.2, 0.5, 0.15, # Trait 2 variance & covariance
                    0.3, 0.15, 0.6),# Trait 3 variance & covariance
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

### Bounded dristibutiom

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

# Generate mu values directly (bounded in (0,1))
mu1 <- rbeta(1, shape1 = 10, shape2 = 10)        # symmetric around 0.5
mu2 <- rbeta(1, shape1 = 5, shape2 = 2)        # skewed toward 1
mu3 <- rbeta(1, shape1 = 2, shape2 = 5)        # skewed toward 0
mu1 <- rep(mu1,n)
mu2 <- rep(mu2,n)
mu3 <- rep(mu3,n)


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


Y <- rnorta(R = 1, cor.matrix = cov2cor(Sigma_total),
            distr = invcdfnames, qparameters = qparameters)

Y1 <- Y[1:926]
Y2 <- Y[927:1852]
Y3 <- Y[1853:2778]


data <- data.frame("Y1" = Y1, "Y2" =  Y2 , "Y3" = Y3)

hist(data$Y1)
hist(data$Y2)
hist(data$Y3)


my_cols <- c("#00AFBB", "#E7B800" )
pairs(data[,1:3], pch = 19,  cex = 0.5,
      lower.panel=NULL)


mat <- mt_grm_1(n=1000, grm = GRM, n_resp = 3, model = "AE", data = NULL)

mat <- mt_rsvd(n=926, grm = GRM, n_resp = 3, model = "AE", data = data)

data.model <- as.data.frame(mat$data)


form_Y1 <- c(Y1 ~ 1)
form_Y2 <- c(Y2 ~ 1)
form_Y3 <- c(Y3 ~ 1)

link = rep("logit", 3)
variance = rep("binomialP", 3)



res <- mglm4twin(linear_pred = c(form_Y1, form_Y2, form_Y3),
                 matrix_pred = c(mat$matrices),
                 link = link,
                 variance = variance,
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




