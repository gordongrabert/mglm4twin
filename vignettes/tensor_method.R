#### Test Tensor Method ####

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
library(copula)

# Required packages
library(Matrix)
library(abind)
library(rTensor)
library(MASS)      # for mvrnorm
library(tensorA)   # optional for unfolding/rotation utilities


#### Load Simulated Data ####

### read data from epigen ###
data <- read_json("/Users/gordonplri/Documents/epigen/epigen/sim/0_1_ASW.json")

### data extraction

df <- t(data.frame(matrix(unlist(data$genotype), nrow=length(data$genotype), byrow=TRUE)))

n = 2000

## 1. Subset data
selected_rows <- sample(nrow(df), n)
df_selected <- df[selected_rows, ]

calcGRM <- function(genoMat,
                    methodGRM = "addNOIA",
                    subpop = NULL,
                    kernel.h = "tuned",
                    returnWMat = FALSE,
                    probaa = NULL,
                    probAa = NULL,
                    batchSize = NULL,
                    n.core = 1) {
  supportedMethods <- c("addNOIA", "domNOIA", "A.mat", "linear",
                        "gaussian", "exponential", "correlation")
  stopifnot(methodGRM %in% supportedMethods)
  
  genoMatUniq <- sort(unique(c(genoMat)), decreasing = FALSE)
  genoMatUniqLen <- length(genoMatUniq)
  
  if (genoMatUniqLen == 2) {
    isScoring1 <- all(genoMatUniq == c(-1, 1)) | all(genoMatUniq == c(-1, 0))
    isScoring2 <- all(genoMatUniq == c(0, 2)) | all(genoMatUniq == c(0, 1))
  } else {
    if (genoMatUniqLen == 3) {
      isScoring1 <- all(genoMatUniq == c(-1, 0, 1))
      isScoring2 <- all(genoMatUniq == c(0, 1, 2))
    } else {
      stop("Something wrong with your genotype data!!")
    }
  }
  
  if (isScoring1) {
    genoMat <- genoMat
  } else {
    if (isScoring2) {
      genoMat <- genoMat - 1
    } else {
      stop("Genotype data should be scored with (-1, 0, 1) or (0, 1, 2)!!")
    }
  }
  
  nInd <- nrow(genoMat)
  nMarkers <- ncol(genoMat)
  mrkNames <- colnames(genoMat)
  
  methodNOIA <- stringr::str_detect(string = methodGRM,
                                    pattern = "NOIA")
  if ((!methodNOIA) & (!is.null(subpop))) {
    message("`subpop` information is utilized only when you use `NOIA` methods ('addNOIA' & 'domNOIA') !")
  }
  if (methodNOIA) {
    if (is.null(subpop)) {
      if (is.null(probaa)) {
        probaa <- apply(genoMat == -1, 2, mean)
      }
      if (is.null(probAa)) {
        probAa <- apply(genoMat == 0, 2, mean)
      }
      if (methodGRM == "addNOIA") {
        replaceaa <- - (2 - probAa - 2 * probaa)
        replaceAa <- - (1 - probAa - 2 * probaa)
        replaceAA <- - (- probAa - 2 * probaa)
      } else if (methodGRM == "domNOIA") {
        probAA <- 1 - probaa - probAa
        denominator <- probAA + probaa - (probAA - probaa) ^ 2
        replaceaa <- - 2 * probAA * probAa / denominator
        replaceAa <- 4 * probAA * probaa / denominator
        replaceAA <- - 2 * probaa * probAa / denominator
      }
      
      HMat <- sapply(1:nMarkers, function(mrkNo) {
        HMatEachMrk <- genoMat[, mrkNo]
        HMatEachMrk[HMatEachMrk == -1] <- replaceaa[mrkNo]
        HMatEachMrk[HMatEachMrk == 0] <- replaceAa[mrkNo]
        HMatEachMrk[HMatEachMrk == 1] <- replaceAA[mrkNo]
        
        return(HMatEachMrk)
      })
      colnames(HMat) <- mrkNames
    } else {
      stopifnot(length(subpop) == nInd)
      
      probaa <- do.call(
        what = rbind,
        args = sapply(X = unique(subpop),
                      FUN = function(subpopEach) {
                        genoMatSubpop <- genoMat[subpop %in% subpopEach, , drop = FALSE]
                        probaaSubpop <- apply(genoMatSubpop == -1, 2, mean)
                      },
                      simplify = FALSE)
      )[subpop, , drop = FALSE]
      rownames(probaa) <- rownames(genoMat)
      
      probAa <- do.call(
        what = rbind,
        args = sapply(X = unique(subpop),
                      FUN = function(subpopEach) {
                        genoMatSubpop <- genoMat[subpop %in% subpopEach, , drop = FALSE]
                        probAaSubpop <- apply(genoMatSubpop == 0, 2, mean)
                      },
                      simplify = FALSE)
      )[subpop, , drop = FALSE]
      rownames(probAa) <- rownames(genoMat)
      
      
      if (methodGRM == "addNOIA") {
        replaceaa <- - (2 - probAa - 2 * probaa)
        replaceAa <- - (1 - probAa - 2 * probaa)
        replaceAA <- - (- probAa - 2 * probaa)
      } else if (methodGRM == "domNOIA") {
        probAA <- 1 - probaa - probAa
        denominator <- probAA + probaa - (probAA - probaa) ^ 2
        replaceaa <- - 2 * probAA * probAa / denominator
        replaceAa <- 4 * probAA * probaa / denominator
        replaceAA <- - 2 * probaa * probAa / denominator
      }
      
      HMat <- genoMat
      HMat[HMat == -1] <- replaceaa[HMat == -1]
      HMat[HMat == 0] <- replaceAa[HMat == 0]
      HMat[HMat == 1] <- replaceAA[HMat == 1]
    }
    
    if (is.null(batchSize)) {
      HHt <- tcrossprod(HMat)
    } else {
      batchIds <- 1:ncol(HMat) %/% batchSize + 1
      
      HHtList <- parallel.compute(
        vec = 1:max(batchIds),
        func = function(batchNo) {
          HHtBatch <- tcrossprod(HMat[, batchIds == batchNo])
          
          return(HHtBatch)
        },
        n.core = n.core,
        count = FALSE
      )
      
      HHt <- Reduce(
        f = `+`,
        x = HHtList
      )
    }
    
    GRM <- HHt * nInd / sum(diag(HHt))
  } else if (methodGRM == "A.mat") {
    GRM <- rrBLUP::A.mat(X = genoMat)
  } else if (methodGRM == "linear") {
    if (is.null(batchSize)) {
      HHt <- tcrossprod(genoMat)
    } else {
      batchIds <- 1:ncol(genoMat) %/% batchSize + 1
      
      if (n.core == 1) {
        HHtList <- lapply(
          X = 1:max(batchIds),
          FUN = function(batchNo) {
            HHtBatch <- tcrossprod(genoMat[, batchIds == batchNo])
            
            return(HHtBatch)
          }
        )
      } else {
        HHtList <- parallel::mclapply(
          X = 1:max(batchIds),
          FUN = function(batchNo) {
            HHtBatch <- tcrossprod(genoMat[, batchIds == batchNo])
            
            return(HHtBatch)
          },
          mc.cores = n.core
        )
      }
      HHt <- Reduce(
        f = `+`,
        x = HHtList
      )
    }
    
    GRM <- HHt * nInd / sum(diag(HHt))
  } else if (methodGRM == "gaussian") {
    distMat <- Rfast::Dist(x = genoMat) / sqrt(ncol(genoMat))
    rownames(distMat) <- colnames(distMat) <- rownames(genoMat)
    if ("character" %in% class(kernel.h)) {
      hinv <- median((distMat ^ 2)[upper.tri(distMat ^ 2)])
      h <- 1 / hinv
    } else if ("numeric" %in% class(kernel.h)) {
      h <- kernel.h
    }
    
    GRM <- exp(- h * distMat ^ 2)
  } else if (methodGRM == "exponential") {
    distMat <- Rfast::Dist(x = genoMat) / sqrt(ncol(genoMat))
    rownames(distMat) <- colnames(distMat) <- rownames(genoMat)
    if ("character" %in% class(kernel.h)) {
      hinv <- median((distMat ^ 2)[upper.tri(distMat ^ 2)])
      h <- 1 / hinv
    } else if ("numeric" %in% class(kernel.h)) {
      h <- kernel.h
    }
    
    GRM <- exp(- h * distMat)
  } else if (methodGRM == "correlation") {
    GRM <- cor(t(genoMat))
  }
  
  
  
  if (methodNOIA & returnWMat) {
    WMat <- HMat * sqrt(nInd / sum(HMat * HMat))
    rownames(WMat) <- rownames(genoMat)
    return(WMat)
  } else {
    rownames(GRM) <- colnames(GRM) <- rownames(genoMat)
    return(GRM)
  }
}


## recode genotypes

recode_genotypes <- function(genoMat, check_valid = TRUE, maf_threshold = NULL) {
  genoMat <- as.matrix(genoMat)
  
  # Validate input values
  if (check_valid) {
    valid_vals <- na.omit(as.vector(genoMat))
    if (!all(valid_vals %in% c(0, 1, 2))) {
      stop("Input contains values outside the expected 0/1/2 range.")
    }
  }
  
  # Recode 0/1/2 to -1/0/1
  recoded <- ifelse(is.na(genoMat), NA, genoMat - 1)
  
  # Apply MAF filter if requested
  if (!is.null(maf_threshold)) {
    maf <- colMeans(genoMat, na.rm = TRUE) / 2  # MAF based on 0/1/2 coding
    maf <- pmin(maf, 1 - maf)
    keep <- maf >= maf_threshold
    recoded <- recoded[, keep, drop = FALSE]
  }
  
  return(recoded)
}

## genoMat_noia 

# genoMat_noia <- recode_genotypes(df_selected, maf_threshold = 0.05)
# 
# A <- calcGRM(genoMat = genoMat_noia, methodGRM = "addNOIA")
# #A <- nearPD(A)$mat
# D <- calcGRM(genoMat = genoMat_noia, methodGRM = "domNOIA")
# #D <- nearPD(D)$mat
# #AD <- A*D

# ## 2. Compute GRM and make PD
 A <- Gmatrix(SNPmatrix = df_selected, missingValue = -9,
                maf = 0.05, method = "VanRaden")
 A <- nearPD(A)$mat
 D <- Gmatrix(SNPmatrix = df_selected, missingValue = -9,
              maf = 0.05, method = "Vitezica")
 D <- nearPD(D)$mat
I <- diag(n)




Sigma_total <- 0.5*A + 0.3*D + 0.2*I



pheno = as.data.frame(t(rmvn(n = 1, mu = rep(0, nrow(Sigma_total)), sigma = Sigma_total, ncores = 9)))


# Prepare for mglm4twin

# Fixed effects

# When there are no fixed effects

linear_pred_1 <- V1 ~ 1

A_sparse <- Matrix(A, sparse = F)
A_dsC <- as(A_sparse, "dsCMatrix")
D_sparse <- Matrix(D, sparse = F)
D_dsC <- as(D_sparse, "dsCMatrix")
I <- as(I, "dtCMatrix")
Mat <- list(I = I, A = A_dsC, D = D_dsC)

res <- mglm4twin(linear_pred = c(linear_pred_1),
                 matrix_pred = c(Mat),
                 data = as.data.frame(pheno))

#### Joint Diagonalization of Real Matrices ###



## 1) D whiten (Ridge + Trunkierung für Stabilität)
De <- eigen(D_dense, symmetric = TRUE)
s  <- pmax(De$values, 0)
tau <- 1e-10 * max(s)
keep <- s > tau
U  <- De$vectors[, keep, drop = FALSE]
s  <- s[keep]
W  <- U %*% diag(1/sqrt(s)) %*% t(U)          # D^{-1/2} auf dem Trägerraum

## 2) A im D-whitened Raum diagonalisieren
B <- crossprod(W, A_dense %*% W)                   # W^T A W
Be <- eigen(B, symmetric = TRUE)
Q <- Be$vectors
Lambda <- Be$values

## 3) Gesamttransformation
V <- W %*% Q                                  # liefert: V^T D0 V = I, V^T A0 V = diag(Lambda)

## 4) Rotierte Größen/Kerne
pheno$V1.rot <- as.vector(crossprod(V, pheno$V1))         # y*
A_diag <- Diagonal(x = Lambda)                # A' = diag(Lambda)
D_I    <- Diagonal(x = rep(1, length(Lambda)))# D' = I
E_star <- t(V) %*% V                       # E' = V^T V  (Residuum in rotem Raum)

## 5) In ein LMM/REML geben
# Falls dein Fitter diagonale + allgemeine Kerne akzeptiert:
Mat <- list( E = as(E_star, "dgCMatrix"),
            A = as(A_diag, "dtCMatrix"),
            D = as(D_I, "dtCMatrix"))


linear_pred_1 <- V1.rot ~ 1


res <- mglm4twin(linear_pred = c(linear_pred_1),
                 matrix_pred = c(Mat),
                 data = as.data.frame(pheno))
res


###### Truncated HOSVD ###

library(abind)
library(rTensor)
library(Matrix)
library(ggplot2)

# 1. Combine A and D into a tensor of shape [n_individuals x n_individuals x 2]
g_tensor <- abind(as.matrix(A), as.matrix(D), along = 3) %>% as.tensor()

# 2. Perform Tucker decomposition with specified low rank
r <- 2000  # initial rank for truncation
tensor_decomp <- tucker(g_tensor, ranks = c(r, r, 2))

# 3. Extract core and factor matrices
core <- tensor_decomp$Z
U1 <- tensor_decomp$U[[1]]
U2 <- tensor_decomp$U[[2]]
U3 <- tensor_decomp$U[[3]]

# 4. Truncate to a smaller core for further approximation
r_trunc <- 1000
core_trunc <- core[1:r_trunc, 1:r_trunc, ]
U1_trunc <- U1[, 1:r_trunc]
U2_trunc <- U2[, 1:r_trunc]

# Apply soft-thresholding to induce sparsity in the truncated core
soft_thresh <- function(x, lambda) sign(x) * pmax(abs(x) - lambda, 0)
core_tensor_sparse <- soft_thresh(core_trunc@data, lambda = 1e-2)

# Extract absolute values of the slices
lowrank_A <- as.matrix(abs(core_tensor_sparse[,,1]))
lowrank_D <- as.matrix(abs(core_tensor_sparse[,,2]))

# Convert to sparse symmetric matrices
A_mat <- as(forceSymmetric(Matrix(lowrank_A, sparse = TRUE)), "dsCMatrix")
D_mat <- as(forceSymmetric(Matrix(lowrank_D, sparse = TRUE)), "dsCMatrix")
I_mat <- Diagonal(n = nrow(A_mat))

# Extract only diagonals from the core slices
A_diag <- diag(core_trunc@data[,,1])
D_diag <- diag(core_trunc@data[,,2])

# Construct diagonal matrices and convert to sparse symmetric format
A_core_diag <- abs(diag(A_diag))
D_core_diag <- abs(diag(D_diag))

A_mat <- as(forceSymmetric(Matrix(A_core_diag, sparse = TRUE)), "dsCMatrix")
D_mat <- as(forceSymmetric(Matrix(D_core_diag, sparse = TRUE)), "dsCMatrix")
I_mat <- Diagonal(n = nrow(A_mat))

# Visualize the first two principal components from U1
tensor.pc <- as.data.frame(U1)
ggplot(tensor.pc, aes(x = V1, y = V2)) +
  geom_point(alpha = 0.7) +
  labs(title = "Mode-1 Tucker Components",
       x = "PC 1", y = "PC 2") +
  theme_minimal()

# Project phenotype matrix (pheno) into individual component space
P_rot <- as.data.frame(t(U1_trunc) %*% as.matrix(pheno))
colnames(P_rot) <- paste0("V", seq_len(ncol(P_rot)))


# Define formula (can include more PCs)
linear_pred_1 <- V1 ~ 1

# Fit multivariate twin model using the rotated phenotype data
res <- mglm4twin(
  linear_pred = list(linear_pred_1),
  matrix_pred = list(I = I_mat, A = A_mat, D = D_mat),
  data = P_rot
)

summary <- summary(res, model = "ADE", biometric = "TRUE")

disp <- summary$Dispersion


#### Loop ####
library(mglm4twin)
library(Matrix)

results_list <- list()
dispersion_summary <- data.frame(
  n = integer(),
  dispersion = numeric(),
  time_tucker = numeric(),
  time_mglm4twin = numeric()
)

# Initialize storage outside the loop
dispersion_summary <- data.frame()
results_list <- list()

for (n in seq(500, 3000, by = 500)) {
  cat("Running sample size n =", n, "\n")
  
  n = 2500
  
  # 1. Subset and prepare genotype data
  selected_rows <- sample(nrow(df), n)
  df_selected <- df[selected_rows, ]
  genoMat_noia <- recode_genotypes(df_selected, maf_threshold = 0.05)
  
  # 2. Compute GRMs
  A <- calcGRM(genoMat = genoMat_noia, methodGRM = "addNOIA")
  D <- calcGRM(genoMat = genoMat_noia, methodGRM = "domNOIA")
  I <- diag(n)
  
  # 3. Create total covariance (true weights: 0.5 A + 0.3 D + 0.2 I)
  Sigma_total <- 0.5 * A + 0.3 * D + 0.2 * I
  
  # 4. Simulate phenotype (one replicate, n-dim vector)
  pheno <- as.data.frame(t(rmvn(n = 1, mu = rep(0, n), sigma = Sigma_total, ncores = 9)))
  
  # 5. Stack GRMs into tensor and do Tucker decomposition
  g_tensor <- abind(as.matrix(A), as.matrix(D), along = 3) %>% rTensor::as.tensor()
  
  time_tucker <- system.time({
    tensor_decomp <- tucker(g_tensor, ranks = c(n, n, 2))
  })["elapsed"]
  
  # Extract core and factor matrices
  core <- tensor_decomp$Z
  U1 <- tensor_decomp$U[[1]]
  U2 <- tensor_decomp$U[[2]]
  U3 <- tensor_decomp$U[[3]]
  

  # Extract diagonals of core slices and create sparse diagonal matrices
  A_diag <- diag(core@data[, , 1])
  D_diag <- diag(core@data[, , 2])
  
  A_mat <- as(forceSymmetric(Matrix(diag(A_diag), sparse = TRUE)), "dsCMatrix")
  D_mat <- as(forceSymmetric(Matrix(diag(D_diag), sparse = TRUE)), "dsCMatrix")
  I_mat <- Diagonal(n = nrow(A_mat))
  
  # Project phenotype onto factor U1
  P_rot <- as.data.frame(t(U1) %*% as.matrix(pheno))
  colnames(P_rot) <- paste0("V", seq_len(ncol(P_rot)))
  
  # Fit mglm4twin model (timed)
  time_mglm4twin <- system.time({
    linear_pred_1 <- V1 ~ 1
    
    fit <- mglm4twin(
      linear_pred = list(linear_pred_1),
      matrix_pred = list(I = I_mat,A = A_mat, D = D_mat),
      data = P_rot
    )
  })["elapsed"]
  
  # Extract dispersion and store results
  fit_summary <- summary(fit, model = "ADE", biometric = TRUE)
  dispersion_value <- fit_summary$Dispersion
  
  dispersion_summary <- rbind(dispersion_summary, data.frame(
    n = n,
    dispersion = dispersion_value,
    time_tucker = time_tucker,
    time_mglm4twin = time_mglm4twin
  ))
  
  results_list[[paste0("n_", n)]] <- list(
    fit = fit,
    summary = fit_summary,
    time_tucker = time_tucker,
    time_mglm4twin = time_mglm4twin
  )
}


