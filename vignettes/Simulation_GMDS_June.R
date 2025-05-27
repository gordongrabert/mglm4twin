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
library(copula)

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
## Groundtruth

# Reorder E and A for lower triangle filling:
# Order needed: E1, E12, E2, E13, E23, E3
E_lt <- c(E[1], E[4], E[5], E[2], E[6], E[3])
A_lt <- c(A[1], A[4], A[5], A[2], A[6], A[3])

# Heritability 
h2_vals <- A_lt / (A_lt + E_lt)

h2_matrix <- matrix(0, 3, 3)
h2_matrix[lower.tri(h2_matrix, diag = TRUE)] <- h2_vals
h2_matrix <- h2_matrix + t(h2_matrix) - diag(diag(h2_matrix))

# Environmentality
e2_vals <- E_lt / (A_lt + E_lt)

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

#Age <- rbeta(n, shape1 = 0.3*2, shape2 = 0.7*2)*20 + 70 ## from Wagner
Age <- rnorm(n, mean = 70, sd = 10)
# Optional: truncate values to keep them within a reasonable range (e.g., 50 to 90)
Age <- pmin(pmax(Age, 50), 90)
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
hist(data$age_std)


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
  
  
  
  # Inverse logit
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  # Safe logit
  safe_logit <- function(p) {
    eps <- 1e-6
    log(pmin(pmax(p, eps), 1 - eps) / (1 - pmin(pmax(p, eps), 1 - eps)))
  }
  
  # Detect numeric & probabilistic columns
  data_numeric <- data[, sapply(data, is.numeric)]
  is_prob_column <- function(x) all(x > 0 & x < 1)
  
  data_prob    <- data_numeric[, sapply(data_numeric, is_prob_column)]
  data_nonprob <- data_numeric[, !sapply(data_numeric, is_prob_column)]
  
  # Transform
  data_logit <- as.data.frame(lapply(data_prob, safe_logit))
  data_nonprob_scaled <- scale(data_nonprob)
  data_transformed <- cbind(data_logit, data_nonprob_scaled)
  
  # GRM decomposition
  grm_sparse <- Matrix(grm, sparse = FALSE)
  res <- eigen(grm_sparse)
  A <- diag(res$values)
  Q <- as.matrix(res$vectors)
  
  # Projection
  P <- t(Q) %*% as.matrix(data_transformed)
  P <- as.data.frame(P)
  
  # Inverse-transform probabilistic variables
  P[ colnames(data_logit) ] <- lapply(P[ colnames(data_logit) ], inv_logit)
  colnames(P) <- colnames(data_numeric)
  P$sex <- data$sex

  # 
  # inv_logit <- function(x) 1 / (1 + exp(-x))
  # data_numeric <- data[, sapply(data, is.numeric)]
  # 
  # grm_sparse <- Matrix(grm, sparse = FALSE)
  # A <- as(grm_sparse, "dsCMatrix")
  # 
  # res <- eigen(A)
  # 
  # A <- diag(res$values)
  # 
  # Q <- as.matrix(res$vectors)
  # 
  # is_prob_column <- function(x) all(x > 0 & x < 1)
  # 
  # data_prob   <- data_numeric[, sapply(data_numeric, is_prob_column)]
  # data_nonprob <- data_numeric[, !sapply(data_numeric, is_prob_column)]
  # 
  # # Apply logit only where valid
  # data_logit <- log(data_prob / (1 - data_prob))
  # # Combine again
  # data_transformed <- cbind(data_logit, data_nonprob)
  # 
  # 
  # # Project
  # P <- t(Q) %*% as.matrix(data_transformed)
  # P <- as.data.frame(P)
  # P[1:3] <- inv_logit(P[1:3])
  # colnames(P) <- colnames(data_numeric)
  # P$sex <- data$sex

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
mt_copula <- function(n, grm, n_resp, model, formula = NULL, data = NULL){
  
  ####################################################################
  ## Prepare matrices ################################################
  ####################################################################
  
  
  # Detect numeric & probabilistic columns
  data_numeric <- data[, sapply(data, is.numeric)]
  is_prob_column <- function(x) all(x > 0 & x < 1)
  
  data_prob    <- data_numeric[, sapply(data_numeric, is_prob_column)]
  data_nonprob <- data_numeric[, !sapply(data_numeric, is_prob_column)]
  
  # Step 2: Gaussian Copula transform function
  gaussian_copula_transform <- function(x) {
    u <- rank(x, ties.method = "average") / (length(x) + 1)
    qnorm(u)
  }
  
  
  # Transform
  # data_copula<- as.data.frame(lapply(data_prob, gaussian_copula_transform))
  # data_nonprob_scaled <- scale(data_nonprob)
  # data_transformed <- cbind(data_copula, data_nonprob_scaled)
  
  data_copula <- as.data.frame(lapply(data_numeric, gaussian_copula_transform))
  
  # GRM decomposition
  grm_sparse <- Matrix(grm, sparse = FALSE)
  res <- eigen(grm_sparse)
  A <- diag(res$values)
  Q <- as.matrix(res$vectors)
  
  # Projection
  P <- t(Q) %*% as.matrix(data_copula)
  P <- as.data.frame(P)
  colnames(P) <- colnames(data_numeric)
  P$sex <- data$sex
  
  # 
  # inv_logit <- function(x) 1 / (1 + exp(-x))
  # data_numeric <- data[, sapply(data, is.numeric)]
  # 
  # grm_sparse <- Matrix(grm, sparse = FALSE)
  # A <- as(grm_sparse, "dsCMatrix")
  # 
  # res <- eigen(A)
  # 
  # A <- diag(res$values)
  # 
  # Q <- as.matrix(res$vectors)
  # 
  # is_prob_column <- function(x) all(x > 0 & x < 1)
  # 
  # data_prob   <- data_numeric[, sapply(data_numeric, is_prob_column)]
  # data_nonprob <- data_numeric[, !sapply(data_numeric, is_prob_column)]
  # 
  # # Apply logit only where valid
  # data_logit <- log(data_prob / (1 - data_prob))
  # # Combine again
  # data_transformed <- cbind(data_logit, data_nonprob)
  # 
  # 
  # # Project
  # P <- t(Q) %*% as.matrix(data_transformed)
  # P <- as.data.frame(P)
  # P[1:3] <- inv_logit(P[1:3])
  # colnames(P) <- colnames(data_numeric)
  # P$sex <- data$sex
  
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
mt_copula_2 <- function(n, grm, n_resp, model, formula = NULL, data = NULL, marginals = NULL, backtransform = FALSE) {
  library(fitdistrplus)
  library(copula)
  library(Matrix)
  library(sn)  # for skew-normal
  
  

# 
#   marginals <- list(
#     Y1 = "beta",
#     Y2 = "beta",
#     Y3 = "beta",
#     age_std = "sn"  # skew-normal
#   )
# 
# 
# 
#     n = nrow(data_select)
#     grm = GRM
#     n_resp = 3
#     model = "AE"  # or "E", depending on what you want
#     formula = NULL
#     data = data_select
#     marginals = marginals
#     backtransform = TRUE  # or FALSE if you want to stay in the copula-transformed space

  
    
      
      # Detect numeric columns and separate sex
      data_numeric <- data[, sapply(data, is.numeric)]
      sex <- data$sex
      
      # Fit marginals if not provided
      if (is.null(marginals)) {
        marginals <- lapply(data_numeric, function(x) {
          if (all(x > 0 & x < 1)) {
            "beta"
          } else if (all(x >= 0 & floor(x) == x)) {
            "pois"
          } else if (abs(skewness(x)) > 1) {
            "sn"
          } else {
            "norm"
          }
        })
      }
      
      # Gaussian copula transformation
      gaussian_copula_transform <- function(x, dist_name) {
        fit <- switch(dist_name,
                      beta = fitdist(x, "beta", start = list(shape1 = 1, shape2 = 1)),
                      sn   = selm(x ~ 1, family = "SN"),
                      norm = fitdist(x, "norm"),
                      stop("Unsupported distribution: ", dist_name))
        
        u <- switch(dist_name,
                    beta = pbeta(x, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]),
                    sn   = psn(x, xi = coef(fit)[1], omega = coef(fit)[2], alpha = coef(fit)[3]),
                    norm = pnorm(x, mean = fit$estimate["mean"], sd = fit$estimate["sd"]))
        
        qnorm(pmin(pmax(u, 1e-10), 1 - 1e-10))
      }
      
      data_numeric <- data[, names(marginals)]
      data_copula <- as.data.frame(mapply(function(x, dist) gaussian_copula_transform(x, dist),
                                          data_numeric, marginals, SIMPLIFY = FALSE))
      
      # GRM decomposition
      grm_sparse <- Matrix(grm, sparse = FALSE)
      res <- eigen(grm_sparse)
      A <- diag(res$values)
      Q <- as.matrix(res$vectors)
      
      # Projection
      P <- t(Q) %*% as.matrix(data_copula)
      P <- as.data.frame(P)
      colnames(P) <- colnames(data_numeric)
      P$sex <- sex
      
      # Optional back-transformation
      if (backtransform) {
        inverse_transform <- function(z, dist_name, fit) {
          u <- pnorm(z)
          switch(dist_name,
                 norm = qnorm(u, mean = fit$estimate["mean"], sd = fit$estimate["sd"]),
                 beta = qbeta(u, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]),
                 sn   = qsn(u, xi = coef(fit)[1], omega = coef(fit)[2], alpha = coef(fit)[3], solver = "RFB"),
                 stop("Unsupported distribution for inverse transformation")
          )
        }
        
        P[names(data_numeric)] <- as.data.frame(
          Map(function(z, dist_name, varname) {
            fit <- switch(dist_name,
                          beta = fitdist(data_numeric[[varname]], "beta", start = list(shape1 = 1, shape2 = 1)),
                          sn   = selm(data_numeric[[varname]] ~ 1, family = "SN"),
                          norm = fitdist(data_numeric[[varname]], "norm"),
                          stop("Unsupported distribution")
            )
            inverse_transform(z, dist_name, fit)
          }, P[names(data_numeric)], marginals, names(data_numeric))
        )
      }

      # Identity matrix for E
      E <- diag(nrow(A))
      
      # Multivariate twin structure
      output <- list()
      if (n_resp > 1) {
        Z_struc <- mglm4twin:::mt_struc(n_resp = n_resp)
        ind_A <- lapply(Z_struc, function(x) kronecker(x, A))
        ind_E <- lapply(Z_struc, function(x) kronecker(x, E))
      }
      
      # Model selection
      if (n_resp == 1) {
        output$matrices <- if (model == "AE") list(ind_E, ind_A) else list(ind_E)
        output$data <- P
      } else if (n_resp > 1 && model == "AE") {
        output$matrices <- c(ind_E, ind_A)
        output$data <- P
      }
      
      # Formula-based transformation
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
mt_copula_3 <- function(n, grm, n_resp, model, formula = NULL, data = NULL, marginals = NULL, backtransform = FALSE) {
  library(fitdistrplus)
  library(copula)
  library(Matrix)
  library(sn)
  
  # Extract numeric data and sex variable
  data_numeric <- data[, sapply(data, is.numeric)]
  sex <- data$sex
  
  # Fit marginal distributions if not provided
  if (is.null(marginals)) {
    marginals <- lapply(data_numeric, function(x) {
      if (all(x > 0 & x < 1)) "beta"
      else if (all(x >= 0 & floor(x) == x)) "pois"
      else if (abs(skewness(x)) > 1) "sn"
      else "norm"
    })
  }
  
  
  adjust_for_beta <- function(x) {
    x[x <= 0] <- 1e-6
    x[x >= 1] <- 1 - 1e-6
    return(x)
  }
  
  mom_beta_start <- function(x) {
    m <- mean(x)
    v <- var(x)
    tmp <- m * (1 - m) / v - 1
    list(shape1 = m * tmp, shape2 = (1 - m) * tmp)
  }
  
  safe_fit_beta <- function(x) {
    x <- adjust_for_beta(x)
    start_vals <- try(mom_beta_start(x), silent = TRUE)
    tryCatch({
      fitdist(x, "beta", start = start_vals)
    }, error = function(e) {
      message("Warning: beta fit failed, switching to normal.")
      fitdist(x, "norm")
    })
  }
  
  # Fit marginal model once
  fit_marginal <- function(x, dist) {
    switch(dist,
           beta = safe_fit_beta(x),
           sn   = selm(x ~ 1, family = "SN"),
           norm = fitdist(x, "norm"),
           stop("Unsupported distribution: ", dist))
  }
  
  cdf_marginal <- function(x, fit, dist) {
    switch(dist,
           beta = pbeta(x, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]),
           sn   = psn(x, xi = coef(fit)[1], omega = coef(fit)[2], alpha = coef(fit)[3]),
           norm = pnorm(x, mean = fit$estimate["mean"], sd = fit$estimate["sd"]),
           stop("Unsupported distribution: ", dist))
  }
  
  data_numeric <- data[, names(marginals)]
  fits <- Map(fit_marginal, data_numeric, marginals)
  data_copula <- as.data.frame(Map(function(x, fit, dist) {
    qnorm(pmin(pmax(cdf_marginal(x, fit, dist), 1e-10), 1 - 1e-10))
  }, data_numeric, fits, marginals))
  
  # GRM decomposition and projection
  grm_eig <- eigen(as.matrix(grm), symmetric = TRUE)
  Q <- grm_eig$vectors
  P <- as.data.frame(t(Q) %*% as.matrix(data_copula))
  colnames(P) <- colnames(data_numeric)
  P$sex <- sex
  
  # Optional back-transformation
  if (backtransform) {
    inverse_transform <- function(z, fit, dist) {
      u <- pnorm(z)
      switch(dist,
             beta = qbeta(u, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]),
             sn   = qsn(u, xi = coef(fit)[1], omega = coef(fit)[2], alpha = coef(fit)[3]),
             norm = qnorm(u, mean = fit$estimate["mean"], sd = fit$estimate["sd"]),
             stop("Unsupported distribution: ", dist))
    }
    
    P[names(data_numeric)] <- as.data.frame(Map(function(z, fit, dist) {
      inverse_transform(z, fit, dist)
    }, P[names(data_numeric)], fits, marginals))
  }
  
  # Matrix structure
  E <- diag(nrow(grm))
  output <- list()
  
  if (n_resp > 1) {
    Z_struc <- mglm4twin:::mt_struc(n_resp = n_resp)
    ind_A <- lapply(Z_struc, function(x) kronecker(x, diag(grm_eig$values)))
    ind_E <- lapply(Z_struc, function(x) kronecker(x, E))
  }
  
  # Select matrices
  if (n_resp == 1) {
    output$matrices <- if (model == "AE") list(ind_E, ind_A) else list(ind_E)
    output$data <- P
  } else {
    output$matrices <- if (model == "AE") c(ind_E, ind_A) else ind_E
    output$data <- P
  }
  
  # Formula-based transformation (optional)
  if (!is.null(formula)) {
    if (length(output) != length(formula)) stop("Formulas don't match matrix components")
    X_list <- lapply(formula, model.matrix, data = data)
    output <- do.call(c, lapply(seq_along(output), function(i) {
      lapply(seq_len(ncol(X_list[[i]])), function(j) {
        X_list[[i]][, j] * output[[i]]
      })
    }))
  }
  
  return(output)
}

marginals <- list(
  Y1 = "beta",
  Y2 = "beta",
  Y3 = "beta",
  age_std = "norm"  # skew-normal
)


mat <- mt_copula_3(
  n = nrow(data_select),
  grm = GRM,
  n_resp = 3,
  model = "AE",  # or "E", depending on what you want
  formula = NULL,
  data = data_select,
  marginals = marginals,
  backtransform = T  # or FALSE if you want to stay in the copula-transformed space
)


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

hist(mat$data$age_std)

form_Y1 <- c(Y1 ~ sex + age_std)
form_Y2 <- c(Y2 ~ sex + age_std)
form_Y3 <- c(Y3 ~ sex + age_std)

form_Y1 <- c(Y1 ~ sex )
form_Y2 <- c(Y2 ~ sex )
form_Y3 <- c(Y3 ~ sex )

link = rep("logit", 3)
variance = rep("binomialP", 3)


res <- mglm4twin(linear_pred = c(form_Y1, form_Y2, form_Y3),
                 matrix_pred = c(mat$matrices),
                 link = link,
                 variance = variance, 
                 data = mat$data)

sum <-summary(res, model = "AE", biometric = T)
sum

initals <- sum$Dispersion$Estimates

control_initial <- mt_initial_values(linear_pred = c(form_Y1, form_Y2, form_Y3), matrix_pred = c(mat$matrices), link = link, variance = variance,data = data_select, Ntrial = NULL)
control_initial$tau <- c(initals)


mat.dense <- mt_grm_1(n=n, grm = GRM, n_resp = 3, model = "AE", data = NULL)


res2 <- mglm4twin(linear_pred = c(form_Y1, form_Y2, form_Y3),
                 matrix_pred = c(mat.dense),
                 link = link, 
                 variance = variance, 
                 control_initial = control_initial,
                 data = data_select)

sum2 <-summary(res2, model = "AE", biometric = T)
sum2



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

### Simulation Pipeline ##
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
n_vals <- c(seq(500, 1000, by = 50), seq(1500, 3000, by = 500))

summaries_list <- list()

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
  #Age <- rbeta(n, shape1 = 0.3*2, shape2 = 0.7*2)*20 + 70 ## from Wagner
  Age <- rnorm(n, mean = 70, sd = 10)
  # Optional: truncate values to keep them within a reasonable range (e.g., 50 to 90)
  Age <- pmin(pmax(Age, 50), 90)
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
      matrix_pred = c(mat),
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
  A_lt <- c(A[1], A[4], A[5], A[2], A[6], A[3])
  E_lt <- c(E[1], E[4], E[5], E[2], E[6], E[3])
  
  # Heritability 
  h2_vals <- A_lt / (A_lt + E_lt)
  
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
save(results, summaries_list, file = "vignettes/simulation_results.Rdata")

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

p1

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
p2

ggsave("vignettes/figures/frob_norm_vs_n.png", plot = p2, width = 6, height = 4, dpi = 600)


#### Copula model ####


# Store results
results.copula <- data.frame(
  n = numeric(),
  runtime_sec = numeric(),
  frob_norm = numeric()
)
summaries_list.copula <- list()

# Set seeds for reproducibility
set.seed(123)

# Sample sizes to loop through
n_vals <- c(seq(500, 1000, by = 50), seq(1500, 3000, by = 500))
#n_vals <- c(seq(500, 1000, by = 50))


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
  #Age <- rbeta(n, shape1 = 0.3*2, shape2 = 0.7*2)*20 + 70 ## from Wagner
  Age <- rnorm(n, mean = 70, sd = 10)
  # Optional: truncate values to keep them within a reasonable range (e.g., 50 to 90)
  Age <- pmin(pmax(Age, 50), 90)
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
  
  marginals <- list(
    Y1 = "beta",
    Y2 = "beta",
    Y3 = "beta",
    age_std = "norm"  # skew-normal
  )
  
  
  mat <- mt_copula_3(
    n = nrow(data),
    grm = GRM,
    n_resp = 3,
    model = "AE",  # or "E", depending on what you want
    formula = NULL,
    data = data,
    marginals = marginals,
    backtransform = T  # or FALSE if you want to stay in the copula-transformed space
  )
  
  
  
  runtime <- system.time({
    res <- mglm4twin(
      linear_pred = list(Y1 ~ sex + age_std, Y2 ~ sex + age_std, Y3 ~ sex + age_std),
      matrix_pred = c(mat$matrices),
      link = rep("logit", 3),
      variance = rep("binomialP", 3),
      data = mat$data
    )
  })["elapsed"]
  
  sum <- summary(res, model = "AE", biometric = TRUE)
  
  ## 6. Reconstruct estimated h2 matrix
  h2_estimate <- matrix(0, 3, 3)
  diag(h2_estimate) <- sum$A_main$Estimates
  h2_estimate[lower.tri(h2_estimate)] <- sum$A_cross$Estimates
  h2_estimate <- h2_estimate + t(h2_estimate) - diag(diag(h2_estimate))
  
  ## 7. Ground truth h2 matrix
  A_lt <- c(A[1], A[4], A[5], A[2], A[6], A[3])
  E_lt <- c(E[1], E[4], E[5], E[2], E[6], E[3])
  
  
  # Heritability 
  h2_vals <- A_lt / (A_lt + E_lt)
  
  h2_matrix <- matrix(0, 3, 3)
  h2_matrix[lower.tri(h2_matrix, diag = TRUE)] <- h2_vals
  h2_matrix <- h2_matrix + t(h2_matrix) - diag(diag(h2_matrix))
  
  ## 8. Frobenius norm
  frob <- norm(h2_matrix - h2_estimate, type = "F")
  
  ## 9. Store results
  results.copula <- rbind(results.copula, data.frame(n = n, runtime_sec = runtime, frob_norm = frob))
  summaries_list.copula[[paste0("n_", n)]] <- sum
}
# View results
print(results.copula)

# Save results to .Rdata
#save(results.copula, summaries_list.copula, file = "vignettes/simulation_results_copula.Rdata")

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
p1 <- ggplot(results.copula, aes(x = n, y = runtime_sec)) +
  geom_line(color = "#1f77b4", linewidth = 1) +
  geom_point(shape = 21, fill = "#1f77b4", size = 2) +
  labs(
    title = "Computation Time vs Sample Size (n)",
    subtitle = "Phenotypes = 3, Type = Multibound",
    x = "Sample Size (n)",
    y = "Computation Time (seconds)"
  ) +
  theme_pub
p1
 

# Plot 2: Frobenius norm
p2 <- ggplot(results.copula, aes(x = n, y = frob_norm)) +
  geom_line(color = "#d62728", linewidth = 1) +
  geom_point(shape = 21, fill = "#d62728", size = 2) +
  labs(
    title = "Estimation Error vs Sample Size (n)",
    subtitle = expression("Frobenius Norm of (" * h^2 * " - " * hat(h)^2 * "), Phenotypes = 3, Type = Multibound"),
    x = "Sample Size (n)",
    y = "Frobenius Norm"
  ) +
  theme_pub
p2


# Combine results

# Combine unique rows by 'n' from results.copula
results.copula_unique <- results.copula %>%
  distinct(n, .keep_all = TRUE)

# Join on 'n' column
combined <- results %>%
  rename(runtime_sec_original = runtime_sec,
         frob_norm_original = frob_norm) %>%
  inner_join(
    results.copula_unique %>%
      rename(runtime_sec_copula = runtime_sec,
             frob_norm_copula = frob_norm),
    by = "n"
  )

# View combined data
print(combined)

# Plot 1: Computation time
p1 <- ggplot(combined, aes(x = n)) +
  geom_line(aes(y = runtime_sec_original, color = "Original"), linewidth = 1) +
  geom_line(aes(y = runtime_sec_copula, color = "Copula-EVD Projection"), linewidth = 1) +
  geom_point(aes(y = runtime_sec_original, fill = "Original"), shape = 21, size = 2) +
  geom_point(aes(y = runtime_sec_copula, fill = "Copula-EVD Projection"), shape = 21, size = 2) +
  labs(
    title = "Computation Time vs Sample Size (n)",
    subtitle = "Phenotypes = 3, Type = Multibound",
    x = "Sample Size (n)",
    y = "Computation Time (seconds)",
    color = "Method",
    fill = "Method"
  ) +
  theme_pub

# Plot 2: Frobenius norm
p2 <- ggplot(combined, aes(x = n)) +
  geom_line(aes(y = frob_norm_original, color = "Original"), linewidth = 1) +
  geom_line(aes(y = frob_norm_copula, color = "Copula-EVD Projection"), linewidth = 1) +
  geom_point(aes(y = frob_norm_original, fill = "Original"), shape = 21, size = 2) +
  geom_point(aes(y = frob_norm_copula, fill = "Copula-EVD Projection"), shape = 21, size = 2) +
  labs(
    title = "Estimation Error vs Sample Size (n)",
    subtitle = expression("Frobenius Norm of (" * h^2 * " - " * hat(h)^2 * "), Phenotypes = 3, Type = Multibound"),
    x = "Sample Size (n)",
    y = "Frobenius Norm",
    color = "Method",
    fill = "Method"
  ) +
  theme_pub

p1
p2

p1 <- ggplot(combined, aes(x = n)) +
  geom_line(aes(y = runtime_sec_original, color = "Original"), linewidth = 1) +
  geom_line(aes(y = runtime_sec_copula, color = "Copula-EVD Projection"), linewidth = 1) +
  geom_point(aes(y = runtime_sec_original, fill = "Original"), shape = 21, size = 2) +
  geom_point(aes(y = runtime_sec_copula, fill = "Copula-EVD Projection"), shape = 21, size = 2) +
  scale_y_log10() +
  labs(
    title = "Computation Time vs Sample Size (n)",
    subtitle = "Log Scale — Phenotypes = 3, Type = Multibound",
    x = "Sample Size (n)",
    y = "Computation Time (seconds, log scale)",
    color = "Method",
    fill = "Method"
  ) +
  theme_pub

p1

### Publication ready #####

library(cowplot)  # for plot_grid

# Base Theme
theme_pub <- theme_minimal(base_size = 16) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

# Define colors manually
method_colors <- c("GRM" = "#1b9e77", "Copula-EVD Projection" = "#d95f02")

# Frobenius Norm Plot
p_frob <- ggplot(combined, aes(x = n)) +
  geom_line(aes(y = frob_norm_original, color = "GRM"), linewidth = 1) +
  geom_line(aes(y = frob_norm_copula, color = "Copula-EVD Projection"), linewidth = 1) +
  geom_point(aes(y = frob_norm_original, fill = "GRM"), shape = 21, size = 2, color = "black") +
  geom_point(aes(y = frob_norm_copula, fill = "Copula-EVD Projection"), shape = 21, size = 2, color = "black") +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Estimation Error vs Sample Size",
    subtitle = expression("Frobenius norm of (" * h^2 * " - " * hat(h)^2 * "), Traits = 3, Type = Multibound"),
    x = "Sample Size (n)",
    y = "Frobenius Norm",
    color = "Method",
    fill = "Method"
  ) +
  theme_pub

# Runtime Plot (log scale)
p_runtime <- ggplot(combined, aes(x = n)) +
  geom_line(aes(y = runtime_sec_original, color = "GRM"), linewidth = 1) +
  geom_line(aes(y = runtime_sec_copula, color = "Copula-EVD Projection"), linewidth = 1) +
  geom_point(aes(y = runtime_sec_original, fill = "GRM"), shape = 21, size = 2, color = "black") +
  geom_point(aes(y = runtime_sec_copula, fill = "Copula-EVD Projection"), shape = 21, size = 2, color = "black") +
  scale_y_log10() +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Computation Time vs Sample Size",
    subtitle = "Traits = 3, Type = Multibound",
    x = "Sample Size (n)",
    y = "Time (seconds, log scale)",
    color = "Method",
    fill = "Method"
  ) +
  theme_pub

# Combine Plots
final_plot <- plot_grid(p_frob, p_runtime, labels = c("A", "B"), ncol = 2, align = "v")
print(final_plot)

ggsave("vignettes/figures/method_comparison.pdf", final_plot, width = 12, height = 6, dpi = 300)

#### Real sparrow data ####
### Heritability Loop ###


library(genio)
library(AGHmatrix)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

# Phenotype and model formula (reuse your existing objects)
phenotypes <- read.table("/Users/gordonplri/Documents/Genomic McGLM/GMDS/doi_10_5061_dryad_hp758sn__v20180716/LundreganEtAl_PhenosAge1.txt", 
                         header = TRUE, sep = "\ ", stringsAsFactors = FALSE)
data.sparrow <- phenotypes %>%
  select(age1billD, age1billL, age1mass, sex, island, hatchyear)


mt_copula_sparrow <- function(n, grm, n_resp, n_covariate = 0, model, formula = NULL, data = NULL, marginals = NULL, backtransform = FALSE) {
  library(fitdistrplus)
  library(copula)
  library(Matrix)
  library(sn)
  
  # Ensure factor encoding for categorical variables
  data[] <- lapply(data, function(col) {
    if (is.character(col) || is.factor(col)) factor(col) else col
  })
  
  # Identify variable types
  response_vars <- names(data)[1:n_resp]
  covariate_vars <- if (n_covariate > 0) names(data)[(n_resp + 1):(n_resp + n_covariate)] else character(0)
  cat_vars <- setdiff(names(data), c(response_vars, covariate_vars))
  
  # Combine response + covariate for transformation
  data_numeric <- data[, c(response_vars, covariate_vars), drop = FALSE]
  
  # Set default marginals if not supplied
  if (is.null(marginals)) {
    marginals <- rep("norm", ncol(data_numeric))
  }
  
  # Helper functions
  fit_marginal <- function(x, dist) {
    switch(dist,
           beta = safe_fit_beta(x),
           sn   = selm(x ~ 1, family = "SN"),
           norm = fitdist(x, "norm"),
           stop("Unsupported distribution: ", dist))
  }
  
  cdf_marginal <- function(x, fit, dist) {
    switch(dist,
           beta = pbeta(x, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]),
           sn   = psn(x, xi = coef(fit)[1], omega = coef(fit)[2], alpha = coef(fit)[3]),
           norm = pnorm(x, mean = fit$estimate["mean"], sd = fit$estimate["sd"]),
           stop("Unsupported distribution: ", dist))
  }
  
  # Fit marginals and transform via copula (probit)
  fits <- Map(fit_marginal, data_numeric, marginals)
  data_copula <- as.data.frame(Map(function(x, fit, dist) {
    qnorm(pmin(pmax(cdf_marginal(x, fit, dist), 1e-10), 1 - 1e-10))
  }, data_numeric, fits, marginals))
  
  colnames(data_copula) <- colnames(data_numeric)
  
  # Combine transformed and categorical data
  data_copula <- cbind(data_copula, data[, cat_vars, drop = FALSE])
  
  # GRM eigendecomposition
  grm_eig <- eigen(as.matrix(grm), symmetric = TRUE)
  Q <- grm_eig$vectors
  
  # Project response variables only
  P_resp <- as.data.frame(t(Q) %*% as.matrix(data_copula[, response_vars]))
  colnames(P_resp) <- response_vars
  
  # Keep covariates (transformed and categorical) as-is
  covariates_all <- data_copula[, setdiff(names(data_copula), response_vars), drop = FALSE]
  P <- cbind(P_resp, covariates_all)
  
  # Optional backtransformation for response vars
  if (backtransform) {
    inverse_transform <- function(z, fit, dist) {
      u <- pnorm(z)
      switch(dist,
             beta = qbeta(u, shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]),
             sn   = qsn(u, xi = coef(fit)[1], omega = coef(fit)[2], alpha = coef(fit)[3]),
             norm = qnorm(u, mean = fit$estimate["mean"], sd = fit$estimate["sd"]),
             stop("Unsupported distribution: ", dist))
    }
    
    P[response_vars] <- as.data.frame(Map(function(z, fit, dist) {
      inverse_transform(z, fit, dist)
    }, P[response_vars], fits[1:n_resp], marginals[1:n_resp]))
  }
  
  # Build variance matrices
  E <- diag(nrow(grm))
  output <- list()
  
  if (n_resp > 1) {
    Z_struc <- mglm4twin:::mt_struc(n_resp = n_resp)
    ind_A <- lapply(Z_struc, function(x) kronecker(x, diag(grm_eig$values)))
    ind_E <- lapply(Z_struc, function(x) kronecker(x, E))
  }
  
  if (n_resp == 1) {
    output$matrices <- if (model == "AE") list(ind_E, ind_A) else list(ind_E)
    output$data <- P
  } else {
    output$matrices <- if (model == "AE") c(ind_E, ind_A) else ind_E
    output$data <- P
  }
  
  # Create fixed effect matrices if formulas are provided
  if (!is.null(formula)) {
    if (length(output$matrices) != length(formula)) stop("Formulas don't match matrix components")
    X_list <- lapply(formula, model.matrix, data = P)
    output$matrices <- do.call(c, lapply(seq_along(output$matrices), function(i) {
      lapply(seq_len(ncol(X_list[[i]])), function(j) {
        X_list[[i]][, j] * output$matrices[[i]]
      })
    }))
  }
  
  return(output)
}

form_billD <- age1billD ~ sex + island + hatchyear + age1mass
form_billL <- age1billL ~ sex + island + hatchyear + age1mass

chromosomes <- paste0("chr", c(1:15, 17:29))  # example; exclude 13 if missing, adjust as needed
#chromosomes <- paste0("chr", c(1:3))  # example; exclude 13 if missing, adjust as needed
base_path <- "/Users/gordonplri/Documents/Genomic McGLM/GMDS/doi_10_5061_dryad_hp758sn__v20180716/plink_by_chr/"

results_list <- list()
runtimes <- data.frame(chromosome = character(), runtime_secs = numeric(), stringsAsFactors = FALSE)

for (chr in chromosomes) {
  cat("Processing", chr, "\n")
  start_time <- Sys.time()
  
  plink_prefix <- file.path(base_path, chr, paste0(chr,"_LD_pruned"))
  # Read PLINK data
  plink_data <- read_plink(plink_prefix)
  geno_mat <- t(plink_data$X)
  geno_mat <- as.matrix(geno_mat)
  rownames(geno_mat) <- plink_data$fam$id
  colnames(geno_mat) <- plink_data$bim$id
  
  # Calculate GRM
  num.SNP <- as.numeric(ncol(Mcheck(geno_mat, missingValue = -9, thresh.maf = 0.05)))
  GRM <- Gmatrix(geno_mat, missingValue = -9, maf = 0.05, method = "VanRaden")
  
  
  mat.output <- mt_copula_sparrow(
    n = nrow(data.sparrow),
    grm = GRM,
    n_resp = 2,
    n_covariate = 1,
    model = "AE",
    formula = NULL,
    data = data.sparrow,
    marginals = rep("norm", 3),  # 2 responses + 1 covariate
    backtransform = TRUE
  )
  
  # ensure factor encoding 
  
  mat.output$data <- mat.output$data %>% 
    mutate(
      sex = factor(sex),
      island = factor(island),
      hatchyear = factor(hatchyear)
    )
  
  
  # Fit model
  res.sparrow <- mglm4twin(
    linear_pred = c(form_billD, form_billL),
    matrix_pred = c(mat.output$matrices),
    link = rep("identity", 2),
    variance = rep("constant", 2),
    data = mat.output$data
  )
  
  sum.sparrow <- summary(res.sparrow, model = "AE", biometric = TRUE)
  
  # Extract and format results for plotting
  df_chr <- bind_rows(
    sum.sparrow$A_main %>%
      mutate(component = rownames(.), type = "main"),
    sum.sparrow$A_cross %>%
      mutate(component = rownames(.), type = "cross")
  )
  
  df_chr$component <- recode(df_chr$component,
                             "h1" = "Bill Depth",
                             "h2" = "Bill Length",
                             "h12" = "Depth × Length")
  df_chr$component <- factor(df_chr$component, levels = c("Bill Depth", "Bill Length", "Depth × Length"))
  df_chr$chromosome <- chr
  df_chr$SNP <- num.SNP
  
  results_list[[chr]] <- df_chr
  
  end_time <- Sys.time()
  runtimes <- rbind(runtimes, data.frame(chromosome = chr, runtime_secs = as.numeric(difftime(end_time, start_time, units = "secs"))))
}


#### Total h^2 heritability ####
# Libraries
library(dplyr)
library(genio)   # for reading PLINK data

# Load phenotype data (sparrow morphology + covariates)
phenotype_raw <- read.table(
  "/Users/gordonplri/Documents/Genomic McGLM/GMDS/doi_10_5061_dryad_hp758sn__v20180716/LundreganEtAl_PhenosAge1.txt",
  header = TRUE, sep = "\ ", stringsAsFactors = FALSE
)

# Select relevant phenotype columns
phenotype_sparrow <- phenotype_raw %>%
  select(age1billD, age1billL, age1mass, sex, island, hatchyear)

# Define model formulas for each trait
formula_bill_depth <- age1billD ~ sex + island + hatchyear + age1mass
formula_bill_length <- age1billL ~ sex + island + hatchyear + age1mass

# Read PLINK data (excluding chromosome 16)
plink_path <- "/Users/gordonplri/Documents/Genomic McGLM/GMDS/doi_10_5061_dryad_hp758sn__v20180716/plink_total/chr_all_no16"
plink_data <- genio::read_plink(plink_path)

# Transpose genotype matrix to [SNPs x Individuals]
genotypes <- t(plink_data$X)
genotypes <- as.matrix(genotypes)
rownames(genotypes) <- plink_data$fam$id
colnames(genotypes) <- plink_data$bim$id

# Compute GRM using VanRaden method
GRM_sparrow <- Gmatrix(
  genotypes,
  missingValue = -9,
  maf = 0.05,
  method = "VanRaden"
)

# Fit multivariate copula AE model
copula_fit <- mt_copula_sparrow(
  n = nrow(phenotype_sparrow),
  grm = GRM_sparrow,
  n_resp = 2,
  n_covariate = 1,
  model = "AE",
  formula = NULL,
  data = phenotype_sparrow,
  marginals = rep("norm", 3),  # 2 responses + 1 covariate
  backtransform = TRUE
)

# Recode covariates as factors
copula_fit$data <- copula_fit$data %>%
  mutate(
    sex = factor(sex),
    island = factor(island),
    hatchyear = factor(hatchyear)
  )

# Fit GREML model using mglm4twin
greml_fit <- mglm4twin(
  linear_pred = c(formula_bill_depth, formula_bill_length),
  matrix_pred = c(copula_fit$matrices),
  link = rep("identity", 2),
  variance = rep("constant", 2),
  data = copula_fit$data
)

# Summarize results
greml_summary <- summary(greml_fit, model = "AE", biometric = TRUE)

# Combine results from main and cross components
heritability_df <- bind_rows(
  greml_summary$A_main %>% mutate(component = rownames(.), type = "main"),
  greml_summary$A_cross %>% mutate(component = rownames(.), type = "cross")
)

# Label components for clarity
heritability_df$component <- recode(heritability_df$component,
                                    "h1" = "Bill Depth",
                                    "h2" = "Bill Length",
                                    "h12" = "Depth × Length")
heritability_df$component <- factor(heritability_df$component,
                                    levels = c("Bill Depth", "Bill Length", "Depth × Length"))

# Add chromosome label (for consistency in downstream plots)
heritability_df$chromosome <- "all_no16"

# Output
print(heritability_df)


### Visualize results
library(ggplot2)
library(dplyr)
library(ggrepel)

# Define colors for components
cbf_colors <- c(
  "Bill Depth" = "#0072B2",
  "Bill Length" = "#E69F00",
  "Depth × Length" = "#D55E00"
)

# Prepare data: compute CI and SNP in Mbp
results_df <- results_df %>%
  mutate(
    SNP_Mbp = SNP / 1e6,
    CI_lower = Estimates - 1.96 * std.error,
    CI_upper = Estimates + 1.96 * std.error
  )

# Get slopes (heritability trend) for annotation
slopes <- results_df %>%
  group_by(component, type) %>%
  summarize(
    slope = coef(lm(Estimates ~ SNP_Mbp))[2],
    .groups = "drop"
  )

# Identify chromosomes with highest heritability per component and type
max_h2 <- results_df %>%
  group_by(component) %>%
  arrange(desc(Estimates)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(label = chromosome) %>%
  mutate(label = as.character(readr::parse_number(chromosome)))


annotation_df <- tibble::tibble(
  component = factor(
    c("Bill Depth", "Bill Length", "Depth × Length"),
    levels = c("Bill Depth", "Bill Length", "Depth × Length")
  ),
  x = 0.0005,    # x position in SNP_Mbp (adjust as needed)
  y = 0.44,      # y position for all labels (adjust if needed)
  label = c(
    "h² = 0.46 ± 0.06 (95% CI)",
    "h² = 0.50 ± 0.07 (95% CI)",
    "h²_biv = 0.47 ± 0.18 (95% CI)"  )
)

# Plot
scatter_heritability <- ggplot(results_df, aes(x = SNP_Mbp, y = Estimates,
                       ymin = CI_lower, ymax = CI_upper,
                       color = component)) +
  
  # CI range
  geom_pointrange(position = position_dodge(width = 0.5), 
                  fatten = 1.2, linewidth = 0.7, size = 1.5) +
  
  # LM fit with color-matched ribbon
  geom_smooth(aes(fill = component, group = interaction(component, "dashed"), linetype = "dashed"),
              method = "lm", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  
  # Horizontal line at h2 = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  
  # Highlight top heritability chromosomes
  
  # Add number-in-circle points and text overlays
  geom_point(data = max_h2,
             aes(x = SNP / 1e6, y = Estimates),
             size = 5.5,
             shape = 21,            # Filled circle
             fill = "white",
             color = "black",
             stroke = 1.2,
             position = position_nudge(y = 0)) +
  
  geom_text(data = max_h2,
            aes(x = SNP / 1e6, y = Estimates, label = label),
            fontface = "bold",
            size = 4,
            color = "black",
            position = position_nudge(y = 0)) +
  
  # Facet by component
  facet_wrap(~ component) +
  
  # Scales and labels
  scale_color_manual(values = cbf_colors) +
  scale_fill_manual(values = cbf_colors) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(-0.0009, 0.009)
  ) + scale_y_continuous(
    limits = c(min(results_df$CI_lower-0.003, na.rm = TRUE), 0.45),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Number of SNPs (Mbp)",
    y = expression(h[c]^2 ~ "(± 95% CI)"),
    color = "Component",
    fill = "Component",
    linetype = "Effect Type"
  ) +
  # Add annotation for overall h² estimates (multivariate GREML excluding chr16)
   geom_label(data = annotation_df,
               aes(x = x, y = y, label = label),
               inherit.aes = FALSE,
               size = 6,
               hjust = 0,
               vjust = 1,
               fill = "white",
               color = "black",
               label.size = 0) +
  # Theme
  theme_bw(base_size = 16) +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


ggsave(
  filename = "~/Documents/Genomic McGLM/GMDS/figures/figure_sparrow_heritability.pdf",
  plot = scatter_heritability,
  device = "pdf",
  width = 17.5,
  height = 7.5,
  units = "in"
)




# Confidence Intervalls:

# Define chromosome order without "chr" prefix
chrom_order <- setdiff(1:29, 16)

# Prepare data
df_all <- df_all %>%
  mutate(
    chrom_num = as.numeric(factor(chromosome, levels = paste0("chr", chrom_order))),
    chrom_label = chrom_order[chrom_num]
  )



# Color-blind friendly colors (Okabe-Ito)

df_all$component <- factor(df_all$component, levels = c("Bill Depth", "Bill Length", "Depth × Length"))

cbf_colors <- c(
  "Bill Depth" = "#0072B2",    # blue
  "Bill Length" = "#E69F00",   # orange
  "Depth × Length" = "#D55E00" # vermillion/red
)


# Compute y-axis limits for 95% CI
y_min <- floor(min(df_all$Estimates - 1.96 * df_all$std.error) * 10) / 10
y_max <- ceiling(max(df_all$Estimates + 1.96 * df_all$std.error) * 10) / 10

# Plot
p <- ggplot(df_all, aes(x = chrom_num, y = Estimates, group = component)) +
  geom_pointrange(
    aes(
      ymin = Estimates - 1.96 * std.error,
      ymax = Estimates + 1.96 * std.error,
      color = component,
      shape = component
    )
  ) +
  geom_point(aes(color = component, shape = component),
             size = 3.5) +
  facet_wrap(~component, scales = "fixed", ncol = 1) +
  scale_x_continuous(name = "Chromosome",
                     breaks = df_all$chrom_num %>% unique(),
                     labels = chrom_order) +
  #scale_color_jco(name = NULL) +
  scale_shape_manual(name = NULL, values = c(16, 17, 15)) +
  scale_y_continuous(name = expression(h^2~Estimate~"(95% CI)"),
                     limits = c(y_min, y_max)) +
  theme_pubclean(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 18),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.position = "none",
    #panel.spacing = unit(0, "lines"),
    axis.title.x = element_text(size = 16)
  ) + 
  # Remove all horizontal grid lines
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  # Add solid horizontal line at y=0 (x-axis)
  geom_hline(yintercept = 0, color = "black", size = 0.7, alpha =0.7) +
  scale_color_manual(values = cbf_colors) +
  scale_fill_manual(values = cbf_colors)

p

ggsave("~/Documents/Genomic McGLM/GMDS/figures/sparrow_heritability_plot.pdf", plot = p, width = 12, height = 10, units = "in", device = cairo_pdf)

# View runtime summary
print(runtimes)


#### PCA on GRM ####
# Load necessary libraries
library(FactoMineR)
library(factoextra)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(ggplot2)

# Count individuals per island
island_counts <- phenotype_raw %>%
  count(island, name = "n_individuals")

# PCA
res.pca <- PCA(GRM, graph = FALSE)

# Extract PCA scores
pca_coords <- as.data.frame(res.pca$ind$coord[, 1:2])
pca_coords$ID <- rownames(pca_coords)

# Add island group (ensure ID formats match)
group_info <- phenotype_raw %>%
  mutate(id = as.character(id)) %>%
  select(id, island)

# Join island info to PCA coordinates
pca_coords <- pca_coords %>%
  mutate(id = rownames(.)) %>%
  left_join(group_info, by = "id") %>%
  mutate(island = as.factor(island))

# Add island name mapping
island_map <- tibble::tibble(
  island = factor(c(20, 22, 23, 24, 26, 27, 28, 38)),
  island_name = c("Nesøya", "Myken", "Selvær", "Træna",
                  "Gjerøy", "Hestmannøy", "Indre Kvarøy", "Aldra")
)

# Merge readable island names
pca_coords <- pca_coords %>%
  left_join(island_map, by = "island")

# Variance explained
expl_var <- res.pca$eig[1:2, 2]


# Okabe-Ito colorblind-friendly palette (8 colors)
cbf_palette <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7",  # reddish purple
  "#999999"   # gray
)

pca_islands <- ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = island_name)) +
  stat_ellipse(level = 0.95, size = 1) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = cbf_palette) +
  labs(
    title = "PCA of Genomic Relationship Matrix",
    x = paste0("PC1 (", round(expl_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(expl_var[2], 1), "% variance)"),
    color = "Island"
  ) +
  theme_apa(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.size = unit(0.6, "cm"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16)
  )

print(pca_islands)


ggsave("~/Documents/Genomic McGLM/GMDS/figures/pca_islands.pdf", plot = pca_islands, width = 10, height = 8)



### Correlation plot Phenotypes ###

# Load required libraries
# Load libraries
library(GGally)
library(ggplot2)
library(dplyr)
library(ggthemes)

# Define colorblind-friendly palette with island names
island_palette <- c(
  "Aldra" = "#E69F00",
  "Gjerøy" = "#56B4E9",
  "Hestmannøy" = "#009E73",
  "Indre Kvarøy" = "#F0E442",
  "Myken" = "#0072B2",
  "Nesøya" = "#D55E00",
  "Selvær" = "#CC79A7",
  "Træna" = "#999999"
)

# Mapping from numeric island codes to names (example: adjust as needed)
island_lookup <- c(
  "20" = "Nesøya",
  "22" = "Myken",
  "23" = "Træna",
  "24" = "Selvær",
  "26" = "Gjerøy",
  "27" = "Hestmannøy",
  "28" = "Indre Kvarøy",
  "38" = "Aldra"
)

# Prepare phenotype data and rename variables
pheno_data <- phenotype_raw %>%
  mutate(island_name = factor(island_lookup[as.character(island)], levels = names(island_palette))) %>%
  select(island_name, age1billD, age1billL) %>%
  rename(
    `Bill Depth` = age1billD,
    `Bill Length` = age1billL,
  )

# Create the ggpairs plot
pheno_plot <- ggpairs(
  data = pheno_data,
  columns = 2:3,  # phenotype columns
  mapping = aes(color = island_name),
  upper = list(continuous = wrap("cor", size = 6)),
  lower = list(continuous = wrap("points", alpha = 0.8, size = 1.5)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5))
) +
  scale_color_manual(values = okabe_ito, name = "Island") +
  theme_apa(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 16),
    strip.text = element_text(face = "bold")
  )

# Show the plot
print(pheno_plot)


# Export the correlation plot to PDF
ggsave(
  filename = "~/Documents/Genomic McGLM/GMDS/figures/correlation_plot_island.pdf",
  plot = pheno_plot,
  # Ensures good font rendering
  width = 10,               # Width in inches
  height = 8,              # Height in inches
)














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
