# Compute the coverage proportion using classical theory or standard bootstrap 
library(glmhd)

# 1 -- Classical or high-dimensional theory 
model <- "" # set up model: logistic/probit/poisson
covariate <- "" # set up covariate distribution: mvt/iid/arch/small
coef <- "" # set up model coefficient: sparse/notsparse 
p <- 400  # set up number of variables: 400/40
n <- 4000 # set up number of observations: 4000/40

# -- Set up model family
if(model == "logistic"){
  family <- get("binomial")(link = "logit")
}else if(model == "probit"){
  family <- get("binomial")(link = "probit")
}else if(model == "poisson"){
  family <- get("poisson")(link = "log")
}
# -- Set up model coef
if(model == "logistic"){
  if(covariate == "mvt" | covariate == "arch"){
    if(coef == "sparse"){
      beta <- scan("beta_sparse.txt")
    }else{
      beta <- scan("beta_large.txt")
    }
  }else if(covariate == "small"){
    beta <- scan("beta_small.txt")
  }
}else if(model == "probit" | model == "poisson"){
  beta <- scan("beta_probit.txt")
}

# -- Set up output location
if(model == "logistic"){
  if(coef == "sparse"){
    j <- c(84, 83) # individual coordinates to compute coverage
    outloc <- "/output/sparse/"
    covloc <- "/result/sparse/"
  }else if(covariate == "mvt"){
    outloc <- "/output/logistic/"
    covloc <- "/result/logistic/"
  }else if(covariate == "arch"){
    outloc <- "/output/arch/"
    covloc <- "/result/arch/"
  }else if(covariate == "small"){
    outloc <- "/output/small/"
    covloc <- "/result/small/"
  }
}else if(model == "probit"){
  j <- c(264, 272)
  outloc <- "/output/probit/"
  covloc <- "/result/probit/"
}else if(model == "poisson"){
  j <- c(88, 284)
  outloc <- "/output/poisson/"
  covloc <- "/result/poisson/"
}
# -- Set up covariance matrix
if(covariate == "mvt"){
  rho <- 0.5
  x <- rho^(c(0:(p/2), (p/2-1):1))
  Sigma <- toeplitz(x)
}else if(covariate == "arch" | covariate == "small"){
  Sigma <- matrix(0, p, p)
  diag(Sigma) <- 1
}

gamma <- sqrt(t(beta) %*% (Sigma %*% beta) / p) # signal strength 

# 1.1 -- parameters of HDT
kappa <- p / n
if(model == "logistic"){
  param <- find_param(kappa = kappa, gamma = gamma, intercept = F)
}else if(model == "probit"){
  rho_prime <- function(t) pnorm(t)
  f_prime1 <- function(t) - dnorm(t) / pnorm(t)
  f_prime0 <- function(t) dnorm(t) / pnorm(-t)
  param <- find_param(kappa = kappa, gamma = gamma, intercept = F,
                      rho_prime = rho_prime,
                      f_prime1 = f_prime1,
                      f_prime0 = f_prime0)
}

# load the MLE and std.dev (from glm function)
mle <- read.table(paste0(outloc, "repeated/mle_repeated.txt"))
sd <- read.table(paste0(outloc, "repeated/sd_repeated.txt"))

# coverage proportion
alpha <- c(0.05, 0.1, 0.2)

for(a in alpha){
  cat("alpha = ", a, ": \n")

  tau <- 1 / sqrt(diag(solve(Sigma)))
  sd_high_dim <- param[3] / tau 
  
  lower <- mle - qnorm(1 - a / 2) * sd
  upper <- mle + qnorm(1 - a / 2) * sd
  covered <- beta < t(upper) & beta > t(lower) # classical 
  
  lower_hd <- (t(mle) - qnorm(1 - a / 2) * sd_high_dim) / param[1]
  upper_hd <- (t(mle) + qnorm(1 - a / 2) * sd_high_dim) / param[1]
  covered_hd <- beta < upper_hd & beta > lower_hd # HDT
  
  cat("Classical: \n")
  cat("single: mean =  ", rowMeans(covered)[j] * 100, "; sd = ", apply(covered[j, ], 1, sd) / sqrt(ncol(covered_hd)) * 100, "\n")
  cat("bulk: mean =  ", mean(colMeans(covered)) * 100, "; sd = ", sd(colMeans(covered))/sqrt(ncol(covered_hd)) * 100, "\n")
  cat("High-dim: \n")
  cat("single: mean =  ", rowMeans(covered_hd)[j] * 100, "; sd = ", apply(covered_hd[j, ], 1, sd) / sqrt(ncol(covered_hd)) * 100, "\n")
  cat("bulk: mean =  ", mean(colMeans(covered_hd)) * 100, "; sd = ", sd(colMeans(covered_hd))/sqrt(ncol(covered_hd)) * 100, "\n")
}

# 2 -- Standard bootstrap
# This is computed only for the logistic regression when covariates are from MVT 
model <- "logistic" 
covariate <- "mvt" 
p <- 400 
jnull <- 39
jnonnull <- 233

# scan coef
beta <- scan("/beta_large.txt")

outloc <- "/output/logistic_standard_boot/"
names <- list.files(outloc)
indices <- unique(sapply(names, function(t) strsplit(t, "_")[[1]][2]))
N <- length(indices)


lower_pairs <- matrix(0, p, N)
upper_pairs <-matrix(0, p, N)
lower_param <- matrix(0, p, N)
upper_param <- matrix(0, p, N)
covered_pairs <- NULL
covered_param <- NULL
for(i in 1:N){
  k <-  indices[i]
  
  boot_pairs <-  tryCatch(read.table(paste0(outloc, "PairsBoot_", k)), error = function(e) {return(-1)})
  if(length(boot_pairs) == 1) next;
  
  lower_pairs[,i] <- apply(boot_pairs, 2, function(t) quantile(t, probs = alpha / 2))
  upper_pairs[,i] <- apply(boot_pairs, 2, function(t) quantile(t, probs = 1 - alpha / 2))
  covered_new <- beta < upper_pairs[,i] & beta > lower_pairs[,i]
  covered_pairs <- cbind(covered_pairs, covered_new)
  
  boot_param <-  tryCatch(read.table(paste0(outloc, "ParamBoot_", k)), error = function(e) {return(-1)})
  if(length(boot_param) == 1) next;
  
  lower_param[,i] <- apply(boot_param, 2, function(t) quantile(t, probs = alpha / 2))
  upper_param[,i] <- apply(boot_param, 2, function(t) quantile(t, probs = 1 - alpha / 2))
  covered_new <- beta < upper_param[,i] & beta > lower_param[,i]
  covered_param <- cbind(covered_param, covered_new)
}

covloc <- "/result/logistic/"
if(!dir.exists(covloc)){dir.create(covloc)}
write.table(covered_pairs, file = paste0(covloc, alpha, "_covered_pairs.txt"), row.names = F, col.names = F)
write.table(covered_param, file = paste0(covloc, alpha, "_covered_param.txt"), row.names = F, col.names = F)

alpha <- c(0.05, 0.1)
for(a in alpha){
  cat("alpha = ", a, ": \n")
  # Parametric bootstrap
  covered_param <- read.table(paste0(covloc, a, "_covered_param.txt"))
  # Single MLE covered
  single_covered <- rowMeans(covered_param)
  cat(a, ": single - param (null) = ", single_covered[jnull] * 100, ", sd = ", sqrt(single_covered[jnull] * (1-single_covered[jnull])) / sqrt(ncol(covered_param))  * 100 , "\n")
  cat(a, ": single - param (non-null) = ", single_covered[jnonnull] * 100, ", sd = ", sqrt(single_covered[jnonnull] * (1-single_covered[jnonnull])) / sqrt(ncol(covered_param))  * 100 , "\n")
  
  # Bulk covered
  bulk_covered_param <- colMeans(covered_param)
  cat(a, ": bulk - param = ", mean(bulk_covered_param)  * 100, ", sd = ", sd(bulk_covered_param)/sqrt(length(bulk_covered_param))  * 100, "\n")
  
  # Pairs bootstrap
  covered_pairs <- read.table(paste0(covloc, a, "_covered_pairs.txt"))
  # Single MLE covered
  single_covered <- rowMeans(covered_pairs)
  cat(a, ": single - pairs (null) = ", single_covered[jnull] * 100, ", sd = ", sqrt(single_covered[jnull] * (1-single_covered[jnull])) / sqrt(ncol(covered_pairs))  * 100 , "\n")
  cat(a, ": single - pairs (non-null) = ", single_covered[jnonnull] * 100, ", sd = ", sqrt(single_covered[jnonnull] * (1-single_covered[jnonnull])) / sqrt(ncol(covered_pairs))  * 100 , "\n")
  
  # Bulk covered
  bulk_covered_pairs <- colMeans(covered_pairs)
  cat(a, ": bulk - pairs = ", mean(bulk_covered_pairs)  * 100, ", sd = ", sd(bulk_covered_pairs)/sqrt(length(bulk_covered_pairs))  * 100, "\n")
}

