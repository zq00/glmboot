# Compute the coverage using the resized bootstrap 

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = T)
model <- args[1]  # set up model: logistic/probit/poisson
covariate <- args[2] # set up covariate distribution: mvt/iid/arch/small
coef = args[3] # set up model coefficient: sparse/notsparse 
p <- as.numeric(args[4]) # set up number of variables: 400/40
alpha <- as.numeric(args[5]) # confidence level 

## A -- Set-up
# 1 -- Scan model coef 
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

# 2 -- Set up output location
if(model == "logistic"){
  if(coef == "sparse"){
    outloc <- "/output/sparse/" # output location (data MLE etc.)
    covloc <- "/result/sparse/" # where to store coverage
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
  outloc <- "/output/probit/"
  covloc <- "/result/probit/"
}else if(model == "poisson"){
  outloc <- "/output/poisson/"
  covloc <- "/result/poisson/"
}

# 3 -- Setting the covariance matrix 
if(covariate == "mvt"){
  rho <- 0.5
  x <- rho^(c(0:(p/2), (p/2-1):1))
  Sigma <- toeplitz(x)
}else if(covariate == "arch" | covariate == "small"){
  Sigma <- matrix(0, p, p)
  diag(Sigma) <- 1
}
# Compute gamma  
gamma <- sqrt(t(beta) %*% (Sigma %*% beta) / p)

# List all the files in your output directory (each file corresponds to one repetition)
names <- list.files(outloc)
indices <- unique(sapply(names, function(t) strsplit(t, "_")[[1]][2]))
N <- length(indices)

# B -- Bootstrap t
lower_known <- matrix(0, p, N)
upper_known <-matrix(0, p, N)
lower_unknown <- matrix(0, p, N)
upper_unknown <- matrix(0, p, N)
covered_known <- NULL
covered_unknown <- NULL
for(i in 1:N){
  k <-  indices[i]
  
  beta_hat <- tryCatch(scan(paste0(outloc, "mle_", k)), error = function(e) {return(-1)})
  betas_known <-  tryCatch(scan(paste0(outloc, "betas_known_", k)), error = function(e) {return(-1)})
  boot_known <-  tryCatch(read.table(paste0(outloc, "boot_known_", k)), error = function(e) {return(-1)})
  if(length(beta_hat) == 1 | length(boot_known) == 1) next;
  
  sd_boot <- apply(boot_known, 2, sd)
  alpha_boot <- lm(colMeans(boot_known)~betas_known + 0, weights = 1/sd_boot^2)$coef
  t_known <- (t(boot_known) - alpha_boot * betas_known) / sd_boot # each column represents one experiment
  
  lower_known[,i] <- (beta_hat - sd_boot * apply(t_known, 1, function(t) quantile(t, probs = (1 - alpha / 2)))) / alpha_boot
  upper_known[,i] <- (beta_hat - sd_boot * apply(t_known, 1, function(t) quantile(t, probs = (alpha / 2)))) / alpha_boot
  covered_new <- beta < upper_known[,i] & beta > lower_known[,i]
  covered_known <- cbind(covered_known, covered_new)
  
  # cat(i, ", (known)", mean(covered_new), "\n")
  
  betas_unknown <-  tryCatch(scan(paste0(outloc, "betas_unknown_", k)), error = function(e) {return(-1)})
  boot_unknown <-  tryCatch(read.table(paste0(outloc, "boot_unknown_", k)), error = function(e) {return(-1)})
  if(length(betas_unknown) == 1 | length(boot_unknown) == 1) next;
  boot_unknown <- t(boot_unknown)
  
  sd_boot <- apply(boot_unknown, 2, sd)
  alpha_boot <- lm(colMeans(boot_unknown)~betas_unknown + 0, weights = 1/sd_boot^2)$coef
  t_unknown <- (t(boot_unknown) - alpha_boot * betas_unknown) / sd_boot # each column represents one experiment
  
  lower_unknown[,i] <- (beta_hat - sd_boot * apply(t_unknown, 1, function(t) quantile(t, probs = (1 - alpha / 2)))) / alpha_boot
  upper_unknown[,i] <- (beta_hat - sd_boot * apply(t_unknown, 1, function(t) quantile(t, probs = (alpha / 2)))) / alpha_boot
  
  covered_new <- beta < upper_unknown[,i] & beta > lower_unknown[,i] 
  covered_unknown <- cbind(covered_unknown, covered_new)
  
  cat(i, ", (unknown)", mean(covered_new), "\n")
}

if(!dir.exists(covloc)){dir.create(covloc)}
write.table(covered_known, file = paste0(covloc, alpha, "_covered_known.txt"), row.names = F, col.names = F)
write.table(covered_unknown, file = paste0(covloc, alpha, "_covered_unknown.txt"), row.names = F, col.names = F)

# C -- Bootstrap g
# here we only use the first 100 bootstrap samples to compute the bias and std.
g_lower_known <- matrix(0, p, N)
g_upper_known <-matrix(0, p, N)
g_lower_unknown <- matrix(0, p, N)
g_upper_unknown <- matrix(0, p, N)
g_covered_known <- NULL
g_covered_unknown <- NULL
for(i in 1:N){
  k <-  indices[i]
  
  beta_hat <- tryCatch(scan(paste0(outloc, "mle_", k)), error = function(e) {return(-1)})
  betas_known <-  tryCatch(scan(paste0(outloc, "betas_known_", k)), error = function(e) {return(-1)})
  boot_known <-  tryCatch(read.table(paste0(outloc, "boot_known_", k)), error = function(e) {return(-1)})
  if(length(beta_hat) == 1 | length(boot_known) == 1) next;
  boot_known <- boot_known[1:100, ] # for Gaussian CI, only use the first 100 bootstrap samples
  
  sd_boot <- apply(boot_known, 2, sd)
  alpha_boot <- lm(colMeans(boot_known)~betas_known + 0, weights = 1/sd_boot^2)$coef
  
  g_lower_known[,i] <- (beta_hat - sd_boot * qnorm(1-alpha/2)) / alpha_boot
  g_upper_known[,i] <- (beta_hat - sd_boot * qnorm(alpha/2)) / alpha_boot
  
  covered_new <- beta < g_upper_known[,i] & beta > g_lower_known[,i]
  g_covered_known <- cbind(g_covered_known, covered_new)
  
  cat(i, ", (known - g)", mean(covered_new), "\n")
  
  betas_unknown <-  tryCatch(scan(paste0(outloc, "betas_unknown_", k)), error = function(e) {return(-1)})
  boot_unknown <-  tryCatch(read.table(paste0(outloc, "boot_unknown_", k)), error = function(e) {return(-1)})
  if(length(betas_unknown) == 1 | length(boot_unknown) == 1) next;
  boot_unknown <- t(boot_unknown)
  boot_unknown <- boot_unknown[1:100,]
  
  sd_boot <- apply(boot_unknown, 2, sd)
  alpha_boot <- lm(colMeans(boot_unknown)~betas_unknown + 0, weights = 1/sd_boot^2)$coef
  
  g_lower_unknown[,i] <- (beta_hat - sd_boot * qnorm(1-alpha/2)) / alpha_boot
  g_upper_unknown[,i] <- (beta_hat - sd_boot * qnorm(alpha/2)) / alpha_boot
  
  covered_new <- beta < g_upper_unknown[,i] & beta > g_lower_unknown[,i] 
  g_covered_unknown <- cbind(g_covered_unknown, covered_new)
  
  cat(i, ", (unknown - g)", mean(covered_new), "\n")
}


if(!dir.exists(covloc)){dir.create(covloc)}
# store the coverage proportion
write.table(g_covered_known, file = paste0(covloc, alpha, "gaussian_covered_known.txt"), row.names = F, col.names = F)
write.table(g_covered_unknown, file = paste0(covloc, alpha, "gaussian_covered_unknown.txt"), row.names = F, col.names = F)

# D -- Compute the coverage proportion 
alpha <- c(0.05, 0.1, 0.2)

# Known or Unknown
for(a in alpha){
  cat("alpha = ", a, ": \n")
  # Known
  covered_known <- read.table(paste0(covloc, a, "_covered_known.txt"))
  # Single MLE covered
  single_covered <- rowMeans(covered_known)
  cat(a, ": single - known (null) = ", single_covered[jnull] * 100, ", sd = ", sqrt(single_covered[jnull] * (1-single_covered[jnull])) / sqrt(ncol(covered_known))  * 100 , "\n")
  cat(a, ": single - known (non-null) = ", single_covered[jnonnull] * 100, ", sd = ", sqrt(single_covered[jnonnull] * (1-single_covered[jnonnull])) / sqrt(ncol(covered_known))  * 100 , "\n")
  
  # Bulk covered
  bulk_covered <- colMeans(covered_known)
  cat(a, ": bulk - known = ", mean(bulk_covered)  * 100, ", sd = ", sd(bulk_covered)/sqrt(length(bulk_covered))  * 100, "\n")
  
  # UnKnown
  covered_unknown <- read.table(paste0(covloc, a, "_covered_unknown.txt"))
  # Single MLE covered
  single_covered <- rowMeans(covered_unknown)
  cat(a, ": single - unknown (null) = ", single_covered[jnull] * 100, ", sd = ", sqrt(single_covered[jnull] * (1-single_covered[jnull])) / sqrt(ncol(covered_unknown))  * 100 , "\n")
  cat(a, ": single - unknown (non-null) = ", single_covered[jnonnull] * 100, ", sd = ", sqrt(single_covered[jnonnull] * (1-single_covered[jnonnull])) / sqrt(ncol(covered_unknown))  * 100 , "\n")
  
  # Bulk covered
  bulk_covered <- colMeans(covered_unknown)
  cat(a, ": bulk - unknown = ", mean(bulk_covered)  * 100, ", sd = ", sd(bulk_covered)/sqrt(length(bulk_covered))  * 100, "\n")
  
  ## Coverage -- Gaussian
  covered_g_known <- read.table(paste0(covloc, a, "gaussian_covered_known.txt"))
  # Single covered
  single_covered <- rowMeans(covered_g_known)
  cat(a, ": single (Gaussian) - known (null) = ", single_covered[jnull] * 100, ", sd = ", sqrt(single_covered[jnull] * (1-single_covered[jnull])) / sqrt(ncol(covered_g_known))  * 100, "\n")
  cat(a, ": single (Gaussian) - known (nonnull) = ", single_covered[jnonnull] * 100, ", sd = ", sqrt(single_covered[jnonnull] * (1-single_covered[jnonnull])) / sqrt(ncol(covered_g_known))  * 100, "\n")
  
  # Bulk covered
  bulk_covered <- colMeans(covered_g_known)
  cat(a, ": bulk (Gaussian) - known = ", mean(bulk_covered) * 100, ", sd = ", sd(bulk_covered)/sqrt(length(bulk_covered))  * 100, "\n")
  
  # Unknown coverage
  covered_g_unknown <- read.table(paste0(covloc, a, "gaussian_covered_unknown.txt"))
  # Single covered
  single_covered <- rowMeans(covered_g_unknown)
  cat(a, ": single (Gaussian) - unknown (null) = ", single_covered[jnull] * 100, ", sd = ", sqrt(single_covered[jnull] * (1-single_covered[jnull])) / sqrt(ncol(covered_g_unknown))  * 100, "\n")
  cat(a, ": single (Gaussian) - unknown (nonnull) = ", single_covered[jnonnull] * 100, ", sd = ", sqrt(single_covered[jnonnull] * (1-single_covered[jnonnull])) / sqrt(ncol(covered_g_unknown))  * 100, "\n")
  
  # Bulk covered
  bulk_covered <- colMeans(covered_g_unknown)
  cat(a, ": bulk (Gaussian) - known = ", mean(bulk_covered) * 100, ", sd = ", sd(bulk_covered)/sqrt(length(bulk_covered))  * 100, "\n")
}





