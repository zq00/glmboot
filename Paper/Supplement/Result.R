# Simulation results 
# Parameters (for running the code on the cluster)
args <- commandArgs(trailingOnly = T)
setting <- args[1] # simulation setting
jnull <- as.numeric(args[2]) # coordinate of a null or non-null variable
jnonnull <- as.numeric(args[3])

# Obtaining gamma and Sigma
dir_param <- "/param/" 
file_setting <- paste0("/setting/", setting, ".R")
source(file_setting, echo = T)
source('/sim_setup.R') 

# Report signal strength
gamma <- sqrt(beta %*% (Sigma %*% beta) / p)
cat("Signal strength equals to: ", gamma, "\n")
kappa <- p/n
cat("Problem proportion equals to: ", kappa, "\n")

# Print the null and nonnull coordinate
cat("A null coordinate:", jnull, "\n")
cat("A NON null coordinate:", jnonnull, "; the effect size is ",
    beta[jnonnull],"\n")

if(intercept) {
  p <- p + 1
  jnull <- jnull + 1
  jnonnull <- jnonnull + 1
  if(setting == "logistic_18") beta <- c(2, beta)
  if(setting == "logistic_18s") beta = c(1.5, beta) # add an intercept term
}

# -- 1. Inflation and standard deviation of a single coordinate 
# Load the estimated inflation and std.dev
dir_path <- "" # path to store results
# Note: fit an intercept when beta0 != 0 and do NOT fit an intercept if beta0 = 0

alpha_names <- list.files(path = paste0(dir_path, dir_output), 
                          pattern = "alpha_")
indices <- sapply(alpha_names, function(t) strsplit(t, "_")[[1]][2]) # get the indices for each simulation

N <- length(alpha_names)

# The matrix alpha_result has four columns:
# 1- alpha_known: estimated alpha assuming known signal strength
# 2- alpha_unknown: resized bootstrap estimate
# 3- gammahat: estimated signal strength
# 4- prop: proportions of positive class
alpha_result <- matrix(NA, nrow = N, ncol = 4)
for(i in 1:N){
  alpha_result[i, ] <- scan(paste0(dir_path,dir_output,"/alpha_",indices[i]))
}


mle <- matrix(NA, ncol = N, nrow = p) # each row represents one variable
sdR<- matrix(NA, ncol = N, nrow = p)
sd_known<- matrix(NA, ncol = N, nrow = p)
sd_unknown<- matrix(NA, ncol = N, nrow = p)

for(i in 1:N){
  cat(i, ",")
  # The matrix MLE has four columns (each file has four columns)
  # 1- mle: the MLE coef
  # 2- sdR: the estimated std.dev using R
  # 3 - sd_known: the estimated std.dev using known signal strength
  # 4 - sd_unknown: the estimated std.dev assuming unknown (and thus estimated) signal strength 
  newdata <- read.table(paste0(dir_path,dir_output,"/mle_",indices[i]))
  
  mle[, i] <- newdata[,1]
  sdR[, i] <- newdata[,2]
  sd_known[, i] <- newdata[,3]
  sd_unknown[, i] <- newdata[,4]
}

# 1.1 Estimated gamma versus true gamma
# estimated gammas
gamma_hat <- alpha_result[,3]

cat("True signal strength is ", gamma, "\n")
cat("Estimated signal strength: Average = ", mean(gamma_hat),
    "; Std.Dev (of the mean) = ", sd(gamma_hat)/sqrt(length(gamma_hat)),"\n")

# 1.2 Estimated inflation

# Compute the estimated inflation and std.dev from high-dimensional theory
rho_prime_probit <- function(t) pnorm(t)
fprime_1_logistic <- function(x) -1 / (1 + exp(x))
fprime_0_logistic <- function(x) 1 / (1 + exp(-x))
fprime_1_probit <- function(t) - dnorm(t) / pnorm(t)
fprime_0_probit <- function(t) dnorm(t) / pnorm(-t)
# The code below computes the parameters in section 5 (when we fit a logistic regression for probit model)
# params <- find_param(kappa = kappa, gamma = gamma, intercept = F,
#                      rho_prime = rho_prime_probit, f_prime1 = fprime_1_logistic,
#                      f_prime0 = fprime_0_logistic)
params <- find_param(kappa = kappa, gamma = gamma, intercept = F,
                     rho_prime = rho_prime_probit, f_prime1 = fprime_1_probit,
                     f_prime0 = fprime_0_probit, verbose = T)
# params <- find_param(kappa = kappa, gamma = gamma, intercept = F) # Use this code for a logistic regression (true model is also logistic)
# params <- find_param(kappa = kappa, gamma = gamma, beta0=1.5, intercept = T) # Use this for the setting when there is an intercept in the model

sd_highdim <- params[3] * sqrt(diag(solve(Sigma)))

cat("A nonnull coordinate is", jnonnull, "\n")
cat("The effect size is", beta[jnonnull], "\n")

cat("The observed inflation is ", mean(mle[jnonnull,])/beta[jnonnull], "\n")
cat("Estimated inflation using High-dim Theory is ", params[1], "\n")
cat("Estimated inflation using known signal strength is ", mean(alpha_result[,1]), 
    "and its std.dev is ", sd(alpha_result[,1])/sqrt(nrow(alpha_result)), "\n")
cat("Estimated inflation using resized bootstrap is ", mean(alpha_result[,2]), 
    "and its std.dev is ", sd(alpha_result[,2])/sqrt(nrow(alpha_result)), "\n")

# 1.3 Estimated std.dev
# sd of null 
cat("Empircal std.dev is ", sd(mle[jnull,]), "\n")
cat("Estimated std.dev using R is ", mean(sdR[jnull,]), "\n")
cat("Estimated std.dev using HIGH-DIM THEORY is ", sd_highdim[jnull], "\n")
cat("Estimated STD.DEV using known signal strength is ", mean(sd_known[jnull,]), 
    "and its std.dev is ", sd(sd_known[jnull,])/sqrt(ncol(sd_known)), "\n")
cat("Estimated STD.DEV using resized bootstrap is ", mean(sd_unknown[jnull,]), 
    "and its std.dev is ", sd(sd_unknown[jnull,])/sqrt(ncol(sd_unknown)), "\n")

# sd of nonnull
cat("Empircal std.dev of a NONNULL is ", sd(mle[jnonnull,]), "\n")
cat("Estimated std.dev of a NONNULL using R is ", mean(sdR[jnonnull,]), "\n")
cat("Estimated std.dev of a NONNULLusing HIGH-DIM THEORY is ", sd_highdim[jnonnull], "\n")
cat("Estimated STD.DEV of a NONNULLusing known signal strength is ", mean(sd_known[jnonnull,]), 
    "and its std.dev is ", sd(sd_known[jnonnull,])/sqrt(ncol(sd_known)), "\n")
cat("Estimated STD.DEV of a NONNULL using resized bootstrap is ", mean(sd_unknown[jnonnull,]), 
    "and its std.dev is ", sd(sd_unknown[jnonnull,])/sqrt(ncol(sd_unknown)), "\n")


# -- 2. Coverage proportion 

## Functions to compute the coverage proportion
# covered is the matrix of 0 or 1 meaning covered or not, where each COLUMN stores z-scores from 1 repetition
# j is the index we want to focus on 
compute_coverage <- function(covered, j, alpha){
  cat("Bulk coverage: Mean = ", mean(colMeans(covered)) * 100,
      "STD = ", sd(colMeans(covered))/sqrt(ncol(covered)) * 100, "\n")
  cat("Individual coordinate coverage: Mean = ", mean(covered[j, ]) * 100,"STD = ", sd(covered[j, ])/sqrt(ncol(covered)) * 100, "\n")
  
  return(list(bulk_mean = mean(colMeans(covered)) * 100,
              bulk_sd =  sd(colMeans(covered))/sqrt(ncol(covered))* 100,
              single_mean = mean(covered[j, ])* 100,
              single_sd = sd(covered[j, ])/sqrt(ncol(covered))* 100))
}

compute_coverage_z <- function(z, j, alpha) # compute the coverage proportion from z-scores
{
  covered <- z < qnorm((1+alpha)/2) & z > qnorm((1-alpha)/2)
  compute_coverage(covered, j, alpha)
}

# 2.1 Theoretical
all_alpha <- c(0.95, 0.9, 0.8)
# 2.1.1 Classical
z_classical <- (mle - beta) / sdR # each column is Z score from 1 repetition 

for(alpha in all_alpha){
  compute_coverage_z(z_classical, jnull, alpha)
  compute_coverage_z(z_classical, jnonnull, alpha)
}


# 2.1.2 High-dimensional
alpha_s <- params[1]
z_highdim <- (mle - alpha_s * beta) /  sd_highdim # each column is Z score from 1 repetition 

for(alpha in all_alpha){
compute_coverage_z(z_highdim, jnull, alpha)
compute_coverage_z(z_highdim, jnonnull, alpha)
}

# 2.2 Resized bootstrap
# 2.2.1 Known Gamma
# 2.2.1.1 Boot-g
z_boot_known <- (mle - beta %*% t(alpha_result[,1])) / sd_known

for(alpha in all_alpha){
compute_coverage_z(z_boot_known, jnull, alpha)
compute_coverage_z(z_boot_known, jnonnull, alpha)
}
# 2.2.1.2 Boot-t
# Compute the quantile of the t-statistics 

compute_coverage_t <- function(j, all_alpha, option){
  # compute the lower and upper confidence bounds
  lower <- list() # the lower and upper confidence bounds are lists --> each element in the list stores the lower/upper bound for one alpha 
  upper <- list()
  
  for(a in 1:length(all_alpha)){
    lower[[a]] <- matrix(NA, p, N)
    upper[[a]] <- matrix(NA, p, N)
  }
  if(option == "boot_known"){
    for(i in 1:N){
      new_boot_mle <- read.table(paste0(dir_path,dir_output,"/boot_known/boot_known_", indices[i]))
      # Get beta_s
      if(intercept){
        s <- gamma[1,1] / sqrt(t(mle[-1,i]) %*% (Sigma %*% mle[-1,i]) / p)[1,1]
      }else{
        s <- gamma[1,1] / sqrt(t(mle[,i]) %*% (Sigma %*% mle[,i]) / p)[1,1]
      }
      beta_s <- mle[,i] * s
      
      new_t <- (new_boot_mle - alpha_result[i, 1] * beta_s) / sd_known[,i] # each row of sd represents the estimated std.dev in one experiment
      
      for(a in 1:length(all_alpha)){ # each column represents the lower bound in one experiment
        lower[[a]][,i] <- (mle[,i] - sd_known[,i] * apply(new_t, 1, function(t) quantile(t, probs = (1 + all_alpha[a]) / 2))) / alpha_result[i, 1]
        upper[[a]][,i] <- (mle[,i] - sd_known[,i] * apply(new_t, 1, function(t) quantile(t, probs = (1 - all_alpha[a]) / 2))) / alpha_result[i, 1]
      }
    }
  }
  
  if(option == "boot_unknown"){
    for(i in 1:N){
      new_boot_mle <- read.table(paste0(dir_path,dir_output,"/boot_unknown/boot_unknown_", indices[i]))
      beta_s <- read.table(paste0(dir_path,dir_output,"/boot_unknown/betas_unknown_", indices[i]))
      beta_s <- as.matrix(beta_s)
     
      new_t <- (new_boot_mle - alpha_result[i, 2] * beta_s) / sd_unknown[ ,i] # each row of sd represents the estimated std.dev in one experiment
      
      for(a in 1:length(all_alpha)){ # each column represents the lower bound in one experiment
        lower[[a]][,i] <- (mle[,i] - sd_unknown[,i] * apply(new_t, 1, function(t) quantile(t, probs = (1 + all_alpha[a]) / 2))) / alpha_result[i, 2]
        upper[[a]][,i] <- (mle[,i] - sd_unknown[,i] * apply(new_t, 1, function(t) quantile(t, probs = (1 - all_alpha[a]) / 2))) / alpha_result[i, 2]
      }
    }
  }
  
  result <- list()
  for(a in 1:length(all_alpha)){
    covered <- lower[[a]] < beta & upper[[a]] > beta
    result[[a]] <- compute_coverage(covered, j, all_alpha[a])
  }
  
  return(result)
}

cat("Coverage for a NULL coordinate (known gamma) \n")
compute_coverage_t(jnull, all_alpha, option = "boot_known")

cat("Coverage for a NONNULL coordinate (known gamma) \n")
compute_coverage_t(jnonnull, all_alpha, option = "boot_known")

# 2.2.2 Unknown Gamma
# 2.2.2.1 Boot-g
z_boot_unknown <- (mle - beta %*% t(alpha_result[,2])) / sd_unknown

for(alpha in all_alpha){
compute_coverage_z(z_boot_unknown, jnull, alpha)
compute_coverage_z(z_boot_unknown, jnonnull, alpha)
}
# 2.2.2.2 Boot-t
cat("Coverage for a NULL coordinate (UNKNOWN gamma) \n")
compute_coverage_t( jnull, all_alpha, option = "boot_unknown")
cat("Coverage for a NONNULL coordinate (UNKNOWN gamma) \n")
compute_coverage_t( jnonnull, all_alpha, option = "boot_unknown")






