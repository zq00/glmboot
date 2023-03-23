# 1/6/2023
args <- commandArgs(trailingOnly = T)
setting <- args[1]
jnull <- as.numeric(args[2])
jnonnull <- as.numeric(args[3])


# Obtaining gamma and Sigma
dir_param <- "/param/" 
file_setting <- paste0("/setting/", setting, ".R")
source(file_setting, echo = T)
source('/sim_setup.R') # change the file path


# Report signal strength
gamma <- sqrt(beta %*% (Sigma %*% beta) / p)
cat("Signal strength equals to: ", gamma, "\n")
kappa <- p/n
cat("Problem proportion equals to: ", kappa, "\n")

# User input a null and a nonnull coordiate
cat("A null coordinate:", jnull, "\n")
cat("A NON null coordinate:", jnonnull, "; the effect size is ",
    beta[jnonnull],"\n")

if(intercept) {
  p <- p + 1
  jnull <- jnull + 1
  jnonnull <- jnonnull + 1
  if(setting == "logistic_18") beta <- c(2, beta)
  if(setting == "logistic_18s") beta <- c(1.5, beta)
}

# -- 1. Inflation and standard deviation of a single coordinate 
# Load the estimated inflation and std.dev
dir_path <- paste0("/scratch/users/qzhao1/logistic/boot_revision/output/n",n,"/") # path to store results
# create a directory to store results for n = 1000 (upper level directory)
# I always fit an intercept when beta0 != 0 and I do not fit an intercept if beta0 = 0
if(setting == "logistic_51"){
  dir_output <- "logistic_51-correct"
}else if(setting == "logistic_9"){
  dir_output <- "logistic_9-correct"
}else{
  dir_output <- setting
}


alpha_names <- list.files(path = paste0(dir_path, dir_output), 
                          pattern = "alpha_")
indices <- sapply(alpha_names, function(t) strsplit(t, "_")[[1]][2])

N <- length(alpha_names)
# get the indices

alpha_result <- matrix(NA, nrow = N, ncol = 4)
for(i in 1:N){
  alpha_result[i, ] <- scan(paste0(dir_path,dir_output,"/alpha_",indices[i]))
}
# The matrix alpha_result has four columns:
# 1- alpha_known: estimated alpha assuming known signal strength
# 2-alpha_unknown: resized bootstrap estimate
# 3-gammahat: estimated signal strength
# 4- prop: proportions of positive class
mle <- matrix(NA, ncol = N, nrow = p)
sdR<- matrix(NA, ncol = N, nrow = p)
sd_known<- matrix(NA, ncol = N, nrow = p)
sd_unknown<- matrix(NA, ncol = N, nrow = p)

# The matrix MLE has four columns (each file has four columns)
# 1- mle: the first columns is the MLE
# 2- sdR: the second columns is the estimated std.dev using R
# 3 - sd_known: the third column is the estimated std.dev using known signal strength
# 4 - sd_unknown: the fourth column is the estimated std.dev assuming unknown signal strenth 
for(i in 1:N){
  cat(i, ",")
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

# -- 2. Coverage proportion 

## Functions to compute the coverage proportion
# covered is the matrix of 0 or 1 meaning covered or not, where each COLUMN stores z-scores from 1 repetition
# j is the index we want to focus on 
compute_coverage <- function(covered, j, alpha){
  cat("Bulk coverage: Mean = ", mean(colMeans(covered)) * 100,"STD = ", sd(colMeans(covered))/sqrt(ncol(covered)) * 100, "\n")
  cat("Individual coordinate coverage: Mean = ", mean(covered[j, ]) * 100,"STD = ", sd(covered[j, ])/sqrt(ncol(covered)) * 100, "\n")
  
  return(list(bulk_mean = mean(colMeans(covered)) * 100,
              bulk_sd =  sd(colMeans(covered))/sqrt(ncol(covered))* 100,
              single_mean = mean(covered[j, ])* 100,
              single_sd = sd(covered[j, ])/sqrt(ncol(covered))* 100))
}

compute_coverage_z <- function(z, j, alpha) # compute the covarege proportion from z-scores
{
  covered <- z < qnorm((1+alpha)/2) & z > qnorm((1-alpha)/2)
  compute_coverage(covered, j, alpha)
}

# 2.1 Theoretical
# 2.1.1 Classical


# 2.2 Resized bootstrap
# 2.2.1 Known Gamma

# 2.2.1.2 Boot-t
# Compute the quantiles of the t-statistics 
# Input all of the alphas here
all_alpha <- c(0.8, 0.9, 0.95)
# for each alpha, I need to compute quantile from the standardized t-statistics
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
      cat(i, ",")
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
      cat(i, ",")
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
# 2.2.2.2 Boot-t
cat("Coverage for a NULL coordinate (UNKNOWN gamma) \n")
compute_coverage_t( jnull, all_alpha, option = "boot_unknown")
cat("Coverage for a NONNULL coordinate (UNKNOWN gamma) \n")
compute_coverage_t( jnonnull, all_alpha, option = "boot_unknown")






