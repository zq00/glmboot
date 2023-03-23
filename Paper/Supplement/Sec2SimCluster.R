# Supplement Section 2
# Simulation (on the cluster)

# -- 1. Input parameters 
args <- commandArgs(trailingOnly = T)
setting <- args[1]
k <- as.numeric(args[2]) # k goes from 1-100
jobid <- as.numeric(args[3]) # jobid goes from 1-100 in each job
index <- (k-1) * 100 + jobid

# Parameter: number of bootstrap samples
b_boot <- 1000

params <- NULL 
dir_param <- "/param/" 
file_setting <- paste0("/setting/", setting, ".R")
source(file_setting, echo = T)

# -- 2. Load the bootstrap functions 
source_code_dir <- "/R/R122022"  #The directory where all source code files are saved.
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# -- 3. Functions to sample X and Y
source('/sim_setup.R') # change the file path

# Report signal strength
gamma <- sqrt(beta %*% (Sigma %*% beta) / p)
cat("Signal strength equals to: ", gamma, "\n")


# -- Function to simulate the resized bootstrap when the signal strength is *known*
resized_boot_known <- function(fit, b_boot, s){
  betas <- s * fit$coef
  bootfit <- bootglm(fit$x, betas, family = fit$family, b_boot  = b_boot, verbose = T) # here the first argument is the covariate matrix X and the second argument is the resized bootstrap 
  
  sd <- apply(bootfit, 1, sd)
  alpha <- lm(rowMeans(bootfit) ~ betas + 0, weights = 1/sd^2)$coef # estimated inflation using the bootstrap 
  
  return(list(
    beta_s = betas,
    alpha = alpha,
    sd = sd,
    boot_sample = bootfit 
  )
  )
}

# -- Simulations 

# Create a directory to store results if the directory does not exist yet
dir_path <- paste0("/scratch/users/qzhao1/logistic/boot_revision/output/n",n,"/") 
dir_output <- setting
if(!dir.exists(dir_path)){
  dir.create(dir_path)
}
if(!dir.exists(paste0(dir_path, dir_output))){
  dir.create(paste0(dir_path, dir_output))
  dir.create(paste0(dir_path, dir_output, "/boot_known/"))
  dir.create(paste0(dir_path, dir_output, "/boot_unknown/"))
}

# sample X
X <- sample_x()
# sample Y
Y <- sample_y(X, family_model, beta, beta0, option = option, params = params)
cat("proportion of positive class is", mean(Y), "\n")
# compute the MLE and the std.dev
# I always fit an intercept when beta0 != 0 and I do not fit an intercept if beta0 = 0
if(intercept){
  fit <- glm(Y ~ X, family = family_fit, x = T, y = T)
}else{
  fit <- glm(Y ~ X + 0, family = family_fit, x = T, y = T)
}
# compute the known bootstrap
if(intercept){
  s <- gamma / sqrt(t(fit$coef[-1]) %*% (Sigma %*% fit$coef[-1]) / p)[1,1]
}else{
  s <- gamma[1,1] / sqrt(t(fit$coef) %*% (Sigma %*% fit$coef) / p)[1,1]
}

# compute the unknown bootstrap
known_boot <- resized_boot_known(fit, b_boot, s)
# Store gamma-hat, alpha-hat, sigma-hat, and the bootstrap samples from both known/unknown situations
unknown_boot <- glm_boot(fit, s_interval = 0.02, b_var = 5, b_boot = b_boot, robust_est = FALSE, verbose = TRUE, filename = NA)

mle <- fit$coef
sdR <- summary(fit)$coef[,2]
alpha_known <- known_boot$alpha
sd_known <- known_boot$sd
alpha_unknown <- unknown_boot$alpha
sd_unknown <- unknown_boot$sd
gammahat <- unknown_boot$gamma_hat
prop <- mean(Y)

# Store results 
write.table(known_boot$boot_sample, file = paste0(dir_path, dir_output, "/boot_known/boot_known_", index, ".txt"), 
            row.names = F, col.names = F)
write.table(unknown_boot$boot_sample, file = paste0(dir_path, dir_output, "/boot_unknown/boot_unknown_", index, ".txt"), 
            row.names = F, col.names = F)
write.table(unknown_boot$beta_s, file = paste0(dir_path, dir_output, "/boot_unknown/betas_unknown_", index, ".txt"),
            row.names = F, col.names = F)

write.table(cbind(mle, sdR,sd_known,sd_unknown), file = paste0(dir_path, dir_output, "/mle_",index,".txt"), 
            row.names = F, col.names = F)
write.table(c(alpha_known,alpha_unknown, gammahat, prop), 
            file = paste0(dir_path, dir_output, "/alpha_",index,".txt"), 
            row.names = F, col.names = F)




