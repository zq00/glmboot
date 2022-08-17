# Computes Table 1 (bias and std.dev of the MLE)

model = "" # Set model: logistic/probit/poisson
covariate =  "" # Set covariate distribution: mvt/arch/small
coef = "" # Set sparse/notsparse coefficient

# A --  Set up
# A.1 Load model coefficients
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

# A.2 Output location
if(model == "logistic"){
  if(coef == "sparse"){
    jnull <- ""; jnonnull <- "" # set null/nonnull variable coordinate
    fileloc <- "/output/sparse/"
    outloc <- "/result/sparse/"
  }else if(covariate == "mvt"){
    jnull <- ""; jnonnull <- "" # set null/nonnull variable coordinate
    fileloc <- "/output/logistic/"
    outloc <- "/result/logistic/"
  }else if(covariate == "arch"){
    jnull <- ""; jnonnull <- "" # set null/nonnull variable coordinate
    fileloc <- "/output/arch/"
    outloc <- "/result/arch/"
  }else if(covariate == "small"){
    fileloc <- "/output/small/"
    outloc <- "/result/small/"
  }
}else if(model == "probit"){
  jnull <- ""; jnonnull <- "" # set null/nonnull variable coordinate
  fileloc <- "/output/probit/"
  outloc <- "/result/probit/"
}else if(model == "poisson"){
  jnull <- ""; jnonnull <- "" # set null/nonnull variable coordinate
  fileloc <- "/output/poisson/"
  outloc <- "/result/poisson/"
}

# A.3 List all the outputs 
names <- list.files(fileloc)
indices <- unique(sapply(names, function(t) strsplit(t, "_")[[1]][2]))
N <- length(indices)

# B -- Compute values in the table
# B.1 known parameters
alpha_known <- NULL
sigma_known <- NULL

for(i in 1:N){
  beta_s <- tryCatch(read.table(paste0(fileloc, "betas_known_", indices[i]))[,1], error = function(e) {return(-1)})
  boot <- tryCatch(read.table(paste0(fileloc, "boot_known_", indices[i])), error = function(e) {return(-1)})
  
  if(length(beta_s) == 1 | length(boot) == 1 ) next;
  
  # compute std
  std_boot <- apply(boot, 2, sd)
  alpha_boot <- lm(colMeans(boot) ~ beta_s + 0, weights = 1/std_boot^2)$coef
  
  alpha_known <- c(alpha_known, alpha_boot)
  sigma_known <- rbind(sigma_known, std_boot[c(jnull, jnonnull)])
  
  cat(i, ",")
}

write.table(alpha_known, paste0(outloc, "alpha_known.txt"))
write.table(sigma_known, paste0(outloc, "sigma_known.txt"))

# unknown parameters
alpha_unknown <- NULL
sigma_unknown <- NULL

for(i in 1:N){
  beta_s <- tryCatch(read.table(paste0(fileloc, "betas_unknown_", indices[i]))[,1], error = function(e) {return(-1)})
  boot <- tryCatch(read.table(paste0(fileloc, "boot_unknown_", indices[i])), error = function(e) {return(-1)})
  
  if(length(beta_s) == 1 | length(boot) == 1 ) next;
  
  # compute std
  std_boot <- apply(boot, 1, sd)
  alpha_boot <- lm(rowMeans(boot) ~ beta_s + 0, weights = 1/std_boot^2)$coef
  
  alpha_unknown <- c(alpha_unknown, alpha_boot)
  sigma_unknown <- rbind(sigma_unknown, std_boot[c(jnull, jnonnull)])
  
  cat(i, ",")
}

write.table(alpha_unknown, paste0(outloc, "alpha_unknown.txt"))
write.table(sigma_unknown, paste0(outloc, "sigma_unknown.txt"))

# C -- Compute tables 
alpha_known <- read.table(paste0(outloc, "alpha_known.txt"))
colMeans(alpha_known); sd(alpha_known[,1])

alpha_unknown <- read.table(paste0(outloc, "alpha_unknown.txt"))
colMeans(alpha_unknown); sd(alpha_unknown[,1])

sigma_known <- read.table(paste0(outloc, "sigma_known.txt"))
colMeans(sigma_known); apply(sigma_known, 2, sd)

sigma_unknown <- read.table(paste0(outloc, "sigma_unknown.txt"))
colMeans(sigma_unknown); apply(sigma_unknown, 2, sd)

# D -- Theoretical and empirical bias and std.dev

mle <- read.table("mle_repeated.txt") # Read MLE from repeated simulations
sd <- read.table("sd_repeated.txt") # read estimated std.dev using glm function (classical theory)

mean(mle[,jnull]); sd(mle[,jnull])
mean(mle[,jnonnull]); sd(mle[,jnonnull]) # empirical std.dev
mean(sd[,jnull]); mean(sd[,jnonnull]) # classical theory

# HDT
gamma <- sqrt(t(beta) %*% beta / p) # Load beta first if needed 
param <- find_param(kappa = 0.1, 
                    gamma = gamma, 
                    beta0=0,
                    intercept = F) # Logistic 

rho_prime <- function(t) pnorm(t)
f_prime1 <- function(t) - dnorm(t) / pnorm(t)
f_prime0 <- function(t) dnorm(t) / pnorm(-t)
param <- find_param(kappa = 0.1, 
                    gamma = gamma, 
                    beta0=0,
                    intercept = F,
                    rho_prime = rho_prime,
                    f_prime1 = f_prime1,
                    f_prime0 = f_prime0) # Probit











