# Simulations in Section 3 
# args <- commandArgs(trailingOnly = T)

model <- "" # set up model: logistic/probit/poisson
coef <- "" # set up model coefficient: sparse/notsparse 
covariate <- "" # set up covariate distribution: mvt/iid/arch/small

## A -- Set-up
# 1 -- Set up model family
if(model == "logistic"){
  family <- get("binomial")(link = "logit")
}else if(model == "probit"){
  family <- get("binomial")(link = "probit")
}else{
  family <- get("poisson")(link = "log")
}

# 2 -- Scan model coef 
if(model == "logistic"){
  if(coef == "sparse" & covariate == "arch"){
    n <- 4000; p <- 400
    beta <- scan("beta_sparse.txt")
  }else if(covariate == "mvt" | covariate == "arch" | covariate == "iid"){
    n <- 4000; p <- 400
    beta <- scan("beta_large.txt")
  }else if(covariate == "small"){
    n <- 400; p <- 40
    beta <- scan("beta_small.txt")
  }
}else if(model == "probit" | model == "poisson"){
  n <- 4000; p = 400
  beta <- scan("beta_probit.txt")
}

cat("model = ", model , ",covariate = ", covariate, ", p = ", p,"\n") # check model is correctly specified

# 3 -- Set up functions to sample covariates 
if(covariate == "mvt"){
  rho <- 0.5
  x <- rho^(c(0:(p/2), (p/2-1):1))
  Sigma <- toeplitz(x) # set the covariance matrix to be a circular matrix
  R <- chol(Sigma)
  
  nu <- 8 # DOF of the MVT distribution
  
  sample_x <- function(){
    X <- matrix(rnorm(n*p, 0, 1), n, p) 
    chi <- rchisq(n, df = nu) / (nu - 2)
    X <- X %*% R / sqrt(chi) / sqrt(p)
    X
  }
}

if(covariate == "arch"){
  # parameters for the modified arch model
  alpha0 <- 0.6
  alpha1 <- 0.4
  nu <- 8
  Sigma <- matrix(0, p, p)
  diag(Sigma) <- 1
  
  sample_x <- function(){
    # sample error terms
    eps <- matrix(0, n, 2 * p)
    for(i in 1:n){
      eps[i,1] <- rnorm(1, 0, sd = sqrt(alpha0 / (1-alpha1)))
      for(j in 2:(2*p)){
        eps[i, j] <- rnorm(1, 0, sd = sqrt(alpha0 + alpha1 * eps[i, j-1]^2))
      }
    }
    eps <- eps[, (p+1):(2*p)]
    
    chi <- sqrt(rchisq(n, df = nu)) 
    
    X <- eps / chi / sqrt(1/6) / sqrt(p)
    X
  }
}

rpareto <- function(n, df){
  mean <- df / (df - 1)
  sd <- sqrt(df / (df - 2) / (df - 1)^2) # standard deviation
  u <- runif(n, 0, 1)
  ((1/(1-u))^(1/df) - mean) / sd
}
if(covariate == "iid" | covariate == "small"){
  Sigma <- matrix(0, p, p)
  diag(Sigma) <- 1
  sample_x <- function(){
    X <- matrix(0, n, p)
    X <- matrix(rpareto(n*p, 5), nrow =n)
    X <- scale(X) / sqrt(p)
    X
  }
}

# Optional: Compute the signal strength parameter 
gamma <- sqrt(t(beta) %*% (Sigma %*% beta) / p) 

# 4 -- Set up function to sample response
sample_y <- function(X, family, beta){
  eta <- X %*% beta
  if(family$family == "binomial"){
    if(family$link == "logit"){
      mu <- 1 / ( 1 + exp(-eta))
      Y <- rbinom(n, 1, mu)
    }else if(family$link == "probit"){
      mu <- pnorm(eta)
      Y <- rbinom(n, 1, mu)
    }
  }else if(family$family == "poisson"){
    mu <- exp(eta)
    Y <- rpois(n, mu)
  }
  return(Y)
}


# 5 -- Set up output location
if(model == "logistic"){
  if(coef == "sparse" & covariate == "arch"){
    outloc <- "sparse/"
  }else if(covariate == "mvt"){
    outloc <- paste0("logistic/")
  }else if(covariate == "arch"){
    outloc <- paste0("arch/")
  }else if(covariate == "small"){
    outloc <- "small/"
  }else if(covariate == "iid"){
    outloc <- "iid/"
  }
}else if(model == "probit"){
  outloc <- "probit/"
}else if(model == "poisson"){
  outloc <- "poisson/"
}

## B -- Repeated simulation 
B <- 10000 
mle <- matrix(0, B, p)
sd <- matrix(0, B, p)

for(b in 1:B){
  X <- sample_x()
  Y <- sample_y(X, family, beta)
  
  fit <- glm(Y ~ X + 0, family = family)
  mle[b, ] <- fit$coef
  sd[b, ] <- summary(fit)$coef[,2]
  
  if(b %% 50 == 0) cat(b, ",")
  if(b %% 1000 == 0) cat("\n")
}

# Store results
write.table(mle, file = paste0(outloc, "mle_repeated.txt"), col.names = F, row.names = F)
write.table(sd, file = paste0(outloc, "sd_repeated.txt"), col.names = F, row.names = F)

## NOTE The remaining are based on the same observations, so the observations are simulated only once

X <- sample_x()
Y <- sample_y(X, family, beta)

# Store X, Y and beta_hat 

fit <- glm(Y ~ X + 0, family = family, x = T, y = T) # Compute MLE
beta_hat <- fit$coef

# C -- HDT 
if(model != "poisson"){
  adjusted_fit <- adjust_binary(fit, verbose = T)
  # Store estimated bias, std.dev and estimated gamma 
}

# D -- Resized bootstrap 
B <- 10000 # number of repetitions

# D.1 -- Known parameters
s <- gamma / sqrt(t(beta_hat) %*% (Sigma %*% beta_hat) / p)
beta_s <- s[1,1] * beta_hat

# Generate bootstrap samples 
boot_known <- matrix(0, B, p)
for(b in 1:B){
  yb <- sample_y(X, family, beta_s)
  fitb <- glm(yb ~ X + 0, family = family)
  boot_known[b, ] <- fitb$coef
  
  if(b %% 50 == 0) cat(b, ",")
  if(b %% 500 == 0) cat("\n")
}

# Store beta_s and the bootstrap samples 

# D.2 -- Unknown parameters 
boot_fit <- glm_boot(fit, s_interval = 0.02, b_var = 5, b_boot = B, robust_est = FALSE, verbose = TRUE, filename = NA)

# Store the resized MLE, estimated gamma, bootstrap samples

# E -- Standard bootstrap

# E.1 -- Pairs bootstrap 
mle_nonparam <- matrix(0, B, p)
for(b in 1:B){
  index <- sample(1:n, n, replace = T)
  Xb <- X[index, ]; yb <- Y[index]
  fitb <- glm(yb ~ Xb + 0, family = family)
  mle_nonparam[b, ] <- fitb$coef
  
  if(b %% 10 == 0) cat(b, ",") ; if(b %% 500 == 0) cat("\n") # keep track of progress
}

# Store bootstrap samples 

# E.2 -- Parametric bootstrap 
mle_param <- matrix(0, B, p)
for(b in 1:B){
  yb <- sample_y(X, family, beta_hat)
  fitb <- glm(yb ~ X + 0, family = family)
  mle_param[b, ] <- fitb$coef
  
  if(b %% 10 == 0) cat(b, ","); if(b %% 500 == 0) cat("\n") # keep track of progress
}

# Store bootstrap samples


