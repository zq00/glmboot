# Simulated example in the Introduction 

# Import the resized bootstrap functions 

# Parameters 
n <- 4000 # number of observations 
p <- 400  # number of variables

rho <- 0.5
x <- rho^(c(0:(p/2), (p/2-1):1))
Sigma <- toeplitz(x)
R <- chol(Sigma) # covariance matrix

nu <- 8 # degree of freedom of the MVT distribution

beta <- scan("/scratch/users/qzhao1/logistic/boot/param/beta_large.txt") # load coef
k <- 194 # Pick a non-null coordinate
cat(beta[k]) # check the coef of the nonnull variable

# -- 1 Repeated samples
B <- 10000 # Number of repetitions
betahat_repeated <- matrix(0, B, p) # Store MLE 
for(i in 1:B){
  # Obtain one sample
  X <- matrix(rnorm(n*p, 0, 1), n, p) 
  chi <- rchisq(n, df = nu) / (nu - 2)
  X <- X %*% R / sqrt(chi) / sqrt(p)
  eta <- X %*% beta; mu <- 1 / ( 1 + exp(-eta))
  Y <- rbinom(n, 1, mu)
  fit <- glm(Y ~ X + 0, family = binomial)
  beta_hat <- fit$coef
  betahat_repeated[i, ] <- beta_hat
  
  if(i %% 10 == 0) cat(i, ",") # monitor progress
  if(i %% 200 == 0) cat("\n")
}

fileloc = "path_to_save_your_results"
write.table(betahat_repeated, file = paste0(fileloc, "intro_repeated.txt"))


# NOTE Use the same sample for nonparametric, parametric and resized bootstrap shown in Figure 1
B <- 10000 # number of bootstrap samples

X <- matrix(rnorm(n*p, 0, 1), n, p) 
chi <- rchisq(n, df = nu) / (nu - 2)
X <- X %*% R / sqrt(chi) / sqrt(p)
eta <- X %*% beta; mu <- 1 / ( 1 + exp(-eta))
Y <- rbinom(n, 1, mu)
fit <- glm(Y ~ X + 0, family = binomial)
beta_hat <- fit$coef # The MLE coef

# -- 2 Nonparametric bootstrap
mle_nonparam <- matrix(0, B, p)
for(b in 1:B){
  index <- sample(1:n, n, replace = T)
  Xb <- X[index, ]; yb <- Y[index]
  fit <- glm(yb ~ Xb + 0, family = binomial)
  mle_nonparam[b, ] <- fit$coef
  
  if(b %% 10 == 0) cat(b, ",")
  if(b %% 500 == 0) cat("\n")
}

# -- 3 Parametric bootstrap
mle_param <- matrix(0, B, p)
for(b in 1:B){
  eta <- X %*% beta_hat; mu <- 1 / ( 1 + exp(-eta))
  yb <- rbinom(n, 1, mu)
  fit <- glm(yb ~ X + 0, family = binomial)
  mle_param[b, ] <- fit$coef
  
  if(b %% 10 == 0) cat(b, ",")
  if(b %% 500 == 0) cat("\n")
}

# -- 4 Resized bootstrap
fit <- glm(Y ~ X + 0, family = binomial, x = T, y = T)
est_boot <- glm_boot(fit, b_boot = B)

fileloc = "path_to_save_your_results"
save.image(file = paste0(fileloc, "IntroEx.RData"))

