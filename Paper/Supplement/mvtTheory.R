# Supplement section 1 
# Parameters alpha_star and sigma_star when the covariates are from a MVT
library(cubature)
library(pracma)

# Input the problem dimension kappa and the signal strength theta1
args <- commandArgs(trailingOnly = T)
kappa <- as.numeric(args[1]) 
theta1 <- as.numeric(args[2])
nu <- 8 # degree of freedom of the chi-squared distribution

# Function to compute the proximal operator 
prox_op <- function(f_prime, lambda, x){
  f <- function(z)  z + lambda * f_prime(z) - x
  uniroot(f, interval = c(-20,20), extendInt = "yes")$root
}

rho_prime_logistic <- function(x) 1 / (1 + exp(-x))
f_prime1_logistic <- function(x) -1 / (1 + exp(x))
f_prime0_logistic <- function(x) 1 / (1 + exp(-x))
f2 <- function(t) 1/(1 + exp(t)) / (1+exp(-t)) # f''

# At each alpha, compute the values of the last two equations 
g <- function(y, a){
  sigma <- y[1]
  lambda <- y[2]
  f <- function(x){
    z1 <- x[1]
    z2 <- x[2]
    t <- x[3]
    
    u <- sqrt((nu - 2) / t) # t is from chi-squared with nu degree of freedom
    # probability Y = 1
    rho <- 1/(1 + exp(- u * z1 * theta1)) # probability of y = 1 
    w <- u^2 * lambda
    eta1 <- prox_op(f_prime1_logistic, w, u * z1 * a * theta1 + u * sigma * z2)
    eta0 <- prox_op(f_prime0_logistic, w, u * z1 * a * theta1 + u * sigma * z2)
    
    integrand_1 <- lambda^2 * (f_prime1_logistic(eta1)^2 * u^2 * rho + f_prime0_logistic(eta0)^2 * u^2 * (1 - rho))
    integrand_2 <- 1 / (1 + f2(eta1) * u^2 * lambda) * rho + 
      1 / (1 + f2(eta0) * u^2 * lambda) * (1 - rho)
    
    c(integrand_1 * dnorm(z1, 0, 1) * dnorm(z2, 0, 1) * dchisq(t, nu),
      integrand_2 * dnorm(z1, 0, 1) * dnorm(z2, 0, 1) * dchisq(t, nu))
  }
  
  val <- hcubature(f, lowerLimit = c(-8, -8, 0), upperLimit = c(8, 8, nu * 5), 
                   fDim = 2, tol = 1e-5)$integral
  
  cat("sigma = ", sigma, ", lambda = ", lambda, ", function values are", round(sigma^2 * kappa - val[1], 2),
      "and", round(1 - kappa - val[2],2) ,"\n")
  c(sigma^2 * kappa - val[1], 
    1 - kappa - val[2])
}

# compute the first equation (setting the first eqn = 0 solves alpha_star)
compute_score <- function(alpha, sigma0, lambda0){
  # solve sigma, lambda for an alpha based on teh
  vals <- fsolve(g, c(sigma0, lambda0), a = alpha, tol = 0.001)$x
  # compute the first equation
  sigma <- vals[1] 
  lambda <- vals[2]
  f <- function(x){
    z1 <- x[1]
    z2 <- x[2]
    t <- x[3]
    
    u <- sqrt((nu - 2) / t) # t is from chi-squared with nu degree of freedom
    # probability Y = 1
    rho <- 1/(1 + exp(- u * z1 * theta1))
    w <- u^2 * lambda
    eta1 <- prox_op(f_prime1_logistic, w, u * z1 * alpha * theta1 + u * sigma * z2)
    eta0 <- prox_op(f_prime0_logistic, w, u * z1 * alpha * theta1 + u * sigma * z2)
    
    integrand1 <- u * z1 * f_prime1_logistic(eta1) * rho 
    integrand2 <- u * z1 * f_prime0_logistic(eta0) * (1 - rho)
    
    c(
      integrand1 * dnorm(z1, 0, 1) * dnorm(z2, 0, 1) * dchisq(t, nu),
      integrand2 * dnorm(z1, 0, 1) * dnorm(z2, 0, 1) * dchisq(t, nu)
    )
  }
  
  score <- hcubature(f, lowerLimit = c(-8, -8,0), upperLimit = c(8, 8, nu * 5), 
                     fDim = 2, tol = 1e-4)$integral
  list(score = score[1] + score[2], sigma = sigma, lambda = lambda)
}

# Compute the first eqn for a sequence of alpha values
alphas <- seq(1.6, 3.0, by = 0.01)
score <- numeric(length(alphas))
sigma0 <- 2;lambda0 <- 2 # initialize sigma and lambda
results <- matrix(NA, nrow = length(alphas), ncol = 4)
for(i in 1:length(alphas)){
  val <- compute_score(alphas[i], sigma0, lambda0)
  score[i] <- val$score
  sigma0 <- val$sigma
  lambda0 <- val$lambda
  cat(c(alphas[i], score[i]), "\n")
  results[i, ] <- c(alphas[i], sigma0, lambda0, score[i])
}

# store results 
dir <- "/theory/"
write.table(results, 
            file = paste0(dir, kappa, "_", theta1, ".txt"),
            row.names = F, col.names = c("a", "s", "l", "score"))
