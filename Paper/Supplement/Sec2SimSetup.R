# -- 1 Sample covariates
if(distribution == "mvt"){
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

if(distribution == "gaussianint"){
  rho <- 0.5
  x <- rho^(c(0:((p - nInt)/2), ((p - nInt)/2-1):1))
  SigmaSmall <- toeplitz(x) # set the covariance matrix to be a circular matrix
  R <- chol(SigmaSmall)
  
  # Compute the Sigma matrix
  N <- 10000
  X <- matrix(rnorm(N* (p-nInt), 0, 1), N, (p-nInt)) 
  X <- X %*% R
  XInt <- apply(indices, 1, function(t) X[,t[1]] * X[,t[2]])
  X <- cbind(X, XInt)
  Sigma <- t(X) %*% X / N
  
  sample_x <- function(){
    X <- matrix(rnorm(n* (p-nInt), 0, 1), n, (p-nInt)) 
    X <- X %*% R
    XInt <- apply(indices, 1, function(t) X[,t[1]] * X[,t[2]])
    X <- cbind(X, XInt) / sqrt(p)
    
    X
  }
  
}

if(distribution == "arch"){
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

if(distribution == "pareto"){
  rho <- 0.5
  x <- rho^(c(0:(p/2), (p/2-1):1))
  Sigma <- toeplitz(x) # set the covariance matrix to be a circular matrix
  R <- chol(Sigma)
  
  sample_x <- function(){
    X <- matrix(rpareto(n*p, 5), nrow =n) %*% R / sqrt(p)
    X
  }
}

# -- 2 Sample responses
# INPUTS
# X -- n by p covariate matrix
# family -- a family object with family and link
# beta -- model coefficients
# beta0 -- intercept term
# option -- "none", "interaction" or "noise"
# ... 
sample_y <- function(X, family, beta, beta0, option = "none", ...){
  params <- list(...)
  
  if(option == "none"){
    eta <- X %*% beta + beta0
  }else if(option == "interaction"){
    if(is.null(params$betaInt)){stop("Please specify coefficients for interaction terms")}
    if(is.null(params$indices)){stop("Please specify coordinates for interaction terms ")}
    betaInt <- param$betaInt # coefficients for interactions
    IndicesInt <- param$indices # a matrix of 2 columns, storing the indices of interactions
    XInt <- apply(IndicesInt, 1, function(t) X[,t[1]] * X[,t[2]])
    
    eta <- X %*% beta + XInt %*% betaInt + beta0
  }else if(option == "noise"){
    if(is.null(params$sigma)){stop("Please specify std.dev of the additional noise")}
    sigma <- param$sigma
    eta <- X %*% beta + rnorm(n, 0, sigma) + beta0
  }
  
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

# -- 3 Set up model family and loss function (they can be different)
if(model == "logistic"){
  family_model <- get("binomial")(link = "logit")
}else if(model == "probit"){
  family_model <- get("binomial")(link = "probit")
}else if(model == "poisson"){
  family_model <- get("poisson")(link = "log")
}

if(loss == "logistic"){
  family_fit <- get("binomial")(link = "logit")
}else if(loss == "probit"){
  family_fit <- get("binomial")(link = "probit")
}else if(loss == "poisson"){
  family_fit <- get("poisson")(link = "log")
}
