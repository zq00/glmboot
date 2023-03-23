# Setting 4
# Setting up model parameters for the simulation
filename <- "setting4.txt"

n <- 4000
kappa <- 0.2
p <- n * kappa

beta0 <- 1.5

sparsity <- 0.10

mu <- 6
sd <- 2
ps <- 0.8
option <- "none"

distribution <- "pareto"
loss <- "logistic"
model <- "logistic"
intercept <- TRUE

# load parameters
if(file.exists(paste0(dir_param, filename))){
  beta <- scan(paste0(dir_param, filename))
}else{
  nonnull <- sample(1:p, sparsity * p, replace = F)
  beta <- numeric(p)
  beta[nonnull] <- rnorm(length(nonnull), mean = mu, sd = sd) * sample(c(1, -1), length(nonnull), replace = T, prob = c(ps, 1-ps)) 
  write.table(beta, file = paste0(dir_param, filename), col.names = F, row.names = F)
}

