# Probit model
# Setting up model parameters for the simulation
filename <- "probit.txt"

n <- 4000
kappa <- 0.1
p <- n * kappa

beta0 <- 0

sparsity <- 0.125

mu <- 3
sd <- 1
ps <- 0.5 # proportion of positive beta
option <- "none"

distribution <- "arch"
loss <- "probit"
model <- "probit"
intercept <- FALSE

# load parameters
if(file.exists(paste0(dir_param, filename))){
  beta <- scan(paste0(dir_param, filename))
}else{
  nonnull <- sample(1:p, sparsity * p, replace = F)
  beta <- numeric(p)
  beta[nonnull] <- rnorm(length(nonnull), mean = mu, sd = sd) * sample(c(1, -1), length(nonnull), replace = T, prob = c(ps, 1-ps)) 
  write.table(beta, file = paste0(dir_param, filename), col.names = F, row.names = F)
}
