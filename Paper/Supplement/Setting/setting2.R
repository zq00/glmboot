# Setting 2
# Setting up model parameters for the simulation
filename <- "setting2.txt"

n <- 4000
kappa <- 0.3
p <- n * kappa

beta0 <- 0

sparsity <- 0.25

mu <- 3
sd <- 2
ps <- 0.5
option <- "none"

distribution <- "arch"
loss <- "logistic"
model <- "logistic"
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

