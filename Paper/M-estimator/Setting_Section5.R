# LOGISTIC 51
# Setting up model parameters for the simulation
filename <- "logistic_51.txt"

n <- 4000
kappa <- 0.1
p <- n * kappa

beta0 <- 0

sparsity <- 0.25

mu <- 0
sd <- 0.25 # Note: in this simulation I did not standardize X_j, thus if we consider
# X_j to be standardized, the sd = 0.25 * sqrt(400) = 5
ps <- 0.5 # proportion of positive beta
option <- "none"

distribution <- "gaussianint"
loss <- "logistic"
model <- "probit"
intercept <- FALSE

nInt <- p / 2 # Number of interactions

# load parameters
if(file.exists(paste0(dir_param, filename))){
  beta <- scan(paste0(dir_param, filename))
  indices <- read.table(paste0(dir_param, "int_",filename), 
                        header = F)
}else{
  indices <- matrix(sample(1:(p-nInt), 2000,replace = T), ncol = 2)
  indices <- t(apply(indices, 1, sort))
  indices <- unique(indices)[1:nInt, ]
  
  nonnull <- sample(1:p, sparsity * p, replace = F)
  beta <- numeric(p)
  beta[nonnull] <- rnorm(length(nonnull), mean = mu, sd = sd) * sample(c(1, -1), length(nonnull), replace = T, prob = c(ps, 1-ps)) 
  write.table(beta, file = paste0(dir_param, filename), col.names = F, row.names = F)
  write.table(indices, file = paste0(dir_param, "int_",filename), col.names = F, row.names = F)
}
