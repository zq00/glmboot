# Generate model coefficients

# For logistic/ARCH model
p <- 400
nonnull <- sample(1:p, 50, replace = F)
beta <- numeric(p)
beta[nonnull] <- rnorm(50, mean = 5, sd = 1) * sample(c(-1, 1), 50, replace = T)
write.table(beta, "beta_large.txt",  col.names = F, row.names = F)

# Small
beta <- numeric(40)
nonnull = sample(1:p, 20, replace = F)
beta[nonnull] <- rnorm(20, mean = 5, sd = 1) * sample(c(-1, 1), 20, replace = T)
write.table(beta, "beta_small.txt", col.names = F, row.names = F)

# Probit/Poisson
p <- 400
nonnull <- sample(1:p, 50, replace = F)
beta <- numeric(p)
beta[nonnull] <- rnorm(50, mean = 3, sd = 1) * sample(c(-1, 1), 50, replace = T)
write.table(beta, "beta_probit.txt", col.names = F, row.names = F)

# Sparse Coefficients 
p <- 400
nonnull <- sample(1:p, 10, replace = F)
beta <- numeric(p)
beta[nonnull] <- rep(10, 10) * sample(c(-1, 1), 10, replace = T)
write.table(beta, "beta_sparse.txt", col.names = F, row.names = F)




















