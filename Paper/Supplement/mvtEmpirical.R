# Supplement section 1
# Compute the empirical parameters when covariates are MVT
args <- commandArgs(trailingOnly = T)
kappa <- as.numeric(args[1])
gamma <- as.numeric(args[2]) 
nu <- as.numeric(args[3])  # DOF of the MVT distribution

n <- 4000
p <- n * kappa
sample_obs <- function(){
  X <- matrix(rnorm(n*p, 0, 1), n, p) 
  chi <- rchisq(n, df = nu) / (nu - 2)
  X <- X / sqrt(chi) 
  Y <- rbinom(n, 1, 1 / ( 1 + exp(-X[,1] * gamma)))
  
  return(list(X = X, Y = Y))
}
f2 <- function(t) 1/(1 + exp(t)) / (1+exp(-t))

B <- 500
vals <- matrix(NA, B, 3)
for(b in 1:B){
  obs <- sample_obs()
  fit <- glm(obs$Y~obs$X + 0, family = binomial)
  vals[b, 1] <- fit$coef[1] / gamma
  vals[b, 2] <- sqrt(sum(fit$coef[-1]^2)) 
  D <- matrix(0, n, n)
  diag(D) <- f2(predict(fit))
  G <- t(obs$X) %*% D %*% obs$X 
  vals[b, 3] <- sum(diag(solve(G)))
  
  cat(b, ":", vals[b,], "\n")
}

dir <- "/theory-empirical/"
write.table(vals, paste0(dir, nu, "_",kappa, "_", gamma, ".txt"),
            row.names = F, col.names = c("a", "s", "l"))

cat("inflation:", mean(vals[,1]), "\n")
cat("sigma:", mean(vals[,2]), "\n")
cat("lambda:", mean(vals[,3]), "\n")





