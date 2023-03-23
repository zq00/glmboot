# Figure 6
library(tidyverse)
# Computes the theoretical relationship between gamma and var(X beta_hat)
# Here we use simulation to compute this curve - this curve can be computed theoretically and maybe i can do it later
n <- 4000
p <- 400
nInt <- p/2

# load the model parameters
beta <- scan("beta.txt") 
beta <- beta * sqrt(p) # Note: multiply beta by sqrt(p) and standardize Xj to have std. 1/sqrt(p)
indices <- read.table("int.txt",
                      header = F)

# 1. Compute the empirical curve 
# Sample X
rho <- 0.5
x <- rho^(c(0:((p - nInt)/2), ((p - nInt)/2-1):1))
SigmaSmall <- toeplitz(x) # set the covariance matrix to be a circular matrix
R <- chol(SigmaSmall)

# Compute the Sigma matrix
N <- 50000
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

sample_y <- function(X, beta){
  mu <- pnorm(X %*% beta)
  Y <- rbinom(n, 1, mu)
  Y
}

# choose a sequence of gamma
gamma0 <- sqrt(t(beta) %*% (Sigma %*% beta) / p)
gamma <- seq(0, 4, by = 0.1); ngamma <- length(gamma)

# fit the MLE each time, compute the variance
B <- 100

vals <- numeric(ngamma)
for(i in 1:ngamma){
  cat("\n", gamma[i], ":")
  # scale beta
  s <- gamma[i] / gamma0
  beta_new <- s[1,1] * beta
  new_val <- 0
  for(b in 1:B){
    X <- sample_x()
    Y <- sample_y(X, beta_new)
    fit <- glm(Y ~ X + 0, family = binomial)
    beta_hat <- fit$coef
    new_val <- new_val + sqrt(beta_hat %*% (Sigma %*% beta_hat) / p)
    
    if(b %% 10 == 0) cat(b, " ")
  }
  # take the average
  vals[i] <- new_val / B
  cat("eta = ", vals[i], "\n")
}

empirical_vals <- data.frame(gamma = gamma,
                             vals = vals)

# Store the empirical curve
write.table(empirical_vals, "curve-empirical.txt",  col.names = T,row.names = F)

# 2. Estimated curve 
# Load functions
source_code_dir <- "/home/users/qzhao1/logistic/glmhd/R010623"  #The directory where all source code files are saved.
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# function to simuulate Y from a probit model
simulate_fun_probit <- function(X, beta){
  t <- X %*% beta
  probs <- pnorm(t)
  rbinom(length(t), 1, probs)
}

# Sample a new obs. 
X <- sample_x()
Y <- sample_y(X, beta)
fit <- glm(Y ~ X + 0, family = binomial)
fit$family$simulate_fun <- simulate_fun_probit

beta_hat <- fit$coef

# Compute sd(X beta)
cat("Gamma using beta-hat is: ", sd(X %*% beta_hat), "\n")
# Observed eta
eta_obs <- estimate_eta(X, Y, fit$coef, fit$family)
cat("Observed eta is: ", eta_obs, "\n")

s_seq <- seq(0, 0.8, by = 0.02)
B <- 5 # repeat 5 times at each s
eta_hat <- matrix(NA, length(s_seq), B)
for(i in 1:length(s_seq)){
  eta_hat[i, ] <- estimate_variance(X, beta = s_seq[i] * fit$coef, fit$family, b_var = B)
  cat(i,":" , eta_hat[i, ], "\n")
}

# Store the estimated curve
result <- cbind(s_seq, eta_hat)
write.table(result, "curve-estimated.txt", row.names = F, col.names = F)

# 3. Plot the curve 
# Load the empirical curve 
curve_empirical <- read.table("curve-empirical.txt", header = T)
# Load the estimated curve
# Gamma using beta-hat is 6.31 
# Observed eta is 6.33
curve_estmated0 <- as.matrix(read.table("curve-estimated.txt"))
curve_estimated <- data.frame(
  gamma = rep(curve_estmated0[,1] * 6.31, 5),
  vals = as.vector(curve_estmated0[,-1]) 
)

# Compute a smoothed curve of the estimated eta-hat 
# Fit a loess curve
eta_obs <- 6.33
curve <- loess(vals ~ gamma, data = curve_estimated[curve_estimated$vals < 10, ])
gamma_new <- seq(0, max(curve_estimated[curve_estimated$vals < 10, ]$gamma), by = 0.01)
eta_new <- predict(curve, gamma_new) # estimated sd_hat on the smoothed loess curve
diff <- abs(eta_new - eta_obs)
gamma_hat <- gamma_new[which.min(diff)] # True gamma = 2.71

loessfit_empirical <- loess(vals ~ gamma, data = curve_empirical)

data <- tibble(
  gamma = c(gamma_new, gamma_new),
  vals = c(predict(loessfit_empirical, gamma_new), eta_new),
  Type = rep(c("Empirical", "Estimated"), each = length(gamma_new))
)
# plot the *smoothed* empirical curve
curve <- ggplot() + 
  geom_line(aes(x = gamma_new, y = predict(loessfit_empirical, gamma_new)), color = "darkorange3", size = 1) + 
  geom_point(aes(x = gamma, y = vals), data = curve_estimated, size = 0.5) +
  # geom_line(aes(x = gamma_new, y = eta_new), color = "black", size = 1) +  
  geom_abline(intercept = eta_obs, slope = 0, 
              color = "black", linetype = "dashed") + 
  xlim(c(0,3.5)) + ylim(c(0, 10)) + 
  xlab(expression(gamma[s])) + 
  ylab(expression(eta)) +
  geom_segment(aes(x = 2.8, y = 4, xend = 2.6, yend = 6),
               arrow = arrow(length = unit(0.2, "cm"))) + 
  annotate(geom = "text", x = 3, y = 3.4, 
           label = expression(paste(tilde(eta), "=6.33"))) + 
  annotate(geom = "text", x = 3, y = 2.4, 
           label = expression(paste(hat(gamma), "=2.56"))) + 
  theme_bw() + 
  theme(text = element_text(size = 18))

ggsave(plot = curve, 
       filename = paste0(figloc, "curve_incorrect.png"),
       width = 14, height = 10, units = "cm", dpi = 300)









