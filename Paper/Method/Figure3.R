# Figure 3: relationship between eta and gamma, accuracy of the estimated eta
source_code_dir <- "/R010623"  #The directory where all source code files are saved.
source_code_path <- list.files(source_code_dir, full.names = T)
for(file in source_code_path){source(file)}

# Problem setting
n <- 4000
p <- 800

# Sample coefficients 
nonnull <- sample(1:p, 50, replace = F)
beta <- numeric(p)
beta[nonnull] <- rnorm(50, mean = 5, sd = 1) * sample(c(-1, 1), 50, replace = T)

# Sample X and Y
nu <- 5
rho <- 0.5
x <- rho^(c(0:(p/2), (p/2-1):1))
cov_mat <- toeplitz(x)
R <- chol(cov_mat)

sample_x <- function(){
  X <- matrix(rnorm(n * p, 0, 1), n, p) / sqrt(p)
  chi <- rchisq(n, df = nu) / (nu - 2)
  X <- X %*% R / sqrt(chi) 
  X 
}

sample_y <- function(X, beta){
  mu <- 1 / (1 + exp(- X %*% beta))
  Y <- rbinom(n, 1, mu)
  Y
}

# True signal strength
gamma0 <- sqrt(beta %*% (cov_mat %*% beta) / p)

# 1. Empirical curve 
# A sequence of gamma
gamma <- seq(0, 3, by = 0.2); ngamma <- length(gamma)

# Fit the MLE each time, compute the variance
B <- 100
empirical_vals <- numeric(ngamma)

for(i in 1:ngamma){
  cat("\n", gamma[i], ":")
  # scale beta
  s <- gamma[i] / gamma0
  beta_new <- s * beta
  new_val <- 0
  for(b in 1:B){
    X <- sample_x()
    Y <- sample_y(X, beta_new)
    fit <- glm(Y ~ X + 0, family = binomial)
    beta_hat <- fit$coef
    new_val <- new_val + sqrt(beta_hat %*% (cov_mat %*% beta_hat) / p)
    
    if(b %% 10 == 0) cat(b, " ")
  }
  # take the average
  empirical_vals[i] <- new_val / B
}

empirical_eta <- cbind(gamma, empirical_vals)

# 2. Estimated curve
beta <- beta / gamma0[1,1] * 2 # scale beta so that gamma = 2

B <- 5 # number of repetitions 

# original data
X <- sample_x()
Y <- sample_y(X, beta)
fit <- glm(Y ~ X + 0, family = binomial)
fit$family$simulate_fun <- get_simulate_fun(fit$family) 

beta_hat <- fit$coef

# A sequence of shrinkage
s <- seq(0, 1, by = 0.02); ns <- length(s)
gamma_hat <- sd(X %*% beta_hat)

# Compute sd(X beta)
cat("Gamma using beta-hat is: ", gamma_hat, "\n")
# Observed eta
eta_obs <- estimate_eta(X, Y, fit$coef, fit$family)
cat("Observed eta is: ", eta_obs, "\n")

estimated_vals <- matrix(0, ns, B)
for(i in 1:ns){
  cat("\n", s[i], ":")
  # Scale beta_hat 
  beta_new <- s[i] * beta_hat
  estimated_vals[i, ] <- estimate_variance(X, beta = beta_new, fit$family, b_var = B)

  cat(mean(estimated_vals[i, ]))
}

# save results
dir <- "/output/"
save.image(file = paste0(dir, "curve.RData"))

# 3. Plot results 

# Compute a smoothed curve of the estimated eta-hat 
load("curve.RData")
# Estimated curve 
curve_estimated <- data.frame(
  gamma = rep(s*gamma_hat, B),
  vals = as.vector(estimated_vals)
)

loessfit_estimated <- loess(vals ~ gamma, data = curve_estimated)
gamma_new <- seq(0, 3, by = 0.01)
eta_new <- predict(loessfit_estimated, gamma_new) # estimated sd_hat on the smoothed loess curve
diff <- abs(eta_new - eta_obs)
gamma_est <- gamma_new[which.min(diff)] 

# LOESS fit of the empirical curve
empirical_eta <- as.data.frame(empirical_eta)
loessfit_empirical <- loess(empirical_vals ~ gamma, data = empirical_eta)

# plot the *smoothed* empirical curve
curve <- ggplot() + 
  geom_point(aes(x = gamma, y = vals), data = curve_estimated[curve_estimated$gamma<3, ], size = 0.5) +
  geom_line(aes(x = gamma_new, y = eta_new), color = "black", size = 1) +  
  geom_abline(intercept = eta_obs, slope = 0, 
              color = "black", linetype = "dashed") + 
  geom_line(aes(x = gamma_new, y = predict(loessfit_empirical, gamma_new)), color = "darkorange3", size = 1) + 
  xlab(expression(gamma[s])) + 
  ylab(expression(eta)) +
  geom_segment(aes(x = 2.4, y = 2.75, xend = 2, yend = 3.4),
               arrow = arrow(length = unit(0.2, "cm"))) + 
  annotate(geom = "text", x = 2.4, y = 2.5, 
           label = expression(paste(tilde(eta), "=3.48"))) + 
  annotate(geom = "text", x = 2.4, y = 2.0, 
           label = expression(paste(hat(gamma), "=1.92"))) + 
  theme_bw() + 
  theme(text = element_text(size = 18))

ggsave(plot = curve, 
       filename = paste0(figloc, "curve.png"),
       width = 14, height = 10, units = "cm", dpi = 300)








