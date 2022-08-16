# Introduction example results
library(tidyverse)
library(glmhd) # Compute HDT

myloc <- "path_where_you_saved_results"
mle <- read.table(paste0(myloc, "intro_repeated.txt"))
beta <- scan("beta_large.txt")

### Bias and std.dev of the MLE  
j <- 194

# 1 -- True coef
cat("True coef = ", beta[j] , "\n") 

# 2 --  Empirical mean and std.dev of the MLE 
cat("Average MLE = ", mean(mle[,j]), "\n")
cat("sd of MLE = ", sd(mle[,j]), "\n")

# 3 -- Classical theory estimate of the std 
# Estiamte the Fisher information matrix
N <- 10000
cov <- matrix(0, p, p)
n <- 4000; p <- 400
nu <- 8
rho <- 0.5
x <- rho^(c(0:(p/2), (p/2-1):1))
Sigma <- toeplitz(x)
R <- chol(Sigma)
sample_x <- function(){
  X <- matrix(rnorm(n*p, 0, 1), n, p) 
  chi <- rchisq(n, df = nu) / (nu - 2)
  X <- X %*% R / sqrt(chi) / sqrt(p)
  
  return(X)
}
for(i in 1:N){
  Xnew <- sample_x()
  eta_new <- Xnew %*% beta
  D <- matrix(0, n, n)
  diag(D) <- 1/(1+exp(eta_new))/(1+exp(-eta_new))
  cov_new <- solve(t(X) %*% (D %*% X)) 
  cov <- cov + cov_new
  cat(i, ",")
  if(i %% 100 == 0) cat("\n")
}
CovEst <- cov / N 

diag(covEst)[j]

# 4 -- High-dimensional theory
kappa <- 0.1
p <- 400

gamma <- sqrt(beta %*% (Sigma %*% beta) / p) # signal strength
params <- find_param(kappa = kappa, gamma = gamma, intercept = F)
# The calculated theoretical parameters are below 
alpha_s <- 1.15 # HDT bias 
sigma_s <- 0.97 
tau <- 1 / sqrt(diag(solve(Sigma)))[j]
sigma_s / tau  # estimated std.dev 

# 5 -- Parametric bootstrap 
# Load the RData file you stored bootstrap results

cat("Average MLE (parametric) = ", mean(mle_param[,j]), "\n")
cat("sd of MLE (parametric) = ", sd(mle_param[,j]), "\n")

# 6 -- Pairs bootstrap 
cat("Average MLE (pairs) = ", mean(mle_nonparam[,j]), "\n")
cat("sd of MLE (pairs) = ", sd(mle_nonparam[,j]), "\n")

# 7 -- Resized bootstrap 
alpha_hat <- est_boot$alpha
cat("Estimated bias = ", alpha_hat, "\n")
cat("sd of MLE (adjusted parametric) = ", est_boot$sd[j], "\n")

### Figure 1 
density_param <- density(boot[,1])
density_pairs <- density(boot[,2])
x <- seq(0, 15, length = 512)
data <- tibble(
  type = rep(c("parametric", "pairs","resized"), each = 512), 
  loc = c(density_param$x, density_pairs$x, x),
  density = c(density_param$y, density_pairs$y, dnorm(x, mean = beta[j] * alpha_hat, sd = est_boot$sd[j])),
) 

# jpeg(filename = "IntroEx.jpg", width = 5, height = 4, units = "in", res = 150) # Uncomment to store the image
ggplot() + 
  geom_histogram(aes(x = mle[,j], y = ..density..), bins = 40, fill = "grey90", color = "grey30") +
  geom_line(data = data, aes(x = loc, y=density, group = type, color = type), 
            geom = "line", position="identity", adjust = 1.5, size = 0.8) + 
  geom_point(aes(x = beta[j], y = -0.008), shape = 17, size = 3) + 
  geom_path(aes(x = c(mean(mle[,j]), mean(mle[,j])), y = c(0, 0.3)), linetype = "dashed", size = 1) + 
  xlab(expression(hat(beta))) + 
  ylab("Density") + 
  scale_colour_manual(name  ="Method",
                      breaks=c("parametric", "pairs", "resized"),
                      values=c("blue2", "seagreen4", "red"),
                      labels=c("parametric", "pairs", "resized"),
                      guide = guide_legend(override.aes = list(linetype = c(1, 1, 1),
                                                               shape = c(16, NA, NA)))
  ) + 
  theme_bw() + 
  theme(text = element_text(size = 18, color = "black"),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
  ) 
# dev.off()



