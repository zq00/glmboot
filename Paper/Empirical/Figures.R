# Generate figures in Section 4 
library(tidyverse)
library(data.table)


model <- "" # Set up models: logistic/probit/poisson
covariate <- "" # Set up covariate distribution: mvt/arch/small
coef <- "" # sparse/notsparse coefficients
p <- 400 # Number of variables

# -- Set up
# Set up covariance matrix 
if(covariate == "mvt"){
  rho <- 0.5
  x <- rho^(c(0:(p/2), (p/2-1):1))
  Sigma <- toeplitz(x)
  R <- chol(Sigma)
}

if(covariate == "arch" | covariate == "small"){
  Sigma <- matrix(0, p, p)
  diag(Sigma) <- 1
}

# Load model coefficients
if(covariate == "mvt" | covariate == "arch"){
  beta <- scan("beta_large.txt")
}else if(covariate == "small"){
  beta <- scan("beta_small.txt")
}

gamma <- sqrt(t(beta) %*% (Sigma %*% beta) / p) # Compute the signal strength 

# 1 -- Bias of the MLE when covariates are from a modified ARCH model (Figure 4)
mle <- fread("/mle_repeated.txt", data.table = F) 
tb <- tibble(beta = beta, mle = colMeans(mle))
# write.table(tb, file = paste0(result_loc, "/bias_arch.txt"), row.names = F, col.names = F)
# tb <- fread(paste0(result_loc, "bias_arch.txt"), data.table = F)
a <- lm(mle ~ beta + 0, data = tb)$coef
# png(filename = paste0(img_loc, "bias_logistic_arch.png"), width = 5, height = 4, units = "in", res = 3000)
ind <- which(tb$beta!=0)
ggplot(tb[ind, ]) +
  geom_point(aes(x = beta, y = mle)) +
  geom_abline(slope = a, intercept = 0, color = "red")+
  xlab("Coefficients") +
  ylab("Average MLE coefficients") +
  theme_bw() +
  theme(text = element_text(size = 18))
# dev.off()

# 2 -- Normal Q-Q plot of the MLE (MVT covariates, Figure 5, Left) 
# Load the MLE from repeated sampling
mle <- fread("mle_repeated.txt", data.table = F) 
j <- "" # pick a random non-null coord
mle_single <- tibble(mle = mle[,j])

# png(filename = "/qqnorm_mvt.png",width = 5, height = 4, units = "in", res = 3000)
ggplot() + 
  geom_qq(aes(sample = mle_single)) + 
  geom_qq_line(aes(sample = mle_single), color = "red") + 
  xlab("Theoretical quantiles") + 
  ylab("Sample quantiles") + 
  theme_bw() + 
  theme(text = element_text(size = 18)) 
# dev.off()

# write.table(mle_single, file = "mle_single_mvt.txt", row.names = F, col.names = F)
# mle_single <- scan("mle_single_mvt.txt")

# 3 --  Q-Q plot of MLE versus the bootstrap MLE, both are standardized
beta_s_known <- scan("betas_known_4.txt") # Pick one random bootstrap sample 
beta_boot_known <-  fread("boot_known_4.txt", data.table = F) 
sd_boot <- apply(beta_boot_known, 2, sd)
alpha_boot <- lm(colMeans(beta_boot_known)~beta_s_known + 0, weights = 1/sd_boot^2)$coef
mle <- fread("mle_repeated.txt", data.table = F) 
j <- "" # pick a random non-null coord
bias <- mean(mle[,j]) / beta[j] # empirical bias
sd <- apply(mle, 2, sd) # empirical std.dev
qqplot <- tibble(mle = (mle[,j] - bias * beta[j])/sd[j],
                 boot = (beta_boot_known[,j] - alpha_boot * beta_s_known[j]) / sd_boot[j])

# png(filename =  "qqboot_mvt.png",width = 5, height = 4, units = "in", res = 300)
ggplot(qqplot) + 
  geom_point(aes(x = sort(mle), y = sort(boot))) + 
  geom_abline(intercept = 0, slope = 1, color = "red") + 
  xlab("Standardized MLE") + 
  ylab("Standardized bootstrap MLE") + 
  theme_bw() + 
  theme(text = element_text(size = 18)) 
# dev.off()
# write.table(qqplot, file = "boot_single_mvt.txt", row.names = F, col.names = F)
# qqplot <- fread("boot_single_mvt.txt", data.table = F)

# 4 -- Normal Q-Q plot for i.i.d. large  in Appendix (n = 4000 and p=400)
# Load the MLE from repeated sampling
mle <- fread("mle_repeated.txt", data.table = F) 
j <- "" # randomly pick a non-null variable 
# Normal Q-Q plot
mle_single <- tibble(mle = mle[,j])

png(filename = paste0(img_loc, "qqnorm_iid_large.png"),width = 5, height = 4, units = "in", res = 3000)
ggplot() + 
  geom_qq(aes(sample = mle_single)) + 
  geom_qq_line(aes(sample = mle_single), color = "red") + 
  xlab("Theoretical quantiles") + 
  ylab("Sample quantiles") + 
  theme_bw() + 
  theme(text = element_text(size = 18)) 
dev.off()

# 5 -- Normal Q-Q plot when the sample size is small (Figure 6, Left)
n <- 400; p <- 40
rpareto <- function(n, df){
  mean <- df / (df - 1)
  sd <- sqrt(df / (df - 2) / (df - 1)^2) # standard deviation
  u <- runif(n, 0, 1)
  ((1/(1-u))^(1/df) - mean) / sd
}
Sigma <- matrix(0, p, p)
diag(Sigma) <- 1
sample_x <- function(){
  X <- matrix(0, n, p)
  X <- matrix(rpareto(n*p, 5), nrow =n)
  X <- scale(X) / sqrt(p)
  X
}
family <- get("binomial")(link = "logit")
sample_y <- function(X, family, beta){
  eta <- X %*% beta
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

beta <- scan("beta_small.txt")
gamma <- sqrt(t(beta) %*% (Sigma %*% beta) / p) 
B <- 10000 
mle <- matrix(0, B, p)
sd <- matrix(0, B, p)

for(b in 1:B){
  X <- sample_x()
  Y <- sample_y(X, family, beta)
  
  fit <- glm(Y ~ X + 0, family = family)
  mle[b, ] <- fit$coef
  sd[b, ] <- summary(fit)$coef[,2]
  
  if(b %% 50 == 0) cat(b, ",")
  if(b %% 1000 == 0) cat("\n")
}

j <- 7 
mle_single <- mle[,j]
# png(filename = "qqnorm_small.png",width = 5, height = 4, units = "in", res = 3000)
ggplot() + 
  geom_qq(aes(sample = mle_single)) + 
  xlab("Theoretical quantiles") + 
  ylab("Sample quantiles") + 
  theme_bw() + 
  theme(text = element_text(size = 18)) 
# dev.off()

# 6 - Q-Q plot of MLE versus the bootstrap MLE, both are standardized
X <- sample_x()
Y <- sample_y(X, family, beta)
fit <- glm(Y ~ X + 0, family = family)
beta_hat <- fit$coef
s <- gamma / sqrt(t(beta_hat) %*% (Sigma %*% beta_hat) / p)
beta_s <- s[1,1] * beta_hat

# Generate bootstrap samples 
boot_known <- matrix(0, B, p)
for(b in 1:B){
  yb <- sample_y(X, family, beta_s)
  fitb <- glm(yb ~ X + 0, family = family)
  boot_known[b, ] <- fitb$coef
  
  if(b %% 50 == 0) cat(b, ",")
  if(b %% 500 == 0) cat("\n")
}

sd_boot <- apply(boot_known, 2, sd)
alpha_boot <- lm(colMeans(boot_known)~beta_s + 0, weights = 1/sd_boot^2)$coef
bias <- mean(mle[,j]) / beta[j] # empirical bias
sd <- apply(mle, 2, sd) # empirical std.dev
qqplot <- tibble(mle = (mle[,j] - bias * beta[j])/sd[j],
                 boot = (boot_known[,j] - alpha_boot * beta_s[j]) / sd_boot[j])

# png(filename =  "qqboot_small.png",width = 5, height = 4, units = "in", res = 300)
ggplot(qqplot) + 
  geom_point(aes(x = sort(mle), y = sort(boot))) + 
  geom_abline(intercept = 0, slope = 1, color = "red") + 
  xlab("Standardized MLE") + 
  ylab("Standardized bootstrap MLE") + 
  theme_bw() + 
  theme(text = element_text(size = 18)) 
# dev.off()
