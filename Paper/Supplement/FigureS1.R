# Results of MVT theory
# From the data extract a data table where each row has 
# (kappa, gamma, alpha, sigma, lambda)
library(ggplot2)
library(glmhd)
dir <- "/theory/"
figdir <- "/fig/"

# extract filenames
filenames <- list.files(dir)

# function to extract the theoretical parameters
get_param <- function(name){
  # extract kappa and gamma values
  vals <- strsplit(strsplit(name, ".txt")[[1]][1], "_")[[1]]
  if(length(vals) == 2){
    kappa <- as.numeric(vals[1])
    gamma <- as.numeric(vals[2])
    nu <- 8
  }else{
    kappa <- as.numeric(vals[2])
    gamma <- as.numeric(vals[3])
    nu <- 15
  } 
    
  # read file
  data <- read.table(paste0(dir, name), header = T)
  # the *first* time score becomes positive 
  if(sum(which(data$score > 0)) == 0){
    cat(name, ": increase alpha value!")
    return(c(kappa,gamma, 0,0,0,0))
  }
  ind <- min(which(data$score > 0)) 
  # linearly interpolate the a, s, l values
  as <- (data$a[ind - 1] * data$score[ind] - data$a[ind] * data$score[ind - 1]) / (data$score[ind] - data$score[ind-1])
  ss <- (data$s[ind - 1] * data$score[ind] - data$s[ind] * data$score[ind - 1]) / (data$score[ind] - data$score[ind-1])
  ls <- (data$l[ind - 1] * data$score[ind] - data$l[ind] * data$score[ind - 1]) / (data$score[ind] - data$score[ind-1])
  # plot the score function 
  plot(data$a, data$score, pch = 16, xlab = "alpha",ylab = "score")
  abline(h = 0)
  
  return(c(nu,kappa, gamma, as, ss, ls))
}

params <- matrix(NA, nrow = length(filenames), ncol = 6)
for(i in 1:length(filenames)){
  params[i, ] <- get_param(filenames[i])
}
colnames(params) <- c("nu","kappa", "gamma", "alpha", "sigma", "lambda")
params <- as.data.frame(params)
# remove the ones that do not work
params <- params[params$kappa == 0.1, ]
params$nu <- as.factor(params$nu)

# Empirical results
dir_empirical <- "/theory-empirical/"

# extract filenames
filenames <- list.files(dir_empirical)

# function to extract empirical parameters
get_empirical <- function(name){
  # extract kappa and gamma values
  vals <- strsplit(strsplit(name, ".txt")[[1]][1], "_")[[1]]
  nu <- as.numeric(vals[1])
  kappa <- as.numeric(vals[2])
  gamma <- as.numeric(vals[3])
  
  # read file
  data <- read.table(paste0(dir_empirical, name), header = T)
  return(c(nu, kappa, gamma, colMeans(data))) 
}

params_empirical <- matrix(NA, nrow = length(filenames), ncol = 6)
for(i in 1:length(filenames)){
  params_empirical[i, ] <- get_empirical(filenames[i])
}
params_empirical <- as.data.frame(params_empirical)
colnames(params_empirical) <- c("nu","kappa", "gamma", "alpha_e", "sigma_e", "lambda_e")
params_empirical$nu <- as.factor(params_empirical$nu)
params_empirical <- params_empirical[params_empirical$kappa  == 0.1, ] 

# Compute the theoretical parameters when the covariates are MVN
gamma <- seq(0, 5, by = 0.1)
params_mvn <- matrix(0, nrow = length(gamma), ncol = 4)
for(i in 1:length(gamma)){
  cat(i, ",")
  new_param <- find_param(kappa = 0.1, gamma = gamma[i], beta0 = 0, intercept = F,verbose = FALSE)
  params_mvn[i, ] <- c(gamma[i], new_param)
}
params_mvn <- as.data.frame(params_mvn)
params_mvn <- params_mvn[-1, ]
colnames(params_mvn) <- c("gamma", "alpha_s", "lambda_s", "sigma_s")
params_mvn$sigma_s <- params_mvn$sigma_s/sqrt(0.1)

# Figure 1 alpha values
g_alpha <-  ggplot(params) + 
  geom_line(aes(x = gamma, y = alpha, group=nu, color = nu),  size = 1) + 
  geom_point(aes(x= gamma, y = alpha_e, color = nu), data = params_empirical, size = 1) +
  geom_line(aes(x = gamma, y = alpha_s), color = "black", size = 0.8, data = params_mvn) + 
  xlab(expression(gamma)) + 
  ylab(expression(alpha[star])) + 
  theme_bw() + 
  theme(text = element_text(size = 15))+
   scale_color_manual(name="DOF",
                      breaks = c(8, 15),
                      values = c("#E1BE6A", "#40B0A6"),
                      labels=c("nu = 8", "nu = 15")) 


ggsave(plot = g_alpha, filename = paste0(figdir, "alpha.png"),
       width = 10, height = 6, units = "cm")
g_sigma <- ggplot(params) + 
  geom_line(aes(x = gamma, y = sigma, color = nu , group = nu), size = 1) + 
  geom_point(aes(x= gamma, y = sigma_e, color = nu), data = params_empirical, size = 1) +
  geom_line(aes(x = gamma, y = sigma_s), color = "black", size = 0.8, data = params_mvn) + 
  xlab(expression(gamma)) + 
  ylab(expression(sigma[star])) + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  scale_color_manual(name="DOF",
                     breaks = c(8, 15),
                    values = c("#E1BE6A", "#40B0A6"),
                    labels=c("nu = 8", "nu = 15")) 

ggsave(plot = g_sigma, filename = paste0(figdir, "sigma.png"),
       width = 10, height = 6, units = "cm")

ggplot(params) + 
  geom_line(aes(x = gamma, y = lambda )) + 
  geom_point(aes(x= gamma, y = lambda_e), data = params_empirical) +
  xlab(expression(gamma)) + 
  ylab(expression(lambda[star])) + 
  theme_bw() + 
  theme(text = element_text(size = 15)) +
  theme(text = element_text(size = 15)) + 
  scale_color_manual(name="DOF",
                     breaks = c(8, 15),
                     values = c("#E1BE6A", "#40B0A6"),
                     labels=c("nu = 8", "nu = 15")) 
  



