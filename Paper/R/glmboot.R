#' Apply resized bootstrap to a fitted GLM model 
#' 
#' The resized bootstrap method estimates the bias and standard deviation of the MLE for a generalized linear model (GLM) 
#' by applying a parametric bootstrap at a \emph{resized} coefficient \eqn{\beta_{\star}} such that
#' \eqn{gamma^2 \approx var(X^\top \beta_{\star})}. Here, \eqn{\gamma} is the signal strength parameter, i.e., \eqn{gamma^2 = var(X^\top \beta)}.
#' @param glm_fit A `glm` object returned by the `glm` function. it should contain an `x` and `y`, i.e., you should set `glm(.., x = T, y = T)` when running `glm` function
#' @param s_interval Increment of the sequence of shrinkage factors. E.g. when `s=0.1`, the sequence of shrinkage factors are {0, 0.1, 0.2, ..., 1}. By default, `s = 0.02`. 
#' @param b_var The number of bootstrap samples at each `s` to estimate signal strength \eqn{gamma}. By default, `b_var = 5`.
#' @param b_boot The number of bootstrap samples to estimate the bias and std.dev of the MLE.
#' @param robust_est If true, use robust mean and std.dev estimator for the bias and std.dev of the MLE
#' @param filename If provided, save the plot of \eqn{eta} versus `s` when estimating \eqn{gamma} 
#' @param verbose Print progress if True 
#' @returns 
#' \itemize{
#' \item beta_s The resized coefficient used in the parametric bootstrap.
#' \item  gamma_hat The estimated signal strength parameter
#' \item  The average of the MLE 
#' \item alpha The estimated bias of the MLE
#' \item sd The estimated std.dev of the MLE
#' \item boot_sample Bootstrap samples 
#' }
#' NOTES
# Deleted the SLOE argument -- always use SLOE to estimate the signal strength
# Deleted beta_b from the output
glm_boot <- function(glm_fit, s_interval = 0.02, b_var = 5, b_boot = 100, robust_est = FALSE, verbose = TRUE, filename = NA){
  # 1. Extract X, Y and MLE from the glm 
  family <- glm_fit$family # family and link 
  family$simulate_fun <- get_simulate_fun(family) # a function to simulate Y from the linear predictor
  X <- glm_fit$x; Y <- glm_fit$y # Covariate matrix X contains a first column of 1s if the model contains an intercept
  n <- nrow(X); p <- ncol(X) 
  beta_hat <- glm_fit$coef
  
  # 2. Estimate the signal strength gamma
  # estimates eta_tilde 
  sd_obs <- estimate_eta(X, Y, beta_hat, family); if(verbose) cat("Observed eta = ", sd_obs, "\n");if(sd_obs > 1000) stop("Error!") 
  
  # for a sequence of shrinkage factor, estimate var(X beta_hat) and alpha_hat
  s_seq <- seq(0, 1, by = s_interval); ns <- length(s_seq)
  sd_hat <- matrix(0, ns, b_var); i <- 0
  for(s in s_seq){
    new_val <- estimate_variance(X, Y, beta_hat, s, family, b_var)
    # stop if cannot fit glm
    if((is.numeric(new_val) && new_val == -1) || mean(new_val$sd_hat) > 1.5 * sd_obs) break;
    i <- i+1; sd_hat[i, ] <- new_val$sd_hat
    if(verbose){if(i %% 2 == 0){cat(s_seq[i], "\t Estimated std is ", mean(sd_hat[i, ]),"\n") }}
  }
  if(i == 1) {s <- 0;sol <- list(gamma_hat = 0);}else{
    s_seq <- s_seq[1:i]; sd_hat <- sd_hat[1:i, ]
    # find solutions 
    sol <- estimate_gamma(X, s_seq, sd_hat, beta_hat, sd_obs, verbose = verbose, filename = filename)
    s_hat <- sol$s_hat
    beta_s <- beta_hat * s_hat; 
    if(verbose){cat("Estimated gamma is", sol$gamma_hat , "\n")}}
  
  # 3. Using bootstrap to estimate the bias and variance
  mle_boot <- boot(X, beta_s, family, b_boot, verbose)
  
  # 4. estimate alpha and sigma
  if(robust_est){
    sd_boot <- apply(mle_boot, 1, function(t) Qn(t, constant = 2.2219, finite.corr = F))
  }else{
    sd_boot <- apply(mle_boot, 1, sd)
  }
  
  if(s_hat == 0){ 
    alpha_boot <- 1
    cat("No signal! \n")
  }else{ 
    if(robust_est){
      mean_ap <- apply(mle_boot, 1, median)
    }else{
      mean_ap <- rowMeans(mle_boot)
    }
    alpha_boot <- lm(mean_ap ~ beta_s + 0, weights = 1/sd_boot^2)$coef
  }
  
  return(list(beta_s = beta_s, gamma_hat = sol$gamma_hat, alpha = alpha_boot, sd = sd_boot, boot_sample = mle_boot))
}
