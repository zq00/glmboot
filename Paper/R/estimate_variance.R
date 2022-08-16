#' Estimates variance \eqn{Var(X_{\text{new}} \hat{\beta}(s))} for a shrinkage factor `s`
#' 
#' For a shrinkage factor `s`, let \eqn{\beta(s) = s\times \hat{\beta}}, where \eqn{\hat{\beta}}
#' is the MLE. Fix the coefficient at \eqn{\beta(s)} and resample response `Y`, and 
#' compute the MLE in this bootstrap sample. Estimate \eqn{Var(X_{\text{new}} \hat{\beta}(s))}
#' at this \eqn{\beta(s)}. Repeat this process `b_var` times and return all the estimates.
#' 
#' @param X Covariate matrix of size n by p
#' @param y A vector of response variable
#' @param beta_hat A vector of the MLE coefficients
#' @param s A shrinkage factor
#' @param family A `family` object, which includes a GLM family and a link function 
#' @param b_var The number of bootstrap samples to estimate \eqn{Var(X_{\text{new}} \hat{\beta}(s))} 
#' @returns
#' Returns -1 if GLM function reports error in more than 50% of times.
#' Otherwise, return a vector of the estimated \eqn{Var(X_{\text{new}} \hat{\beta}(s))} in `b_var` repetitions. 
estimate_variance <- function(x, y, beta_hat, s, family, b_var){
  beta_s <- beta_hat * s
  sd_hat <- numeric(b_var)
  
  b <- 0; nerror <- 0
  while(b < b_var){
    beta_hat_new <- boot_estimate_eta(x, beta_s, family)
    if(is.numeric(beta_hat_new) && beta_hat_new == -1) {# handling GLM error 
      nerror <- nerror + 1; 
      if(nerror > 0.5 * b_var) {return(-1)}
      next;
    }
    b <- b + 1; 
    sd_hat[b] <- beta_hat_new$eta_hat
  }
  
  return(list(sd_hat  = sd_hat))
}
