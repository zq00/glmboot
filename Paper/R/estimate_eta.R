#' Compute the SLOE estimator
#' 
#' Estimates \eqn{eta} defined as \eqn{\eta^2 = Var(X_\text{new}^\top \hat{\beta})},
#' where \eqn{X_{\text{new}}} is a new observation and \eqn{\hat{\beta}} is the MLE.
#' 
#' @param X Covariate matrix of size n by p
#' @param y A vector of response variable
#' @param beta_hat A vector of the MLE coefficients
#' @param family A `family` object, which includes a GLM family and a link function 
#' @returns A numeric value containing the estimated \eqn{\hat{\eta}}
estimate_eta <- function(X, y, beta_hat, family){
  if(family$family == "binomial") y <- 2 * y - 1
  f <- getg(family)
  g <- f$g
  gprime <- f$gprime
  
  eta_hat <- X %*% beta_hat
  D <- as.vector(gprime(y, eta_hat))
  H <-  t(X) %*% (X * D)
  w <- diag(X %*% (solve(H) %*% t(X)))
  q <- w / (1 - D * w)
  
  eta_tilde <- eta_hat + q * g(y, eta_hat)
  sqrt(mean(eta_tilde^2) - mean(eta_tilde)^2)
}
