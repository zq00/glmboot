#' Estimates \eqn{eta} in one parametric bootstrap sample
#' 
#' Generates new responses `Y` using a given coefficient \eqn{\beta} and estimates 
#' \eqn{\eta^2 = Var(X_\text{new}^\top \hat{\beta})}, where \eqn{\hat{\beta}} is the MLE
#'
#' @param X Covariate matrix of size n by p
#' @param beta Model coefficients that is used to sample response `Y`
#' @param family A `family` object, which includes a GLM family and a link function 
#' @param sloe If True, then compute the SLOE estimator, otherwise, return only the MLE coefficients.
#' @returns 
#' \itemize{
#' \item coef The vector of MLE coefficient \eqn{\hat{\beta}}
#' \item eta_hat The estimated  \eqn{\hat{\eta}} where \eqn{\eta} is defined as \eqn{\hat{\eta}^2 = Var(X_\text{new}^\top \hat{\beta})}
#' }
#' 
boot_estimate_eta <- function(X, beta, family, sloe = T){
  Y <- family$simulate_fun(X%*%beta)
  fit <- tryCatch(error = function(e) {cat("E! "); return(-1)},
                  glm(Y ~ X + 0, family = family))
  if(length(fit) == 1 | !fit$converged){return(-1)}else{
    if(sloe) {return(list(coef = fit$coef, eta_hat = estimate_eta(X, Y, fit$coef, family)))}else{
      return(list(coef = fit$coef))
    }
  }
}

