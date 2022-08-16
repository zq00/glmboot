#' Estimate the signal strength \eqn{\gamma}
#' 
#' For a sequence of shrinkage factor `s`, let \eqn{\beta(s) = s\times \hat{\beta}}, correspondingly, 
#' \eqn{\gamma^2(s) = Var(X \beta(s))}, where \eqn{X} is the observed covariate matrix. Use parametric bootstrap
#' to generate new responses with coefficients equal to \eqn{\beta(s)} and compute the MLE \eqn{\hat{\beta}(s)}.
#' Then, estimate \eqn{\eta^2(s) = Var(X_\text{new}^\top \hat{\beta}(s))}. Estimate \eqn{\gamma} such that
#' \eqn{\eta(s)} equals what is observed in the sample (where `Y` is generated from the unknown model coefficient \eqn{\beta}).
#' 
#' @param x Covariate matrix of size n by p
#' @param s_seq A sequence of shrinkage factors
#' @param sd_hat Estimated \eqn{\hat{\eta}(s)}. `sd_hat` is a matrix where the number of rows is the length of `s_seq`
#' and the number of columns is the number of repetitions at each `s`.
#' @param beta_hat A vector of the MLE coefficients
#' @param sd_obs The observed \eqn{\hat{\eta}} from the observed data
#' @param verbose Print progress if True 
#' @param filename If provided, save the plot of \eqn{eta} versus `s` 
#' 
#' @returns
#' \itemize{
#' \item s_hat A numeric value of the estimated \eqn{\hat{s}}, i.e., the estimated \eqn{\hat{\gamma}} is
#' \eqn{\gamma^2(\hat{s}) = Var(X \beta(\hat{s}))}
#' \item gamma_hat A numeric value of the estimated \eqn{\gamma}.
#' }
#' 
estimate_gamma <- function(x, s_seq, sd_hat, beta_hat, sd_obs, verbose = T, filename = NULL){
  # fit a smooth loess curve 
  data <- data.frame(cbind(rep(s_seq, time = ncol(sd_hat)), as.vector(sd_hat)))
  colnames(data) <- c("s", "val")
  curve <- loess(val ~ s, data = data)
  s_new <- seq(min(s_seq), max(s_seq), by = 0.001)
  sd_new <- predict(curve, s_new) # estimated sd_hat using the smoothed loess curve
  diff <- abs(sd_new - sd_obs)
  s_sol <- s_new[which.min(diff)]
  sd <- rms(x %*% beta_hat)
  gamma <- sd * s_new; gamma_hat <- sd * s_sol # estimated gamma_hat

  if(verbose){
    # if input a file location, save the plot to the file location
    if(!is.na(filename)){png(filename, width = 500, height = 400)}
    plot(gamma, sd_new, type = "l", ylim = range(sd_new), 
         xlab = expression(gamma), ylab = expression(paste(hat(eta))))
    points( sd * data$s, data$val, pch = 16, cex = 0.5)
    abline(h = sd_obs, lty = "dotted")
    mtext(text = paste0("gamma_hat = ", round(gamma_hat, 2)), line = 0.5)
    if(!is.na(filename)) dev.off()
  }
  return(list(s_hat = s_sol, gamma_hat = gamma_hat))
}
 #' Internal function 
rms <- function(t) sqrt(sum(t^2) / length(t))
