#' Compute bootstrap MLE 
#' 
#' Reports error if bootstrap MLE does not exist in more than 20% of the bootstrap samples
#' 
#' @param X Covariate matrix of size n by p
#' @param beta Model coefficients that is used to sample response `Y`
#' @param family A `family` object, which includes a GLM family and a link function 
#' @param b_boot Number of bootstrap samples
#' @param verbose Print progress if True 
#' 
#' @returns A matrix of size p by b_boot containing the bootstrap MLE. Returns error if bootstrap MLE does not exist more than 20% times 
boot <- function(X, beta, family, b_boot, verbose){
  p <- ncol(X)
  b <- 0; nerror <- 0
  mle_boot <- matrix(0, p, b_boot)
  if(verbose) cat("\nNumber of bootstrap samples: \n")
  while(b < b_boot){
    if(nerror / b_boot > 0.2){error("Bootstrap MLE does not exist for shrinked coefficients!")}
    mle_boot_new <- boot_estimate_eta(X, beta, family, sloe = F)
    if(is.numeric(mle_boot_new) && mle_boot_new == -1) {nerror <- nerror + 1; next;}
    b <- b + 1; mle_boot[ ,b] <- mle_boot_new$coef
    
    if(verbose){ if(b %% (b_boot/100) == 0) cat(b, " "); if(b %% (b_boot/10) == 0) cat("\n")}
  }
  
  return(mle_boot)
}

