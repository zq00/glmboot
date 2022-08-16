#' Computes the derivatives of the negative log-likelihood
#' 
#' \eqn{g(t) = f'_y(t)} and \eqn{gprime(t) = f''_y(t)}
#' \eqn{f_y(t)} is the negative log-likelihood when the linear predictor is `t` and 
#' the response is `y`. In case of a logistic regression,
#' \eqn{f_y(t) = \log(1+\exp(-yt)) } for \eqn{y\in\{\pm 1\}}.
#' 
#' The function currently supports binary regression with `logit`, `probit` and `cloglog` links, 
#' or Poisson regression with `log` link.
#' 
#' @param family A GLM family
#' @returns `g` and `gprime` are the derivatives and the second derivative of \eqn{f_y(t)} respectively
getg <- function(family){
  if(family$family == "poisson"){
    if(family$link == "log"){
      g <- function(y, t) -y + exp(t)
      gprime <- function(y, t) exp(t)
    }
  }
  
  if(family$family == "binomial"){
    if(family$link == "logit"){
      g <- function(y, t) -y/(1 + exp(y * t))
      gprime <- function(y, t) 1/(1+exp(t))/(1+exp(-t))
    }else if(family$link == "probit"){
      g <- function(y, t) -y * dnorm(y*t)/pnorm(y * t)
      gprime <- function(y, t) {
        phiprime <- -(y * t)*exp(-t^2 / 2)/sqrt(2*pi)
        dnorm(y*t)^2 / pnorm(y * t)^2 - phiprime / pnorm(y*t)
      }
    }else if(family$link == "cloglog"){
      g <- function(y, t){
        n <- length(y)
        val <- numeric(n)
        for(i in 1:n){
          val[i] <-ifelse(y[i] > 0, - exp(t[i]) / (exp(exp(t[i])) - 1), exp(t[i]))
        }
        val
      }
      gprime <- function(y, t){
        n <- length(y)
        val <- numeric(n)
        for(i in 1:n){
          val[i] <- ifelse(y[i] > 0, if(t[i] < 5) {return((exp(t[i]) * (-exp(exp(t[i])) + exp(t[i] + exp(t[i])) + 1)) / (exp(exp(t[i])) - 1)^2)}else{return(0)}, exp(t[i])) 
        }
        val
      }
    }
  }
  return(list(g = g, gprime = gprime))
}
