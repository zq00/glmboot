#' Generates a function to simulate Y given linear predictor
#' 
#' For a given family and link function, `get_simulate_fun` returns a 
#' function that samples Y by the linear predictor \eqn{X^\top \beta}. 
#' 
#' Currently, the function supports two families: `family = binomial`
#' or `family = poisson`. For example, if `family = binomial` and `link = 'logit' `, then
#' sample Y from 0 or 1 with \eqn{P(Y=1|\eta=X^\top \beta) = 1/(1+\exp(-eta))}.
#' @param family A `family` object, which includes a GLM family and a link function 
#' @returns A function which takes a linear predictor \eqn{\eta = X^\top \beta} as the input and returns 
#' one simulated `Y` given the linear predictor. The input can also be a vector, and in that case the function
#' returns a vector of simulated `Y`.
get_simulate_fun <- function(family){
  if(family$family == "binomial"){
    simulate_fun <- function(t) rbinom(length(t), 1, family$linkinv(t))
  }
  if(family$family == "poisson"){
    simulate_fun <- function(t) rpois(length(t), family$linkinv(t))
  }
  return(simulate_fun)
}
