#' Sigmoid function with mean and scale parameters
#' 
#' @param mu Mean parameter.
#' @param sig Scale parameter.
#' @param x Input to the sigmoid function, a real number.
#' @return Sigmoid transformation of x, a number between 0 and 1.
#' @export

sigm <- function(mu, sig, x){ 
  exp(sig*(x-mu))/(1+exp(sig*(x-mu)))
}
