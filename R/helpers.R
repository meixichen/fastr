#' Transform the correlation coefficient to a real number
#' @param x Correlation coefficient that lies in (-1, 1).
#' @return Transformed x, a real number.
#' @details `x` is first transformed to the scale of (0, 1), and then is applied to a logit function.
#' @export
cor2real <- function(x){
  logit((x+1)/2)
}

#' Transform a real number to a correlation coefficient
#' @param x A real number.
#' @return A number between -1 and 1.
#' @details This function is the inverse function of `cor2real()`.
#' @export
real2cor <- function(x){
  invlogit(x)*2-1
}

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

#' Logit function
#'
#' @param x A number between 0 and 1.
#' @return Logit of x, a real number.
#' @export
logit <- function(x){
  log(x/(1-x))
}

#' Inverse logit function
#'
#' @param x A real number.
#' @return Inverse logit of x, a number between 0 and 1.
#' @export
invlogit <- function(x){
  1/(1+exp(-x))
}

#' Fill a correlation matrix with a vector of its lower triangular elements
#'
#' @param lower_tri Vector of lower triangular elements (by column) of the correlation matrix
#' @param dim Number of rows/columns of the correlation matrix
#' @return A dim x dim correlation matrix
#' @export

fill_cor_mat <- function(lower_tri, dim){
  if (length(lower_tri) != (dim*(dim-1)/2)){
    stop("Length of lower_tri must be dim*(dim-1)/2!")
  }
  cor_mat <- matrix(1, nrow=dim, ncol=dim)
  cor_mat[lower.tri(cor_mat)] <- lower_tri
  cor_mat <- t(cor_mat)
  cor_mat[lower.tri(cor_mat)] <- lower_tri
  return(cor_mat)
}

