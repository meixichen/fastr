#' Inverse logit function
#' 
#' @param x A real number.
#' @return Inverse logit of x, a number between 0 and 1.
#' @export

invlogit <- function(x){
  1/(1+exp(-x))
}