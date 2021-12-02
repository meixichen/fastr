#' Logit function
#' 
#' @param x A number between 0 and 1.
#' @return Logit of x, a real number.
#' @export

logit <- function(x){
  log(x/(1-x))
}