#' Transform a real number to a correlation coefficient
#' @param x A real number.
#' @return A number between -1 and 1.
#' @details This function is the inverse function of `cor2real()`.
#' @export

real2cor <- function(x){
  invlogit(x)*2-1
}