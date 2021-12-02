#' Transform the correlation coefficient to a real number
#' @param x Correlation coefficient that lies in (-1, 1).
#' @return Transformed x, a real number.
#' @details `x` is first transformed to the scale of (0, 1), and then is applied to a logit function.
#' @export

cor2real <- function(x){ 
  logit((x+1)/2)
}