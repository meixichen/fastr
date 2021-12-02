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