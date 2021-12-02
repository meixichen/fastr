#' Normalize the rows of L and Psi
#' 
#' @param L Loading matrix
#' @param psi Vector of diagonal elements of the error term matrix
#' @param prec Precision of the output
#' @return A list containing the normalized L matrix and the diagonal vector of Psi, such that the norm of each row of L plus the diagonal of Psi equals 1
#' @export

norm_mat <- function(L, psi, prec){
  n_cell <- nrow(L)
  for (i in 1:n_cell){
    norm2 <- sum(L[i,]^2) + psi[i]
    L[i,] <- L[i,]/sqrt(norm2)
    psi[i] <- psi[i]/norm2
  }
  return(list(L=round(L, prec), psi=round(psi, prec)))
}