#' Sample from a multivariate normal with sparse precision matrix.
#'
#' @param n Number of random draws.
#' @param mean Mean vector.
#' @param prec Sparse precision matrix, i.e., inheriting from 
#' [Matrix::sparseMatrix] or its Cholesky factor, i.e., 
#' inheriting from [Matrix::CHMfactor].
#'
#' @return A matrix with `n` rows, each of which is a draw from the 
#' corresponding normal distribution.
#'
#' @details If the matrix is provided in precision form, it is converted to 
#' Cholesky form using `Matrix::Cholesky(prec, super = TRUE)`.  Once it is of 
#' form [Matrix::CHMfactor], this function is essentially copied from local 
#' function `rmvnorm()` in function `MC()` defined in [TMB::MakeADFun()].
#' @export
rmvn_prec <- function(n, mean, prec) {
  d <- ncol(prec) # number of mvn dimensions
  if(!is(prec, "CHMfactor")) {
    prec <- Matrix::Cholesky(prec, super = TRUE)
  }
  u <- matrix(rnorm(d*n),d,n)
  u <- Matrix::solve(prec,u,system="Lt")
  u <- Matrix::solve(prec,u,system="Pt")
  u <- t(as(u, "matrix") + mean)
}