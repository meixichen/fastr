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

#' Pairwise correlation based on spike counts
#'
#' @param neuron1 A length `n` vector of neuron 1's binary spiking data in one trial.
#' @param neuron2 A length `n` vector of neuron 2's binary spiking data in one trial.
#' @param block_len Number of bins to combine into blocks in which the number of spikes is counted.
#' @return Correlation between two neurons' spiking counts
#' @export
pairwise_corr <- function(neuron1, neuron2, block_len){
  if (length(neuron1)!=length(neuron2)){
    stop("Data of neuron 1 and 2 must have the same length.")
  }
  n_bins <- length(neuron1)
  n_block <- ceiling(n_bins/block_len)
  start_vec <- seq(from=1, to=n_bins, by=block_len)
  data1 <- rep(0, n_block)
  data2 <- rep(0, n_block)
  for (i in 1:n_block){
    start <- start_vec[i]
    data1[i] <- sum(neuron1[start:min(i*block_len, n_bins)])
    data2[i] <- sum(neuron2[start:min(i*block_len, n_bins)])
  }
  cor(data1, data2)
}

#' Sample from a multivariate normal with sparse precision matrix.
#'
#' @param n Number of random draws.
#' @param mean Mean vector.
#' @param prec Sparse precision matrix, i.e., inheriting from
#' [Matrix::sparseMatrix] or its Cholesky factor, i.e.,
#' inheriting from "CHMfactor".
#'
#' @return A matrix with `n` rows, each of which is a draw from the
#' corresponding normal distribution.
#'
#' @details If the matrix is provided in precision form, it is converted to
#' Cholesky form using `Matrix::Cholesky(prec, super = TRUE)`.  Once it is of
#' form "CHMfactor", this function is essentially copied from local
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

