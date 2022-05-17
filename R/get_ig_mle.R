#' Calculate the MLEs for mu and lambda of an inverse Gaussian distribution
#'
#' @param Y A `q x n x r` or `n x r` array of 0/1 spike trains where `q` is the number of cells, 
#" `n` is the number of time bins, and `r` is the number of trials.
#' @param dt Scalar of length of the bins.
#' @return A list of MLEs of mu and lambda.
#' @details Let T be the vector of all ISIs. MLEs are calculated by 
#' ```
#' \hat{\mu} = \bar{T}, \hat{\lambda} = n/\sum_{i=1}^n(1/T_i - 1/\bar{T}).
#' ```
#' If this IG distribution represents the distribution of first passage time of a
#' Wiener process with threshold k, drift alpha and variance sigma^2, then mu =
#' k/alpha, lambda=k^2/sigma^2.
#' @export
get_ig_mle <- function(Y, dt){
  if (length(dim(Y)) == 3) {
    all_ISI <- apply(Y, 1, 
		     function(yy) {
		       unlist(apply(yy, 2, function(y) { diff(which(y==1)*dt) }))
                     })
    n <- sapply(all_ISI, length)
    mu_hat <- sapply(all_ISI, mean)
    lam_hat <- sapply(1:length(n), function(i) { n[i]/sum(1/all_ISI[[i]] - 1/mu_hat[i]) })
  } 
  else if (length(dim(Y)) == 2) {
    all_ISI <- unlist(apply(Y, 2, function(y) { diff(which(y==1)*dt) }))
    n <- length(all_ISI)
    mu_hat <- mean(all_ISI)
    lam_hat <- n/sum(1/all_ISI - 1/mu_hat)
  }
  else {
    stop("Check the dimension of Y.")
  }
  list(mu=mu_hat, lam=lam_hat)
}
