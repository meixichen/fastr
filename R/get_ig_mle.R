#' Calculate the MLEs for mu and lambda of an inverse Gaussian distribution
#'
#' @param Y A `q x n x r` or `n x r` array of 0/1 spike trains where `q` is the number of cells, 
#' `n` is the number of time bins, and `r` is the number of trials.
#' @param dt Scalar of length of the bins.
#' @return A list of MLEs and Hessian matrix of the (log) parameters c(log_k, log_a).
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
    n_cell <- dim(Y)[1]
    hess <- matrix(0, nrow=2*n_cell, ncol=2*n_cell)
    mu_hat <- rep(NA, n_cell)
    lam_hat <- rep(NA, n_cell)
    log_k_hat <- rep(NA, n_cell)
    log_a_hat <- rep(NA, n_cell)
    for (i in 1:n_cell){
      neuron_i_isi <- unlist(apply(Y[i,,], 2, 
				   function(y) {diff(which(y==1)*dt)} ))
      n <- length(neuron_i_isi)
      mu_hat[i] <- mean(neuron_i_isi)
      lam_hat[i] <- n/sum(1/neuron_i_isi - 1/mu_hat[i])
      log_k_hat[i] <- 0.5*log(lam_hat[i])
      log_a_hat[i] <- log_k_hat[i] - log(mu_hat[i]) 
      temp <- n*exp(log_k_hat[i]+log_a_hat[i])
      hess_i <- matrix(c(-temp-2*n,temp,temp,-temp), ncol=2)
      names <- paste0(c("k", "a"), i)
      colnames(hess_i) <- names
      rownames(hess_i) <- names
      idx <- (2*(i-1)+1):(2*i)
      hess[idx, idx] <- hess_i
    }
    names <- paste0(c("log_k", "log_a"), rep(1:n_cell, each=2))
    colnames(hess) <- names
    rownames(hess) <- names
    neworder <- c(seq(1,(2*n_cell), by=2), seq(2, (2*n_cell), by=2))
    hess <- hess[neworder, neworder] # order by k and then a
  } 
  else if (length(dim(Y)) == 2) {
    all_ISI <- unlist(apply(Y, 2, 
			    function(y) { diff(which(y==1)*dt) }))
    n <- length(all_ISI)
    mu_hat <- mean(all_ISI)
    lam_hat <- n/sum(1/all_ISI - 1/mu_hat)
    log_k_hat <- 0.5*log(lam_hat)
    log_a_hat <- log_k_hat - log(mu_hat) 
    temp <- n*exp(log_k_hat+log_a_hat)
    hess <- matrix(c(-temp-2*n,temp,temp,-temp), ncol=2)
  }
  else {
    stop("Check the dimension of Y.")
  }
  
  return(list(log_k=log_k_hat, log_a=log_a_hat, hess=hess))
}
