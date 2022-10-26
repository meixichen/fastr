#' Calculate the MLEs for parameters of an inverse Gaussian distribution
#'
#' @param Y A `q x n x r` or `n x r` array of 0/1 spike trains 
#' where `q` is the number of cells, 
#' `n` is the number of time bins, and `r` is the number of trials.
#' @param dt Scalar of length of the bins.
#' @return A list of MLEs and Hessian matrix of the (log) parameters 
#' c(log_k, log_a) and the SEs of the IG parameters c(1/mu, lambda).
#' @details Let T be the vector of all ISIs. MLEs are calculated by 
#' ```
#' \hat{\mu} = \bar{T}, \hat{\lambda} = n/\sum_{i=1}^n(1/T_i - 1/\bar{T}).
#' ```
#' If this IG distribution represents the distribution of first passage time of a
#' Wiener process with threshold k, drift alpha and variance sigma^2, then mu =
#' k/alpha, lambda=k^2/sigma^2.
#' @export
get_ig_mle <- function(Y, dt){
  if (length(dim(Y)) == 3) { # More than 1 trial
    n_cell <- dim(Y)[1]
  } else if (length(dim(Y)) == 2) { # only one neuron
    n_cell <- 1
  } else {
    stop("Check the dimension of Y.")
  }
  all_ISI <- lapply(1:n_cell, bin2isi, Y=Y, dt=dt)
  
  hess_ak <- matrix(0, nrow=2*n_cell, ncol=2*n_cell)
  se_ig <- rep(NA, n_cell*2) # not storing the hessian for the ig param because it is diagonal
  mu_hat <- rep(NA, n_cell)
  gamma_hat <- rep(NA, n_cell)
  lam_hat <- rep(NA, n_cell)
  log_k_hat <- rep(NA, n_cell)
  log_a_hat <- rep(NA, n_cell)
  
  for (i in 1:n_cell){
    neuron_i_isi <- all_ISI[[i]] 
    n <- length(neuron_i_isi)
    mu_hat[i] <- mean(neuron_i_isi)
    gamma_hat[i] <- 1/mu_hat[i]
    lam_hat[i] <- n/sum(1/neuron_i_isi - 1/mu_hat[i])
    log_k_hat[i] <- 0.5*log(lam_hat[i])
    log_a_hat[i] <- log_k_hat[i] - log(mu_hat[i]) 
    # Construct Hessian for the a/k parameterization
    temp_ak <- n*exp(log_k_hat[i]+log_a_hat[i])
    hess_ak_i <- matrix(c(-temp_ak-2*n,temp_ak,temp_ak,-temp_ak), ncol=2)
    idx <- (2*(i-1)+1):(2*i)
    hess_ak[idx, idx] <- hess_ak_i
    # Get SEs for the IG parameterizatino with gamma and lambda
    se_ig[(2*i-1): (2*i)] <- sqrt(1/c(n*lam_hat[i]/gamma_hat[i], 
      			    n/(2*lam_hat[i]^2)))
  }
  names <- paste0(c("log_k", "log_a"), rep(1:n_cell, each=2))
  colnames(hess_ak) <- names
  rownames(hess_ak) <- names
  neworder <- c(seq(1,(2*n_cell), by=2), seq(2, (2*n_cell), by=2))
  hess_ak <- hess_ak[neworder, neworder] # order by k and then a

  names_ig <- paste0(c("gamma", "lambda"), rep(1:n_cell, each=2))
  names(se_ig) <- names_ig
  se_ig <- se_ig[neworder] # order by gamma and then lambda
  
  return(list(log_k=log_k_hat, log_a=log_a_hat, hess_ak=hess_ak, 
	      gamma=gamma_hat, lambda=lam_hat, se_ig=se_ig))
}
