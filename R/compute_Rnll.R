#' Compute the negative log-likelihood of the factor model in R
#' @param data A list of a 3d data array `Y`, penalization parameter `lam`, time bin size `dt`, and num of factors `n_factor`.
#' @param param_list A list of parameters including vector `log_k`, vector `log_a`, 3d array `x` (latent paths), 
#' and vector `Lt` (unnormalized lower tri elements of loading matrix)
#' @return Scalar of negative log-likelihood
#' @export

compute_Rnll <- function(data, param_list){
  k <- exp(param_list$log_k)
  alpha <- exp(param_list$log_a)
  Lt <- param_list$Lt
  x <- param_list$x
  dim_x <- dim(x)
  n_cell <- dim_x[1]
  n_bin <- dim_x[2]
  n_trial <- dim_x[3]
  n_factor <- data$n_factor
  Y <- data$Y
  lam <- data$lam
  dt <- data$dt
  #----- fill in the covariance matrix for latent increments -----
  L <- matrix(0, nrow=n_cell, ncol=n_factor) # loading matrix
  psi <- rep(0, n_cell) # diagonal of uniqueness matrix
  indx <- 0
  for (j in 1:n_factor){ # column
    for (i in j:n_cell){ # row
      indx <- indx + 1
      L[i,j] <- Lt[indx]
    }
  }
  for (i in 1:n_cell){ # row
    norm2 <- 1
    for (j in 1:n_factor){
      norm2 <- norm2 + L[i,j]^2
    }
    psi[i] <- 1/norm2
    for (j in 1:n_factor){
      L[i,j] <- L[i,j]/sqrt(norm2)
    }
  }
  Sig <- L%*%t(L) + diag(psi)
  #----- compute nll -----------
  nll <- lam*sum(Lt^2) # l2 penalization term
  for (u in 1:n_trial){
    mu <- x[, 1, u] - alpha*dt
    nll <- nll - mvtnorm::dmvnorm(mu, sigma=Sig*dt, log=T)
    for (w in 2:n_bin){
      mu <- x[, w, u] - x[, w-1, u] - alpha*dt
      nll <- nll - mvtnorm::dmvnorm(mu, sigma=Sig*dt, log=T)
    }
    Nt <- rep(1, n_cell)
    for (j in 1:n_bin){
      for (i in 1:n_cell){
        nll <- nll - dbinom(Y[i, j, u], 1, invlogit(x[i, j, u] - Nt[i] * k[i]), log=T)
        Nt[i] <- Nt[i] + Y[i, j, u] 
      }
    }
  }
  return(nll)
}
