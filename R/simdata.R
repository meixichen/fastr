#' Simulate spike trains given parameters
#' 
#' @param dt Time bin length.
#' @param n_bin Number of time bins.
#' @param n_trial Number of trials.
#' @param alpha Vector of drift parameters.
#' @param k Vector of threshold parameters.
#' @param rho Vector of correlation parameters of length q(q-1)/2 where q is the length of alpha. 
#' If this is provided, the pairwise correlation model is used for simulation.
#' @param L Loading matrix of dimension q x d, where d is the number of latent factors. 
#' If this is provided, the factor model is used.
#' @export

simdata <- function(dt, n_bin, n_trial, alpha, k, rho, L){
  n_cell <- length(alpha)
  if (n_cell != length(k)){
    stop("Lengths of alpha and k must match.")
  }
  
  # Simulate 1 neuron
  if (n_cell == 1){
    dx <- array(0, dim=c(n_bin, n_trial))
    Y <- array(0, dim=c(n_bin, n_trial))
    for (u in 1:n_trial){
      Nt <- 1
      for (t in 1:n_bin){
      	dx[t,u] <- rnorm(1, mean=alpha*dt, sd=sqrt(dt))
      	p_jump <- ifelse(sum(dx[1:t,u])>= Nt*k, 1, 0)
      	Y[t,u] <- rbinom(1, 1, p_jump)
      	Nt <- Nt + Y[t,u]
      }
    }
    x <- array(0, dim = c(n_bin, n_trial))
    for (u in 1:n_trial){
      x[,u] <- cumsum(dx[,u])
    }
    
    return(list(Y=Y, x=x))
  }
  
  # Simulate more than 1 neuron
  if (missing(rho)){
    if (missing(L)){
      stop("Must provide either rho or L")
    }
    cov <- L%*%t(L)
    diag(cov) <- rep(1, n_cell)
    cov <- cov*dt
  }
  else{
    if (length(rho) != n_cell*(n_cell-1)/2){
      stop("Length of rho must be n_cell*(n_cell-1)/2")
    }
    cov <- matrix(dt, nrow=n_cell, ncol=n_cell)
    for (i in 1:(n_cell-1)){
      for (j in (i+1):n_cell){
        cov[i,j] <- rho[n_cell*(i-1) - i*(i-1)/2 + j-i]*dt # fill the cov matrix along the column (from top to bottom, and then from left to right)
        cov[j,i] <- cov[i,j]
      }
    }
  }

  dx <- array(0, dim = c(n_cell,n_bin, n_trial))
  Y <- array(0, dim = c(n_cell,n_bin, n_trial))
  
  for (u in 1:n_trial){
    Nt <- rep(1, n_cell)
    for (t in 1:n_bin){
      dx[,t,u] <- mvtnorm::rmvnorm(1, mean=alpha*dt, sigma=cov)
      for (i in 1:n_cell){
        p_jump <- ifelse(sum(dx[i,1:t,u])>= Nt[i]*k[i], 1, 0)
        Y[i,t,u] <- rbinom(1, 1, p_jump)
        Nt[i] <- Nt[i] + Y[i,t,u]
      }
    }
  }
  x <- array(0, dim = c(n_cell,n_bin, n_trial))
  for (u in 1:n_trial){
    x[,,u] <- t(apply(dx[,,u], 1, cumsum))
  }
  
  return(list(Y=Y, x=x))
}
