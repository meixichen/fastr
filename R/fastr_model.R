#' @rdname fastr_fit
#' @export
fastr_model <- function(data, dt, n_factor, init=NULL, method="2step",
			lam=NULL, nu=15, woodbury=TRUE, left_out_neuron=0,
			integrate_random=TRUE,
			...){
  if (length(dim(data))==3){ # multiple neurons
    n_cell <- dim(data)[1]
    n_bin <- dim(data)[2]
    n_trial <- dim(data)[3]
  } else if (length(dim(data))==2){ # single neuron
    n_cell <- 1
    n_bin <- dim(data)[1]
    n_trial <- dim(data)[2]
    woodbury <- FALSE # DO NOT use the woodbury for single neuron model as it will crash R
    n_factor <- 1 # override whatever n_factor that is specified
    method <- "joint" # canNOT use two-step method for single neuron model
  } else{
    stop("Check the dimension of data, which must be either `n_cell x n_bin x n_trial`, or
         `n_bin x n_trial` for a single neuron.")
  }

  if (n_factor < 1) stop("n_factor must be an integer >= 1.")

  # Get k and a MLE
  all_mle <- get_ig_mle(data, dt)
  log_k_hat <- all_mle$log_k
  log_a_hat <- all_mle$log_a

  if (method == "2step"){
    if (is.null(init)){
      Lt <- rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2)
    }
    else{
      Lt <- init$Lt
    }
    if (is.null(init$x)){
      x <- prop_paths(data, dt, log_k_hat, log_a_hat)
    }
    else{
      x <- init$x
    }
    init_param <- list(log_k = log_k_hat,
                       log_a = log_a_hat,
                       Lt = Lt,
                       x = x)
    map <- list(log_k=rep(factor(NA), n_cell),
                log_a=rep(factor(NA), n_cell))
  }
  else if (method == "joint"){
    if (is.null(init)){
      init_param <- list(log_k = log_k_hat,
                         log_a = log_a_hat,
                         Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                         x = prop_paths(data, dt, log_k_hat, log_a_hat))
    }
    else{
      if (is.null(init$x)){
        x <- prop_paths(data, dt, init$log_k, init$log_a)
      }
      else{
        x <- init$x
      }
      init_param <- list(log_k = init$log_k,
                         log_a = init$log_a,
                         Lt = init$Lt,
                         x = x)
    }
    map <- NULL
  }
  else{
    stop("method must be one of '2step' or 'joint'.")
  }


  if (n_cell == 1){
    model_choice <- "single_model"
    data <- list(model=model_choice, dt=dt, Y=data, nu=as.double(nu))
    init_param <- init_param[-which(names(init_param)=="Lt")] # rm Lt from param list
  } else{
    model_choice <- ifelse(woodbury, "factor_model_eff", "factor_model")
    if (is.null(lam)) lam <- ifelse(n_cell<10, 1, 0.5)
    data <- list(model=model_choice, n_factor=n_factor, dt=dt, Y=data,
                 lam=as.double(lam), nu=as.double(nu), held_out_cell=left_out_neuron)
  }

  if (integrate_random){
    random <- "x"
  } else {
    random <- NULL
  }
  out <- list(data=data, init_param=init_param, map=map, random=random,
	      n_cell=n_cell, n_bin=n_bin, n_trial=n_trial, all_mle=all_mle)
  out
}

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
    # Get SEs for the IG parameterization with gamma and lambda
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

#' Propose an initial path (x) given parameters k and alpha
#'
#' @param Y A n x r array or q x n x r array of data where q is num of neurons, n is num of bins, r is num of trials
#' @param dt Length of time bin
#' @param log_k Length q vector of log threshold parameters
#' @param log_a Length q vector of log drift parameters
#' @noRd
prop_paths <- function(Y, dt, log_k, log_a){
  k <- exp(log_k)
  alpha <- exp(log_a)
  dimY <- dim(Y)
  # If only 1 neuron
  if (length(dimY) == 2){
    if (length(log_k)!=1 | length(log_a)!=1) stop("Assume one neuron. Check dimensions of arguemnts.        ")
    n <- dimY[1]
    r <- dimY[2]
    output <- array(NA, dim = c(n, r))
    for (i in 1:r){
      y <- Y[,i]
      U <- (cumsum(y) + 1)*k # upper bounds for neuron path
      L <- cumsum(y)*k # lower bounds for neuron path
      prop_x <- rep(0, n)
      for (j in 1:n){
        if (j > 1){
          last_x <- prop_x[j-1]
        }
        else{
          last_x <- 0
        }
        prop_x[j] <- truncnorm::rtruncnorm(1, a=L[j], b=U[j], mean=last_x+alpha*dt, sd=sqrt(dt))
      }
      output[ ,i] <- prop_x
    }
    return(output)
  }
  # If there are more than 1 neuron
  q <- dimY[1] # num of neurons
  n <- dimY[2] # num of bins
  r <- dimY[3] # num of trials
  output <- array(NA, dim = c(q, n, r))
  for (i in 1:r){
    for (w in 1:q){
      y <- Y[w,,i]
      U <- (cumsum(y) + 1)*k[w] # upper bounds for neuron path
      L <- cumsum(y)*k[w] # lower bound
      prop_x <- rep(0, n)
      for (j in 1:n){
        if (j > 1){
          last_x <- prop_x[j-1]
        }
        else{
          last_x <- 0
        }
        prop_x[j] <- truncnorm::rtruncnorm(1, a=L[j], b=U[j], mean=last_x+alpha[w]*dt, sd=sqrt(dt))
      }
      output[w, ,i] <- prop_x
    }
  }
  output
}

