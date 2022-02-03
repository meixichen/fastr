#' Propose an initial path (x) given parameters k and alpha
#' 
#' @param Y A n x r array or q x n x r array of data where q is num of neurons, n is num of bins, r is num of trials
#' @param dt Length of time bin
#' @param log_k Length q vector of log threshold parameters
#' @param log_a Length q vector of log drift parameters
#' @export

prop_paths <- function(Y, dt, log_k, log_a){
  k <- exp(log_k)
  alpha <- exp(log_a)
  dimY <- dim(Y)
  # If only 1 neuron
  if (length(dimY) == 2){
    if (length(log_k)!=1 | length(log_a)!=1) stop("Assume one neuron. Check dimensions of arguemnts.") 
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
