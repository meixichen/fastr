#' Propose an initial path (x) given parameters k and alpha
#' 
#' @param Y q x n x r array of data where q is num of neurons, n is num of bins, r is num of trials
#' @param dt Length of time bin
#' @param log_k Length q vector of log threshold parameters
#' @param log_a Length q vector of log drift parameters
#' @export

prop_paths <- function(Y, dt, log_k, log_a){
  k <- exp(log_k)
  alpha <- exp(log_a)
  q <- dim(Y)[1] # num of neurons
  n <- dim(Y)[2] # num of bins
  r <- dim(Y)[3] # num of trials
  output <- array(NA, dim = c(q, n, r))
  for (i in 1:r){
    for (w in 1:q){
      y <- Y[w,,i]
      U <- (cumsum(y) + 1)*k[w] # upper bounds for neuron path
      prop_x <- rep(0, n)
      for (j in 1:n){
        if (j > 1){
          last_x <- prop_x[j-1]
        }
        else{
          last_x <- 0
        }
        prop_x[j] <- truncnorm::rtruncnorm(1, b=U[j], mean=last_x+alpha[w]*dt, sd=sqrt(dt))
      }
      output[w, ,i] <- prop_x
    }
  }
  output
}