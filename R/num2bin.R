#' Convert multiple neurons' data of spike times into the binary form
#'
#' @param data A list of data frames, each being a neuron's data from multiple trials. 
#' Each column of the data frame is a trial containing spike times.
#' @param dt Time bin length
#' @param trial_len Length of a trial (in sec)
#' @return A q x n x r array, where q is the number of neurons, n is the number of time bins, 
#' and r is the number of trials.
#' @export

num2bin <- function(data, dt, trial_len){
  q <- length(data)
  r <- length(data[[1]])
  n <- ceiling(trial_len/dt)
  output <- array(0, dim = c(q, n, r))
  for (i in 1:q){
    all_trials <- data[[i]] # all data of neuron i
    for (j in 1:r){
      trial_num <- all_trials[,j] # single trial data in spike times    
      trial_bin <- rep(0, n) # single trial data in the binary form
      which_spikes <- ceiling(trial_num/dt)
      trial_bin[which_spikes[!is.na(which_spikes)]] <- 1
      output[i, ,j] <- trial_bin
    }
  }
  output
}