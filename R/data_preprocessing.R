#' Convert a binary array of spike indicators to ISIs
#' @param Y A `q x n x r` or `n x r` array of 0/1 spike trains
#' where `q` is the number of cells,
#' `n` is the number of time bins, and `r` is the number of trials.
#' @param dt Scalar of length of the bins.
#' @param i Index of the neuron for which ISI will be obtained. Only need
#' to be provided if Y is 3D.
#' @return A vector of ISIs.
#' @export
bin2isi <- function(Y, dt, i=NULL){
  if (length(dim(Y))==3){ # if more than 1 trial
      if (dim(Y)[3]>1){
        isi <- unlist(apply(Y[i,,], 2,
                                   function(y) {diff(which(y==1)*dt)} ))
      }
      else{ # only one trial
        y_i <- Y[i,,]
        isi <- diff(which(y_i==1)*dt)
      }
  }
  else if (length(dim(Y))==2){
      isi <- unlist(apply(Y, 2,
                    function(y) { diff(which(y==1)*dt) }))
  }
  else{
      stop("Check the dimension of Y.")
  }
  return(isi)
}

#' Convert binary spike data to spike time data
#'
#' @param spikes A vector of 0s and 1s with 1 indicating a spike, or a matrix / data frame
#' of `r` rows of such vectors, with `r` being the number of trials.
#' @param dt A scalar of length of the time bins.
#' @return A vector or list of spike times.
#' @export

bin2num <- function(spikes, dt){
  if (is.vector(spikes)){
    output <- which(spikes==1)*dt
  }
  else if (is.matrix(spikes) | is.data.frame(spikes)){
    n_trials <- nrow(spikes)
    output <- vector(mode = "list", length = n_trials)
    for (i in 1:n_trials){
      output[[i]] <- which(spikes[i,]==1)*dt
    }
  }
  else{
    stop("Spikes must be of type vector, matrix, or data frame.")
  }
  return(output)
}


#' Convert multiple neurons' data of spike times into a binned format
#'
#' @param data A list of data frames, each being a neuron's data from multiple trials.
#' Each column of the data frame is a trial containing spike times.
#' @param dt Time bin length
#' @param trial_len Length of a trial (in sec)
#' @param binary Should the data be binned into a binary format?
#' @return A q x n x r array, where q is the number of neurons, n is the number of time bins,
#' and r is the number of trials. The elements in the array are the number of spikes for the
#' corresponding neuron in the corresponding bin and trial.
#' @export

num2bin <- function(data, dt, trial_len, binary=TRUE){
  q <- length(data)
  r <- length(data[[1]])
  n <- ceiling(trial_len/dt)
  output <- array(0, dim = c(q, n, r))
  for (i in 1:q){
    all_trials <- data[[i]] # all data of neuron i
    for (j in 1:r){
      trial_num <- all_trials[,j] # single trial data in spike times
      trial_bin <- rep(0, n) # single trial data in the binary form
      if (binary) {
        which_spikes <- ceiling(trial_num/dt)
        trial_bin[which_spikes[!is.na(which_spikes)]] <- 1
        output[i, ,j] <- trial_bin
      } else{
        output[i, ,j] <- sapply(1:n,
               function(u) {
                 sum(trial_num/dt <= u & trial_num/dt > u-1, na.rm=TRUE)
               }) # Count the number of spikes in each (coarse) bin
      }
    }
  }
  output
}

