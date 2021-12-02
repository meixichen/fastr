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