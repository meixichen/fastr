#' Raster plot of spike trains
#' 
#' @param data A `r` by `n` matrix / data frame, or a list of `r` vectors.
#' If a matrix / df is provided, each row is a trial discretized into `n`
#' time bins in which a spike is indicated by 1 and 0 otherwise. 
#' If a list is provided, each vector is a trial and contains the spike times.
#' @param ... Additional arguments to pass to `plot()`.
#' @return A raster plot
#' @export

rasterplot <- function(data, ...){
  if (is.matrix(data) | is.data.frame(data)){
    n_trials <- nrow(data)
    n_bins <- ncol(data)
    plot(NULL, NULL, xlim=c(1,n_bins), ylim=c(1,n_trials), 
         xlab="Time (bin index)", ...)
    for (i in 1:n_trials){
      spike_ind <- which(data[i,]==1)
      n_spikes <- length(spike_ind)
      points(spike_ind, rep(i, n_spikes), pch = "|")
    }
  }
  else if (is.list(data)){
    n_trials <- length(data)
    x_max <- max(sapply(data, max)) # The max value of x-axis
    plot(NULL, NULL, xlim=c(0,ceiling(x_max)), ylim=c(1,n_trials), 
         xlab="Time", ...)
    for (i in 1:n_trials){
      spike_times <- data[[i]]
      n_spikes <- length(spike_times)
      points(spike_times, rep(i, n_spikes), pch = "|")
    }
  }
  else{
    stop("Data must be of type matrix, data frame, or list.")
  }
}