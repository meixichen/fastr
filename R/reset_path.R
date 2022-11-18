#' Reset the latent path
#' @param path A vector of the path which has not been reset.
#' @param spks A vector of 0s and 1s for spikes, which must be the same
#' length as `path`.
#' @return A vector of reset path according to the spikes.
#' @export
reset_path <- function(path, spks){
  if (length(path) != length(spks)){ 
	  stop("path and spks must have the same length.")
  }
  n_bin <- length(path)
  out <- path
  for (t in 1:(n_bin-1)){
    if (spks[t] == 1){
      reset_val <- path[t]
      out[(t+1):n_bin] <- path[(t+1):n_bin] - reset_val
    }
  }
  out
}

