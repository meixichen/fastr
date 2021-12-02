#' Pairwise correlation based on spike counts
#' 
#' @param neuron1 A length `n` vector of neuron 1's binary spiking data in one trial.
#' @param neuron2 A length `n` vector of neuron 2's binary spiking data in one trial.
#' @param block_len Number of bins to combine into blocks in which the number of spikes is counted.
#' @return Correlation between two neurons' spiking counts
#' @export

pairwise_corr <- function(neuron1, neuron2, block_len){
  if (length(neuron1)!=length(neuron2)){
    stop("Data of neuron 1 and 2 must have the same length.")
  }
  n_bins <- length(neuron1)
  n_block <- ceiling(n_bins/block_len)
  start_vec <- seq(from=1, to=n_bins, by=block_len)
  data1 <- rep(0, n_block)
  data2 <- rep(0, n_block)
  for (i in 1:n_block){
    start <- start_vec[i]
    data1[i] <- sum(neuron1[start:min(i*block_len, n_bins)])
    data2[i] <- sum(neuron2[start:min(i*block_len, n_bins)])
  }
  cor(data1, data2)
}
