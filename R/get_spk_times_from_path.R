#' Get spike times from the latent path and threshold.
#' @param path Length `N` vector of the path.
#' @param k Scalar threshold value for spike generating when the path crosses
#' the threshold.
#' @return A vector of integers representing at which bins spikes occur.
#' @export

get_spk_times_from_path <- function(path, k){
  temp <- floor(path/k)
  temp[which(temp<0)] <- 0
  spk_bin_ind <- sapply(1:max(temp),
                        function(u){
                          which(temp==u)[1]
                        })
  spk_bin_ind
}