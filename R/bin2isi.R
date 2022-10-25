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

