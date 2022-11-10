#' Perform exploratory factor analysis to choose the number of factors
#' @param data A `q x d x r` array, where `q` is the number of neurons, `d` is 
#' the number of time bins, `r` is the number of trials.
#' @param dt Scalar. Length of each time bin for time discretization.
#' @param d_seq A vector of candidate `n_factor` values. See details.
#' @param fa_obj A `q x dr` matrix of increments from the result of a previous  
#' run of `choose_n_factor()`. This argument is helpful when we want to test a 
#' different set of `d_seq` but don't need to fit all single-neuron models again.
#' @param plot Output a plot of cumulative var explained versus `d`?
#' @param vertical If `plot=TRUE", should the plots be arranged vertically (or
#' horizontally)?
#' @param progress Print the fitting progress?
#' @param ...  Additional arguments supplied to `fastr_fit()` for fitting single-neuron
#' models.
#' @details Given data with `q` neurons, a single-neuron model is fitted to each neuron
#' to obtain `q` estimates of neuronal latent paths x_i for i=1,..., q. Then a factor 
#' analysis model with `d` factor is fitted to the matrix whose column i is the latent
#' path increments dx_ij for neuron i. Using the estimated dx_ij, two methods will be 
#' employed to find the appropriate `d`.
#' 
#' The 'cumvar' method calculates and, optionally, plots, 
#' the cumulative proportion of variance explained (between 0 and 1) for each
#' model with the chosen number of `d`. One can look for the "elbow" in the plot
#' to determine when an increase number of `d` no longer increases the variance 
#' explained. 
#' 
#' The 'minmaxload' method does the following:
#' 1. For each `d`, obtain the estimated loading matrix Lam(d). 
#' 2. For each column (factor) of Lam(d), find the max absolute loading on it. If the 
#' max loading is a large value, this means that the corresponding factor contributes
#' to explaining the activity of one or more neurons. As a result, for each Lam(d), 
#' we can obtain `d` max loading values. 
#' 3. For each `d`, we can obtain the min of all `d` max loading values. If the min
#' is small (say < 0.1), this means that there is at least a factor among the `d` 
#' factors that does not contribute much to explaining the commonality.
#' Since we take the min of the max loadings for each `d`, this is called the 
#' 'minmaxload' method.
#' @return 
#' 1. A vector of length `length(d_seq)` containing 
#' (1) the cumulative proportion of variance explained by the corresponding 
#' number of factors (if method='cumvar'), and (2) the min of the max loadings for the 
#' corresponding number of factors (if method='minmaxload'). 
#' 2. Plot of cumvar(d) versus d, plot of minmaxload(d) versus d.
#' 3. An array or a vector of all the single-neuron model fitted paths.
#' @export

choose_n_factor <- function(data, dt, d_seq, 
			    fa_obj=NULL,
			    plot=TRUE, vertical=TRUE, progress=TRUE, ...){

  # Function that gets the cumulative prop of variance given loading mat
  get_cumvar <- function(lmat){cumsum(colSums(lmat^2)/nrow(lmat))}
  
  if (is.matrix(fa_obj)){
    allpaths <- fa_obj
    n_cell <- dim(fa_obj)[2]
  } else{
    dims <- dim(data)
    n_cell <- dims[1]
    n_bin <- dims[2]
    n_trial <- dims[3]

    allpaths <- matrix(0, nrow=(n_bin-1)*n_trial, ncol=n_cell)

    for (i in 1:n_cell){
      if (progress) cat("Fitting neuron ", i, "...\n")
      Yi <- array(data[i,,], c(n_bin, n_trial))
      fit_i <- fastr_fit(data=Yi, dt=dt, silent=TRUE, ...)
      paths_i <- array(fit_i$paths, c(n_bin, n_trial))
      dpaths_i <- as.vector(apply(paths_i, 2, diff)) # Get the increments and vectorize
      allpaths[, i] <- dpaths_i 
    }
  } 

  cumvars <- rep(NA, length(d_seq))
  minmaxload <- rep(0, length(d_seq))
  for (j in 1:length(d_seq)){
    d <- d_seq[j]
    fa.res <- factanal(allpaths, d)
    # cumvar method
    temp <- get_cumvar(fa.res$loadings)
    cumvars[j] <- temp[d]

    # minmaxload method
    lmat <- fa.res$loadings
    if (!is.matrix(lmat)){
      minmaxload[j] <- max(abs(lmat))
    } else{
      minmaxload[j] <- min(apply(lmat, 2, function(val){max(abs(val))}))
    }
  }
  out <- list(cumvars=cumvars, minmaxload=minmaxload, increments=allpaths)
  
  if (plot){
    if (vertical) mfr <- c(3,1)
    else mfr <- c(1,3)
    par(mfrow=mfr)
    plot(d_seq, cumvars, type="l", ylab="Cumvar(d)",
	 xlab="Number of factors (d)")
    plot(d_seq[1] : tail(d_seq,2)[1], diff(cumvars), type="l",
	 ylab="Cumvar(d+1)-Cumvar(d)", 
	 xlab="Number of factors (d)")
    plot(d_seq, minmaxload, type="l",
	 ylab="Min of max loadings per factor", 
	 xlab="Number of factors (d)")
    return(invisible(out))
  }
  else return(out)

}
