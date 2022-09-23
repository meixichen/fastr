#' Plot the estimated loading matrix for the fitted model
#' @param fit An object of class `fastr_fit`.
#' @param varimax Apply Varimax rotation to Lambda?
#' @param zlim Range of z (matrix) values.
#' @param ylab Y-axis label.
#' @param xlab X-axis label.
#' @param ... Arguments supplied for the `image.plot()` function.
plot.fastr_fit <- function(fit, neuron_lab=NULL,
			   varimax=TRUE, zlim=c(-1,1),
			   ylab="Neuron index", xlab="", ...){
  n_cell <- fit$env$n_cell
  n_factor <- fit$env$n_factor
  if (is.null(neuron_lab)) neuron_lab <- 1:n_cell
  
  if (varimax) {
    mat <- fit$lmat_varimax
  } else{ 
    mat <- fit$lmat_hat
  }
  fields::image.plot(x=1:n_factor, y=1:n_cell, 
  	       z=t(mat), zlim=zlim,
	       ylab=ylab, xlab=xlab, axes=FALSE,...) 
  axis(1, at=1:n_factor, tick=FALSE)
  axis(2, at=1:n_cell, labels=neuron_lab, las=2, tick=FALSE)
}
