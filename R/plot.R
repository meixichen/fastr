#' Plot the estimated loading matrix for the fitted model
#' @param fit An object of class `fastr_fit`.
#' @param varimax Apply Varimax rotation to Lambda?
#' @param zlim Range of z (matrix) values.
#' @param legend Plot legend?
#' @param ylab Y-axis label.
#' @param xlab X-axis label.
#' @param col Color palatte. Default is a color blindness friendly palatte 
#' from the viridis package.
#' @param ... Arguments supplied for the `image.plot()` function.
#' @export
plot.fastr_fit <- function(fit, neuron_lab=NULL,
			   varimax=TRUE, zlim=c(-1,1), legend=TRUE,
			   ylab="Neuron index", xlab="", 
			   col=viridis(12), ...){
  n_cell <- fit$env$n_cell
  n_factor <- fit$env$n_factor
  if (is.null(neuron_lab)) neuron_lab <- 1:n_cell
  
  if (varimax) {
    mat <- fit$lmat_varimax
  } else{ 
    mat <- fit$lmat_hat
  }

  if (legend){
    fields::image.plot(x=1:n_factor, y=1:n_cell, 
		       z=t(mat), zlim=zlim, col=col,
		       ylab=ylab, xlab=xlab, axes=FALSE,...) 
  } else{
    image(x=1:n_factor, y=1:n_cell, 
          z=t(mat), zlim=zlim, col=col,
          ylab=ylab, xlab=xlab, axes=FALSE,...)
  }
  axis(1, at=1:n_factor, tick=FALSE)
  axis(2, at=1:n_cell, labels=neuron_lab, las=2, tick=FALSE)
}
