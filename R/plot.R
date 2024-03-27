#' Plot the estimated loading matrix for the fitted model
#' @param fit An object of class `fastr_fit`.
#' @param varimax Apply Varimax rotation to Lambda?
#' @param zlim Range of z (matrix) values.
#' @param legend Plot legend?
#' @param ylab Y-axis label.
#' @param xlab X-axis label.
#' @param col Color palatte. Default is a color blindness friendly palatte
#' from the viridis package.
#' @param cex.axis.x Font size of the x-axis labels.
#' @param cex.axis.y Font size of the y-axis labels.
#' @param ... Arguments supplied for the `image.plot()` function.
#' @export
plot.fastr_fit <- function(fit, neuron_lab=NULL,
			   varimax=TRUE, zlim=c(-1,1), legend=TRUE,
			   ylab="Neuron index", xlab="",
			   col=viridis::viridis(12),
			   cex.axis.x=1, cex.axis.y=0.9, ...){
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
  axis(1, at=1:n_factor, tick=FALSE, cex.axis=cex.axis.x)
  axis(2, at=1:n_cell, labels=neuron_lab, las=2, tick=FALSE,
       cex.axis=cex.axis.y)
}

#' Plot the resetting latent path of a neuron with the 95% CI
#'
#' @param path A vector of the latent path of neuron. This can be either the non-reset path
#' or the reset path (specify using the `reset` argument).
#' @param se A `length(path)` vector of the point-wise SE of the latent path.
#' @param reset TRUE if `path` is the reset path, i.e. not resetting needs to be performed.
#' @param spks An optional `length(path)` vector of 0s and 1s indicating spike occurrences.
#' Must be provided if `reset=FALSE`.
#' @param k Value of the threshold of excitation. Optional. If provided, a horizontal line will
#' be drawn indicating k.
#' @param line_size Size of the path in the plot.
#' @param alpha Transparency of the shaded credible intervals.
#' @param ci_fill Fill color for the shaded credible intervals.
#' @param ci_col Line color for the credible intervals.
#' @param pt_col Color of the points indicating spike occurrences.
#' @param pt_size Size of the points indicating spike occurrences.
#' @param hline_col Color of the horizontal line indicating k.
#' @param ann_coor A length 2 vector of the coordinate of the annotation "k".
#' @param ann_size Size of the annotation of "k".
#' @return A ggplot plot.
#' @export
plot_path <- function(path, se, reset=FALSE, spks=NULL, k=NULL,
                      line_size=1.2, alpha=0.1, ci_fill="green",
                      ci_col="black", pt_col="red", pt_size=1.2,
                      hline_col="blue", ann_coor=c(1, k+0.1), ann_size=8, ...){
  if (length(path) != length(se)) stop("path and se must have the same length.")
  if (!is.null(spks) & length(path)!=length(spks)) stop("spks and path must have the same length.")
  if (!reset){
    if (is.null(spks)) stop("Must provide `spks` if reset=FALSE.")
    path <- reset_path(path, spks)
  }
  n_bin <- length(path)
  df_path <- data.frame(time=1:n_bin,
                        path=path,
                        path_se=se,
                        spks=spks)
  out <- ggplot2::ggplot(data=df_path, aes(time, path)) +
    geom_line(size = line_size) +
    geom_ribbon(aes(ymin=path-2*path_se, ymax=path+2*path_se),
                alpha=alpha, fill = ci_fill,
                color = ci_col, linetype = "dotted")
  if (!is.null(spks)){
    out <- out +
           geom_point(data=df_path[df_path$spks==1,], aes(time, spks-1),
             color=pt_col, size=pt_size)
  }
  if (!is.null(k)){
    out <- out +
      geom_hline(yintercept=k, color=hline_col, linetype="dashed")+
      annotate(geom="text", x=ann_coor[1], y=ann_coor[2], label="k",
                  color=hline_col, size=ann_size)
  }
  out
}

#' Raster plot of spike trains
#'
#' @param data A `r` by `n` matrix / data frame, or a list of `r` vectors.
#' If a matrix / df is provided, each row is a trial discretized into `n`
#' time bins in which a spike is indicated by 1 and 0 otherwise.
#' If a list is provided, each vector is a trial and contains the spike times.
#' @param pch Point shape for `points()`.
#' @param pt_cex Point size (cex) for `points()`.
#' @param ... Additional arguments to pass to `plot()`.
#' @return A raster plot
#' @export
rasterplot <- function(data, pch = "|", pt_cex = 0.9, ...){
  if (is.matrix(data) | is.data.frame(data)){
    n_trials <- nrow(data)
    n_bins <- ncol(data)
    plot(NULL, NULL, xlim=c(1,n_bins), ylim=c(1,n_trials),
         xlab="Time (bin index)", ...)
    for (i in 1:n_trials){
      spike_ind <- which(data[i,]==1)
      n_spikes <- length(spike_ind)
      points(spike_ind, rep(i, n_spikes), pch = pch, cex=pt_cex)
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
      points(spike_times, rep(i, n_spikes), pch = pch, cex=pt_cex)
    }
  }
  else{
    stop("Data must be of type matrix, data frame, or list.")
  }
}

