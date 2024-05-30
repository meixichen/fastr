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


#' Observed versus predicted cumulative spike counts as as graphical GOF test
#'
#' @param obs Length `n_bin` vector of observed spikes, where each element is
#' a binary indicator of whether a spike has occurred in that bin. 
#' @param pred `n_bin x n_sample` matrix of predicted spikes. Each column 
#' corresponds to a sample from the leave-neuron-out posterior predictive 
#' distribution, obtained using [fastr::leave_one_out_prediction()].
#' @param add_average Add the average of all lines?
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param title Optional. 
#' @param alpha_line Transparency for color of the lines.
#' @param ... Additional arguments to the `plot()` function.
#' @return Plot whose the X-axis corresponds to the observed 
#' cumulative spike count and Y-axis corresponds to the predicted cumulative  
#' spike counts. The counts are joined by lines, with each dot plotting the 
#' predicted cumulative spike count from the LNO-PP distribution against 
#' the observed value.
#' If `add_average=TRUE`, a red line corresponding to average of all lines is 
#' added to the plot.
#' @export
plot_gof_spk_csum <- function(obs, pred, add_average=TRUE,
                              xlab="Observed cumulative spike count",
                              ylab="Predicted cumulative spike count",
                              title="", alpha_line=0.5, ...){
  cum_spk_obs <- cumsum(obs)
  cum_spk_sample <- apply(pred, 2, cumsum)
  n_sample <- dim(pred)[2]
  plim <- range(as.vector(cum_spk_sample), cum_spk_obs)
  plot(cum_spk_obs, cum_spk_sample[,1], 
       xlim=plim, ylim=plim, type="l", 
       col=adjustcolor("gray", alpha.f=alpha_line),
       xlab=xlab, ylab=ylab,main=title, ...)
  avg_rep <- apply(cum_spk_sample, 1, mean)
  for (iter in 2:n_sample){
    lines(cum_spk_obs, cum_spk_sample[,iter],
          col=adjustcolor("gray", alpha.f=alpha_line))
  }
  abline(0, 1, col="green", lty="dashed", lwd=2)
  if (add_average) lines(cum_spk_obs, avg_rep, col="red", lwd=2)
}

#' Observed spike time vs histogram of predicted spike times as a graphical GOF test
#' 
#' @param obs Length `n_bin` vector of observed spikes, where each element is
#' a binary indicator of whether a spike has occurred in that bin. 
#' @param pred `n_bin x n_sample` matrix of predicted spikes. Each column 
#' corresponds to a sample from the leave-neuron-out posterior predictive 
#' distribution, obtained using [fastr::leave_one_out_prediction()].
#' Alternatively, a length `n_sample` list of predicted spike times can be 
#' provided. 
#' @param dt Scalar time scale of the spike trains. I.e., what is the length of 
#' each bin in secs.
#' If spike times are provided for `pred`, make sure they are the same time scale
#' as `dt`.
#' @param spk_ind A vector of integers \eqn{I_1,\ldots,I_N}. The spike times for
#' each of the \eqn{I_n}-th spike are plotted. If not supplied, 4 spikes are 
#' randomly selected.
#' @param plot_type "density" or "hist"?
#' @return A ggplot object of histogram or density plot.
#' @export
plot_gof_spk_times <- function(obs, pred, dt,
                               spk_ind=NULL, plot_type=c("density", "hist")){
  require(ggplot2)
  require(dplyr)
  plot_type <- match.arg(plot_type)
  obs_spk_times <- which(obs==1)*dt
  if (is.matrix(pred)){
    pred_spk_times <- apply(pred, 2, function(x) which(x==1)*dt)
    n_sam <- dim(pred)[2]
  } else if (is.list(pred)){
    pred_spk_times <- pred
    n_sam <- length(pred)
  } else{
    stop("`pred` argument should be either a `n_bin x n_sample`  
         or a length `n_sample` list of spike times.")
  }
  
  min_spk_count <- min(c(length(obs_spk_times), sapply(pred_spk_times, length)))
  if (is.null(spk_ind)){
    n_ind <- min(min_spk_count, 4)
    spk_ind <- sort(sample(1:min_spk_count, n_ind))
  } else{
    if (all(sapply(spk_ind, is.integer)) && (max(spk_ind) < min_spk_count)){
      n_ind <- length(spk_ind)
    } else{
      stop("`spk_ind` must be all integers and its largest value should
           not exceed the minimum spike count in `pred`.")
    }
  }
  
  spk_df <- data.frame(spike_time=obs_spk_times[spk_ind], 
                       spike_index=as.factor(spk_ind),
                       type=rep("Observed", n_ind))
  pred_spk_times_mat <- sapply(spk_ind, 
                              function(spk_ind){
                                sapply(pred_spk_times, function(rep) rep[spk_ind])
                              })
  spk_df <- rbind(spk_df,
                  data.frame(spike_time=as.vector(pred_spk_times_mat),
                             spike_index=as.factor(rep(spk_ind, each=n_sam)),
                             type=rep("Predicted", n_ind*n_sam)))

  spk_df_pred <- spk_df[spk_df$type=="Predicted",]
  spk_df_obs <- spk_df[spk_df$type=="Observed",]
  if (plot_type=="density"){
    ggobj <-
      ggplot(data=spk_df_pred, aes(x=spike_time, color=spike_index)) + 
      geom_density() +
      guides(color=guide_legend(title="Spike index"))
  } else{
    ggobj <- ggplot(data=spk_df_pred, 
                    aes(x=spike_time, color=spike_index, fill=spike_index)) + 
      geom_histogram() +
      guides(color="none", fill=guide_legend(title="Spike index"))
  }
    ggobj +
    geom_vline(data=spk_df_obs, 
               aes(xintercept=spike_time, color=spike_index), linetype="dotdash")+
    facet_grid(spike_index~.) +
    xlab("Spike time (in sec)")
}

