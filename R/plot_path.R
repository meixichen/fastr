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
