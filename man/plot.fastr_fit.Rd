% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.fastr_fit}
\alias{plot.fastr_fit}
\title{Plot the estimated loading matrix for the fitted model}
\usage{
\method{plot}{fastr_fit}(
  fit,
  neuron_lab = NULL,
  varimax = TRUE,
  zlim = c(-1, 1),
  legend = TRUE,
  ylab = "Neuron index",
  xlab = "",
  col = viridis::viridis(12),
  cex.axis.x = 1,
  cex.axis.y = 0.9,
  ...
)
}
\arguments{
\item{fit}{An object of class \code{fastr_fit}.}

\item{varimax}{Apply Varimax rotation to Lambda?}

\item{zlim}{Range of z (matrix) values.}

\item{legend}{Plot legend?}

\item{ylab}{Y-axis label.}

\item{xlab}{X-axis label.}

\item{col}{Color palatte. Default is a color blindness friendly palatte
from the viridis package.}

\item{cex.axis.x}{Font size of the x-axis labels.}

\item{cex.axis.y}{Font size of the y-axis labels.}

\item{...}{Arguments supplied for the \code{image.plot()} function.}
}
\description{
Plot the estimated loading matrix for the fitted model
}
