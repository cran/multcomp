\name{plot.hmtest}
\alias{plot.hmtest}
\title{ Plot Simultaneous Confidence Intervals }
\description{
  This method plots objects of class \code{hmtest}.
}
\usage{
\method{plot}{hmtest}(x, ltycint = 2, ltyzero = 3, main = NULL, 
                      xlab = NULL,...)
}
\arguments{
  \item{x}{an object of class \code{hmtest}}
  \item{ltycint}{integer specifying the line type of the confidence
intervals.}
  \item{ltyzero}{integer specifying the line type of the vertical line at zero.}
  \item{main}{an overall title of the plot.}
  \item{xlab}{a title for the \code{x} axis.}
  \item{\dots}{other parameters to be passed through to plotting functions.}
}
\details{
  The confidence intervals are plotted as horizontal lines.
}
\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\keyword{hplot}
