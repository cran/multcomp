\name{csimint}
\alias{csimint}
\title{ Simultaneous Confidence Intervals Based on Parameter Estimates }
\description{
Computes simultaneous intervals for several multiple procedures based on
precomputed parameter estimates. 
}
\usage{
csimint(estpar, df, covm, cmatrix = NULL, 
        ctype = "user-defined", conf.level = 0.95, 
        alternative = c("two.sided", "less", "greater"),
        asympt = FALSE, eps = 0.001, maxpts = 1e+06)
}
\arguments{
  \item{estpar}{ estimated parameter vector. }
  \item{df}{ degrees of freedom. }
  \item{covm}{ estimated covariance matrix of \code{estpar}. }
  \item{cmatrix}{ contrast matrix. }
  \item{ctype}{ a string decribing the kind of contrast matrix used. Only
used for printing in \code{\link{print.hmtest}}. }
  \item{conf.level}{confidence level. }
  \item{alternative}{ the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.  You can specify just the initial letter. }
\item{asympt}{a logical indicating whether the (exact) t-distribution or
the normal approximation should be used.}
  \item{eps}{absolute error tolerance as double.}
  \item{maxpts}{maximum number of function values as integer.}
}
\details{
  This is the work-horse for \code{\link{simint}}. If only parameter
estimates and estimates of the (co)variances are available, simultaneous
confidence intervals can be computed with this low-level function. See
\code{\link{cholesterol}} for an example.

}
\value{
an object of class \code{hmtest}
}
\references{

  Frank Bretz, Alan Genz and Ludwig A. Hothorn (2001), On the numerical 
  availability of multiple comparison procedures. \emph{Biometrical Journal},
  \bold{43}(5), 645--656.

}
\author{ Frank Bretz <bretz@ifgb.uni-hannover.de> and  \cr 
 Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\keyword{htest}
