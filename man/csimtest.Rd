\name{csimtest}
\alias{csimtest}
\title{ Simultaneous Comparisons Based on Parameter Estimates. }
\description{
Computes multiplicity adjusted p-values based on
precomputed parameter estimates.
}
\usage{
csimtest(estpar, df, covm, cmatrix = NULL, 
         ctype = "user-defined", ttype = c("free", "logical"), 
         alternative = c("two.sided", "less", "greater"), 
         asympt = FALSE, eps = 0.001, maxpts = 1e+06)
}
\arguments{
  \item{estpar}{ estimated parameter vector. }
  \item{df}{ degrees of freedom. }
  \item{covm}{ estimated covariance matrix of \code{estpar}. }
  \item{cmatrix}{ contrast matrix. }
 \item{ctype}{ a string decribing the kind of contrast matrix used. Only
used for printing in \code{\link{print.hmtestp}}. }
  \item{ttype}{Specifies whether the logical contraint method 
of Westfall (1997) will be used, or whether the uncontrained 
method will be used.}
  \item{alternative}{ the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.  You can specify just the initial letter. }
\item{asympt}{a logical indicating whether the (exact) t-distribution or
the normal approximation should be used.}
  \item{eps}{absolute error tolerance as double.}
  \item{maxpts}{maximum number of function values as integer.}
}
\details{
  This is the work-horse for \code{\link{simtest}}. If only parameter
estimates and estimates of the (co)variances are available, adjusted
p-values can be computed with this low-level function. See
\code{\link{cholesterol}} for an example.

}
\value{
an object of class \code{hmtestp}
}
\references{ 
Peter Westfall (1997), Multiple testing of general contrasts using logical
constraints and correlations, \emph{Journal of the American Statistical
Association}, \bold{92}(437), 299--306.

Frank Bretz, Alan Genz and Ludwig A. Hothorn (2001), On the numerical 
availability of multiple comparison procedures. \emph{Biometrical Journal},
\bold{43}(5), 645--656.

}
\author{ Frank Bretz <bretz@ifgb.uni-hannover.de> and  \cr 
 Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\keyword{htest}
