\name{cftest}
\alias{cftest}
\title{ Testing Estimated Coefficients }
\description{
  A convenience function for univariate testing via z- and t-tests of 
  estimated model coefficients
}
\usage{
cftest(model, ...)
}
\arguments{
  \item{model}{a fitted model.}
  \item{\dots}{additional arguments passed to \code{\link{summary.glht}}.}
}
\details{
  The usual z- or t-tests are tested without adjusting for 
  multiplicity.
}
\value{
  An object of class \code{summary.glht}.
}
\seealso{
  \code{\link[lmtest]{coeftest}}
}
\examples{

  ## The function is currently defined as
  function(model, ...)
      summary(glht(model), test = univariate(), ...)

  lmod <- lm(dist ~ speed, data = cars)
  summary(lmod)
  cftest(lmod)
  
}
\keyword{htest}
