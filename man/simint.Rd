\name{simint}
\alias{simint}
\alias{simint.default}
\alias{simint.formula}
\title{ Simultaneous Confidence Intervals }
\description{Computes simultaneous intervals for several multiple procedures.}

\usage{
\method{simint}{default}(y, x=NULL, type=c("Dunnett", "Tukey",
       "Sequen", "AVE", "Changepoint", "Williams", "Marcus",
       "McDermott","Tetrade"), cmatrix=NULL, conf.level=0.95,
       alternative=c("two.sided","less", "greater"), asympt=FALSE,
       eps=0.001, maxpts=1000000, nlevel=NULL, ...)
\method{simint}{formula}(formula, data=list(), subset, na.action, ...)
}

\arguments{
  \item{y}{a numeric vector of responses.}
  \item{x}{a numeric matrix, the design matrix.}
  \item{type}{the type of contrast to be used.}
  \item{cmatrix}{the contrast matrix itself can be specified. If
\code{cmatrix} is defined, \code{type} is ignored.}
  \item{conf.level}{confidence level.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.  You can specify just the initial letter.}
  \item{asympt}{a logical indicating whether the (exact) t-distribution or
the normal approximation should be used.}
  \item{eps}{absolute error tolerance as double.}
  \item{maxpts}{maximum number of function values as integer.}
  \item{nlevel}{a vector containing the number of levels for each factor for
\code{type == "Tetrade"}.}
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{an optional data frame containing the variables in the model.
              By default the variables are taken from
              \code{Environment(formula)}, typically the environment from which
              \code{simint} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be
          used.}
  \item{na.action}{a function which indicates what should happen when the
         data contain \code{NA}'s.  Defaults to
  \code{GetOption("na.action")}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
Computes simultaneous confidence intervals for several multiple comparisons. 
The 
implemented algorithms take the stochastical correlations between the test statistics
into account. Only single step comparisons are performed. The present function allows
for multiple comparisons of generally correlated means in general linear models under
the classical ANOVA assumptions, as well as more general approximate procedures for 
approximately normal and generally correlated parameter estimates. Either multivariate 
normal or multivariate t statistics can be used.
The interface allows the use of the multiple comparison procedures as for
example Dunnett and Tukey. 
The resulting confidence intervals are not associated with the p-values from
\code{\link{simtest}}.
}
\value{ an object of class \code{hmtest} }
\references{ 

Frank Bretz, Alan Genz and Ludwig A. Hothorn (2001), On the numerical 
availability of multiple comparison procedures. \emph{Biometrical Journal},
\bold{43}(5), 645--656.

}
\author{ Frank Bretz <bretz@ifgb.uni-hannover.de> and  \cr 
 Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{
data(recover)

# one-sided simultaneous confidence intervals for Dunnett 
# in the one-way layout
summary(simint(minutes~blanket, data=recover, type="Dunnett", conf.level=0.9, 
       alternative="less",eps=0.0001))

}
\keyword{htest}
