\name{simtest}
\alias{simtest.default}
\alias{simtest.formula}
\alias{simtest}
\title{ Simultaneous Comparisons }
\description{
Computes multiplicity adjusted p-value for several multiple comparisons.
}

\usage{
\method{simtest}{default}(y, x=NULL, type=c("Dunnett", "Tukey",
        "Sequen", "AVE", "Changepoint", "Williams", "Marcus",
        "McDermott","Tetrade"), cmatrix=NULL,
        alternative=c("two.sided","less", "greater"),
        asympt=FALSE, ttype=c("free","logical"), eps=.,
        maxpts=, nlevel=NULL, ...)
\method{simtest}{formula}(formula, data=list(), subset, na.action, ...)
}

\arguments{
  \item{y}{a numeric vector of responses.}
  \item{x}{a numeric matrix, the design matrix.}
  \item{type}{the type of contrast to be used. If type is not specified, cmatrix has to be specified.}
  \item{cmatrix}{the contrast matrix itself can be specified. If
\code{cmatrix} is defined, \code{type} is ignored.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.  You can specify just the initial letter.}
  \item{asympt}{a logical indicating whether the (exact) t-distribution or
the normal approximation should be used.}
  \item{ttype}{Specifies whether the logical contraint method of Westfall (997) will be used, or whether the uncontrained method will be used.}
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
Computes multiplicity adjusted p-value for several multiple comparisons. The 
implemented algorithms take the logical relationships between the hypotheses 
and the stochastical correlations between the test statistics into account. 
Logical information is included via the methods described by Westfall (997).
Stochastic information is included via the \code{\link[mvtnorm]{pmvt}} 
function. The p-values are 
generally the same as the come out in a closed test procedure using max-T-type 
statistics.  The procedure differs in a very subtle way from closed testing, 
but still controls FWE strongly under point null configurations; see Westfall (997).
The present function allows for multiple comparisons of generally correlated 
means in general linear models under the classical ANOVA assumptions, as well 
as more general approximate procedures for approximately normal and generally 
correlated parameter estimates.  Either multivariate normal or multivariate t 
statistics can be used. The interface allows the use of the multiple comparison 
procedures as for example Dunnett and Tukey. The resulting p-values are 
not associated  with the confidence intervals from \code{\link{simint}}.
}
\value{ an object of class \code{hmtestp} }
\references{ 
  Peter Westfall (997), Multiple testing of general contrasts using 
  logical constraints and  correlations. \emph{Journal of the American
  Statistical Association} \bold{92}(437), 299-36

  Frank Bretz, Alan Genz and Ludwig A. Hothorn (2), On the numerical 
  availability of multiple comparison procedures. \emph{Biometrical Journal},
  \bold{43}(), 64--66.

}
\author{ Frank Bretz <bretz@ifgb.uni-hannover.de> and  \cr 
 Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{
data(cholesterol)

# adjusted p-values for all-pairwise comparisons in a onw-way 
# layout (tests for restricted combinations)
simtest(response ~ trt, data=cholesterol, type="Tukey", ttype="logical")

}
\keyword{htest}
