\name{simtest}
\alias{simtest.default}
\alias{simtest.formula}
\alias{simtest.lm}
\alias{simtest.glm}
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
        asympt=FALSE, ttype=c("free","logical"), eps=0.001,
        maxpts=1e+06, nlevel=NULL, nzerocol=c(0,0),...)
\method{simtest}{formula}(formula, data=list(), subset, na.action, whichf, ...)
\method{simtest}{lm}(y, psubset=NULL, cmatrix = NULL, ttype=c("free","logical"),
         alternative=c("two.sided","less","greater"), asympt=FALSE,
         eps=0.001, maxpts=1000000, ...)
}

\arguments{
  \item{y}{a numeric vector of responses or an object of class
           \code{\link{lm}} or \code{\link{glm}}.}
  \item{x}{a numeric matrix, the design matrix.}
  \item{type}{the type of contrast to be used. If type is not specified, cmatrix has to be specified.}
  \item{cmatrix}{the contrast matrix itself can be specified. If
\code{cmatrix} is defined, \code{type} is ignored.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or
    \code{"less"}.  You can specify just the initial letter.}
  \item{asympt}{a logical indicating whether the (exact) t-distribution or
the normal approximation should be used.}
  \item{ttype}{Specifies whether the logical contraint method of Westfall (1997) will be used, or whether the uncontrained method will be used.}
  \item{eps}{absolute error tolerance as double.}
  \item{maxpts}{maximum number of function values as integer.}
  \item{nlevel}{a vector containing the number of levels for each factor for
    \code{type == "Tetrade"}.}
  \item{nzerocol}{ a vector of two elements defining the number of
     zero-columns to add to the contrast matrix from left (the first element,
     usually 1 for the intercept) and right (usually 0 if no covariables are
     specified). \code{nzerocol} is automatically determined by
     \code{simint.formula}.}
  \item{psubset}{a vector of integers or characters indicating for which
                 subset of coefficients of a (generalized) linear model   
                 \code{y}  simultaneous p-values should be
                 computed.}
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
  \item{whichf}{if more than one factor is given in the right hand side of
    \code{formula}, \code{whichf} can be used to defined the factor to compute
    contrasts of.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
Computes multiplicity adjusted p-value for several multiple comparisons. The 
implemented algorithms take the logical relationships between the hypotheses 
and the stochastical correlations between the test statistics into account. 
Logical information is included via the methods described by Westfall (1997).
Stochastic information is included via the \code{\link[mvtnorm]{pmvt}} 
function. The p-values are 
generally the same as the come out in a closed test procedure using max-T-type 
statistics.  The procedure differs in a very subtle way from closed testing, 
but still controls FWE strongly under point null configurations; see Westfall
(1997).
The present function allows for multiple comparisons of generally correlated 
means in general linear models under the classical ANOVA assumptions, as well 
as more general approximate procedures for approximately normal and generally 
correlated parameter estimates.  Either multivariate normal or multivariate t 
statistics can be used. The interface allows the use of the multiple comparison 
procedures as for example Dunnett and Tukey. The resulting p-values are 
not associated  with the confidence intervals from \code{\link{simint}}.

The formula interfaces to \code{simtest} and \code{\link{simint}} are 
able to work with the following situations at the right hand side (the left
hand side is one continuous variable).

As long as the contrasts are specified for one single factor of interest,
any ANOVA or ANCOVA model can be used. If any of the covariables is again a
factor, specify the factor of interest with
the \code{whichf} option. The remaining (zero) columns are added automatically
to the contrast matrix (but you can also specify the number of zero 
columns by hand through
\code{nzerocol}). One exception of supplied contrasts which involve more
than one factor are the Tetrade contrasts for the analysis of two-fold
interactions (see \code{\link{waste}} for an example). 
In this case only the two-way layout model with interactions 
may be specified on the right hand side of `formula' (continuous covariables 
are possible). If a contrast matrix is specified (via \code{cmatrix}) and
\code{whichf} is missing, the complete design matrix is derived from the
right hand side of \code{formula} is used, whenever the their dimensions
match with those of \code{cmatrix}.
Some toy examples are given in the examples section.

In all other cases \code{\link{csimtest}} or \code{\link{csimint}} 
should be used which allow a greater flexibility and more potential
situations of use (e.g. multivariate data, contrasts involving more than 1
factor, non-linear models, ...), also the user has to compute the
estimates, df and covariance matrices on his own.
  
}
  
}
\value{ an object of class \code{hmtestp} }
\references{ 
  Peter Westfall (1997), Multiple testing of general contrasts using 
  logical constraints and  correlations. \emph{Journal of the American
  Statistical Association} \bold{92}(437), 299-36.

  Frank Bretz, Alan Genz and Ludwig A. Hothorn (2001), On the numerical 
  availability of multiple comparison procedures. \emph{Biometrical Journal},
  \bold{43}(5), 64--66.

}
\author{ Frank Bretz <bretz@ifgb.uni-hannover.de> and  \cr 
 Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{
data(cholesterol)

# adjusted p-values for all-pairwise comparisons in a one-way 
# layout (tests for restricted combinations)
simtest(response ~ trt, data=cholesterol, type="Tukey", ttype="logical")

# some examples for the formula interface, statistically non-sense!

# response
y <- rnorm(21)	

# three factors
f1 <- factor(c(rep(c("A", "B", "C"), 7)))
f2 <- factor(c(rep("D", 10), rep("E", 11)))

# and one continuous covariable
x <- rnorm(21)
testdata <- cbind(as.data.frame(y), f1, f2, x)

# one factor only
summary(simtest(y ~ f1))

# one factor only, the same
summary(simtest(y ~ f1, data=testdata))

# and a continuous covariable
summary(simtest(y ~ f1 + x, data=testdata))

# without intercept
summary(simtest(y ~ f1 + x - 1, data=testdata))

# with an additional factor as covariable
# use `whichf' to specify the term in the model to 
# calculate p-values or confidence intervals for
summary(simtest(y ~ f1 + f2 + x - 1, data=testdata, whichf="f1"))

# with interaction terms
summary(simtest(y ~ f1*f2 + x - 1, data=testdata, whichf="f1"))

# inference about the interactions term
# either Tetrade contrasts 
summary(simtest(y ~ f1:f2, data=testdata, type="Tetrade"))

# or a user-defined contrast matrix
# note: this is a contrast matrix for the interactions only, 
# the column for the intercept is added automatically
simtest(y ~ f1:f2, data=testdata, cmatrix=diag(6))

# works too, if the column for the intercept is included
summary(simtest(y ~ f1:f2, data=testdata, cmatrix=cbind(0, diag(6))))

# additional covariable
summary(simtest(y ~ f1:f2 + x, data=testdata, cmatrix=diag(6)))

# again with intercept and covariables in included in cmatrix
summary(simtest(y ~ f1:f2 + x, data=testdata, 
                cmatrix=cbind(0, diag(6), 0)))

}
\keyword{htest}
