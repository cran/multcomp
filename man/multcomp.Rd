\name{multcomp}
\alias{multcomp}
\title{General Information on the multcomp Package}
\non_function{}
\description{
The \code{multcomp} package allows for multiple comparisons of
\eqn{k} groups in general linear models. We use the unifying
representations of multiple contrast tests, which include all
common multiple comparison procedures, such as the many-to-one
comparisons of Dunnett, the all-pairwise comparisons of Tukey and
many other procedures. We provide a list of 9
standard procedures (Dunnett, Tukey, Sequen, AVE, Changepoint,
Williams, Marcus, McDermott, Tetrade), where the user selects the 
comparisons of his interest. In addition, a free input interface 
for the contrast matrix allows for more special comparisons.

The comparisons itself are not restricted to balanced or simple
designs. Instead, the programs are designed to suit multiple
comparisons within the general linear model, thus allowing for
covariates, nested effects, correlated means and missing values.
The program is designed for the normal set-up with a common (possibly unknown)
variance and a known covariance matrix. But instead of using the
usual Bonferroni and Holm procedures, we take the exact
correlations between the test statistics into account by use of
the multivariate t-distribution. The resulting procedures are
therefore more powerful (the Bonferroni and Holm adjusted 
p-values are reported for reference). We also allow the user to perform an
asymptotic analysis based on the multivariate normal
distribution (as required e.g. in multiple comparisons based on
asymptotic rank transformations; assumed asymptotic normality
when comparing binomial parameters; etc.). Two functions will be
provided. The first one computes confidence intervals for the
common single step procedures (\code{\link{simint}}). 
This approach can be uniformly
improved by applying the closed testing principle, what is
implemented in the second function (\code{\link{simtest}}; 
but no confidence intervals
are available for the latter procedure). Use either \code{\link{csimint}} or
\code{\link{csimtest}} if you want to pass the estimates by hand. 

For testing and validation purposes we included some examples from Westfall
et al. (1999).
}
    
\references{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.

  Peter Westfall (1997), Multiple testing of general contrasts using logical
  constraints and correlations, \emph{Journal of the American Statistical
  Association}, \bold{92}(437), 299--306.
 
  Frank Bretz, Alan Genz and Ludwig A. Hothorn (2001), On the numerical
  availability of multiple comparison procedures. \emph{Biometrical Journal},
  \bold{43}(5), 645--656.

}
\author{
 Frank Bretz <bretz@ifgb.uni-hannover.de>,
 Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> and
 Peter Westfall <WESTFALL@ba.ttu.edu>
}
\keyword{misc}
