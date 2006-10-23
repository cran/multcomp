\name{glht-methods}
\alias{summary.glht}
\alias{confint.glht}
\alias{coef.glht}
\alias{vcov.glht}
\alias{univariate}
\alias{adjusted}
\alias{Ftest}
\alias{Chisqtest}
\title{ Methods for General Linear Hypotheses }
\description{
  Simultaneous tests and confidence intervals for general linear 
  hypotheses.
}
\usage{
\method{summary}{glht}(object, test = adjusted(), ...)
\method{confint}{glht}(object, parm, level = 0.95, ...)
\method{coef}{glht}(object, rhs = FALSE, ...)
\method{vcov}{glht}(object, ...)
univariate()
adjusted(type = c("free", "Shaffer", "Westfall", p.adjust.methods), 
         ...)
Ftest()
Chisqtest()
}
\arguments{
  \item{object}{ an object of class \code{\link{glht}}.}
  \item{test}{ a function for computing p values.}
  \item{parm}{ additional parameters, currently ignored.}
  \item{level}{ the confidence level required.}
  \item{rhs}{logical, indicating whether the linear function
              \eqn{K \hat{\beta}} or the right hand side 
              \eqn{m} (\code{rhs = TRUE}) of the linear hypothesis
              should be returned.}
  \item{type}{ the multiplicity adjustment (\code{adjusted}) 
               to be applied. See below and \code{\link{p.adjust}}.}
  \item{...}{ additional arguments, such as \code{maxpts}, 
              \code{abseps} or \code{releps} to
              \code{\link[mvtnorm]{pmvnorm}} in \code{adjusted} or 
              \code{\link[mvtnorm]{qmvnorm}} and \code{adjusted = FALSE} 
              in \code{confint}.}
}
\details{

  The methods for general linear hypotheses as described by objects returned
  by \code{\link{glht}} can be used to actually test the global
  null hypothesis, each of the partial hypotheses and for
  simultaneous confidence intervals for the linear function $K \beta$.
  
  The \code{\link{coef}} and \code{\link{vcov}} methods compute the linear
  function $K \hat{\beta}$ and its covariance, respectively.

  The \code{test} argument to \code{summary} takes a function specifying
  the type of test to be applied. Classical Chisq (Wald test) or F statistics
  for testing the global hypothesis $H_0$ are implemented in functions
  \code{Chisqtest} and \code{Ftest}. Several approaches to multiplicity adjusted p 
  values for each of the linear hypotheses are implemented 
  in function \code{adjusted}. The \code{type}
  argument to \code{adjusted} specifies the method to be applied:
  \code{"free"} implements adjusted p values based on the joint
  normal or $t$ distribution of the linear function, and
  \code{"Shaffer"} and \code{"Westfall"} implement logically constraint 
  multiplicity adjustments (Shaffer, 1986; Westfall, 1997). 
  In addition, all adjustment methods
  implemented in \code{\link{p.adjust}} are available as well.

  Simultaneous confidence intervals for linear functions can be computed
  using method \code{\link{confint}}. Univariate confidence intervals
  can be computed by specifying the additional argument \code{adjusted = FALSE}
  to \code{confint}.

  All simultaneous inference procedures implemented here control
  the family-wise error rate (FWER). Multivariate
  normal and $t$ distributions, the latter one only for models of 
  class \code{\link{lm}}, are evaluated using the procedures
  implemented in package \code{mvtnorm}.

}
\value{

  \code{summary} computes (adjusted) p values for general linear hypotheses,
  \code{\link{confint}} computes (adjusted) confidence intervals. 
  \code{\link{coef}} returns estimates of the linear function $K \beta$
  and \code{\link{vcov}} its covariance. 

}
\references{

   Juliet P. Shaffer (1986), 
   Modified sequentially rejective multiple test procedures. 
   \emph{Journal of the American Statistical Association},
   \bold{81}, 826--831.

   Peter H. Westfall (1997), 
   Multiple testing of general contrasts using logical constraints
   and correlations. \emph{Journal of the American Statistical Association},
   \bold{92}, 299--306.

}
\examples{

  ### set up a two-way ANOVA with interactions
  amod <- aov(breaks ~ wool * tension, data = warpbreaks)

  ### set up all-pair comparisons for factor `tension'
  wht <- glht(amod, linfct = mcp(tension = "Tukey"))

  ### 95\% simultaneous confidence intervals
  plot(print(confint(wht)))

  ### the same (for balanced designs only)
  TukeyHSD(amod, "tension")

  ### corresponding adjusted p values
  summary(wht)


  ### confidence bands for a simple linear model, `cars' data
  plot(cars, xlab = "Speed (mph)", ylab = "Stopping distance (ft)",
       las = 1)

  ### fit linear model and add regression line to plot
  lmod <- lm(dist ~ speed, data = cars)
  abline(lmod)

  ### a grid of speeds
  speeds <- seq(from = min(cars$speed), to = max(cars$speed), 
                length = 10)

  ### linear hypotheses: 10 selected points on the regression line != 0
  K <- cbind(1, speeds)                                                        

  ### set up linear hypotheses
  cht <- glht(lmod, linfct = K)

  ### confidence intervals, i.e., confidence bands, and add them plot
  cci <- confint(cht)
  lines(speeds, cci$confint[,"lwr"], col = "blue")
  lines(speeds, cci$confint[,"upr"], col = "blue")


  ### simultaneous p values for parameters in a Cox model
  if (require("survival") && require("MASS")) {
      data("leuk", package = "MASS")
      leuk.cox <- coxph(Surv(time) ~ ag + log(wbc), data = leuk)

      ### set up linear hypotheses
      lht <- glht(leuk.cox, linfct = diag(length(coef(leuk.cox))))

      ### adjusted p values
      print(summary(lht))
  }

}
\keyword{htest}
