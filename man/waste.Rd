\name{waste}
\alias{waste}
\docType{data}
\title{ Industrial Waste Data Set }
\usage{data("waste")}
\description{
  Industrial waste output in a manufactoring plant.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{temp}{temperature, a factor at three levels: \code{low}, \code{medium}, \code{high}.}
    \item{envir}{environment, a factor at five levels: \code{env1} \dots \code{env5}.}
    \item{waste}{response variable: waste output in a manufacturing plant.}
  }
}
\details{

  The data are from an experiment designed to study the effect of temperature
  (\code{temp}) and environment (\code{envir}) on waste output in a manufactoring plant.
  Two replicate measurements were taken at each temperature / environment combination.

}
\source{

  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc., page 177.

}
\examples{

  ### set up two-way ANOVA with interactions
  amod <- aov(waste ~ envir + temp + envir*temp, data=waste)

  ### set up linear hypotheses for all-pairs of both factors
  wht <- glht(amod, linfct = mcp(temp = "Tukey", envir = "Tukey"))

  ### cf. Westfall et al. (1999, page 181)
  summary(wht, test = adjusted("Shaffer"))

}
\keyword{datasets}
