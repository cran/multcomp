\name{litter}
\alias{litter}
\docType{data}
\title{ Litter Weights Data Set }
\usage{data("litter")}
\description{
  Dose response of litter weights in rats.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{dose}{dosages at four levels: \code{0}, \code{5}, \code{50}, 
                \code{500}.}
    \item{gesttime}{gestation time as covariate.}
    \item{number}{number of animals in litter as covariate.}
    \item{weight}{response variable: average post-birth weights 
                  in the entire litter.}
  }
}
\details{

  Pregnant mice were divided into four groups and the compound in four
  different doses was administered during pregnancy. Their litters
  were evaluated for birth weights.

}
\source{

  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc., page 109.

  P. H. Westfall (1997). Multiple Testing of General Contrasts Using
  Logical Constraints and Correlations. \emph{Journal of the American
  Statistical Association}, \bold{92}(437), 299--306.

}
\examples{

  ### fit ANCOVA model to data
  amod <- aov(weight ~ dose + gesttime + number, data = litter)

  ### define matrix of linear hypotheses for `dose'
  doselev <- as.integer(levels(litter$dose))
  K <- rbind(contrMat(table(litter$dose), "Tukey"),
             otrend = c(-1.5, -0.5, 0.5, 1.5),
             atrend = doselev - mean(doselev),
             ltrend = log(1:4) - mean(log(1:4)))

  ### set up multiple comparison object
  Kht <- glht(amod, linfct = mcp(dose = K), alternative = "less")

  ### cf. Westfall (1997, Table 2)
  summary(Kht, test = univariate())
  summary(Kht, test = adjusted("bonferroni"))
  summary(Kht, test = adjusted("Shaffer"))
  summary(Kht, test = adjusted("Westfall"))
  summary(Kht, test = adjusted("single-step"))

}
\keyword{datasets}
