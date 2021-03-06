\name{detergent}
\alias{detergent}
\docType{data}
\title{ Detergent Durability Data Set }
\usage{data("detergent")}
\description{
  Detergent durability in an incomplete two-way design.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{detergent}{detergent, a factor at levels \code{A}, \code{B}, 
                     \code{C}, \code{D}, and \code{E}.}
    \item{block}{block, a factor at levels \code{B_1}, ..., \code{B_10}.}
    \item{plates}{response variable: number of plates washed before 
                 the foam disappears.}
  }
}
\details{

  Plates were washed with five detergent varieties, in ten blocks. A complete
  design would have 50 combinations, here only three detergent varieties 
  in each block were applied in a balanced incomplete block design. Note that
  there are six observations taken at each detergent level.

}
\source{

  H. Scheffe (1959). \emph{The Analysis of Variance}. 
  New York: John Wiley & Sons, page 189.

  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc., page 189.
}
\examples{

  ### set up two-way ANOVA without interactions
  amod <- aov(plates ~ block + detergent, data = detergent)

  ### set up all-pair comparisons
  dht <- glht(amod, linfct = mcp(detergent = "Tukey"))

  ### see Westfall et al. (1999, p. 190)
  confint(dht)

  ### see Westfall et al. (1999, p. 192)
  summary(dht, test = univariate())
  summary(dht, test = adjusted("Shaffer"))
  summary(dht, test = adjusted("Westfall"))

}
\keyword{datasets}
