\name{waste}
\alias{waste}
\non_function{}
\title{ Industrial Waste Data Set }
\usage{data(waste)}
\description{
  Industrial waste data set taken from Westfall et al. (1999, p. 177). 
Simultaneous confidence intervals for interaction contrasts in 
a two-way layout with interactions.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{temp}{Temperature at 3 levels: \code{low}, \code{medium},
\code{high}.}
    \item{envir}{Environment at 5 levels: \code{env1} \dots \code{env5}.}
    \item{waste}{Response variable: waste output in a manufacturing plant.}
  }
}
\details{
  See Westfall et al. (1999, p. 177)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{
data(waste)
summary(aov(waste ~ envir + temp + envir*temp, data=waste))

summary(simint(waste ~ envir:temp, data=waste,
               type="Tetrade", eps=0.1))
summary(simtest(waste ~ envir:temp, data=waste,
                type="Tetrade", eps=0.1))

}

}
\keyword{datasets}
