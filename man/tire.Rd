\name{tire}
\alias{tire}
\non_function{}
\title{ Tire Wear Data Set }
\usage{data(tire)}
\description{
  Operating costs per mile for two brands of tires; data set 
  taken from Westfall et al. (1999, p. 109). Confidence bands for 
  difference of regression functions.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{make}{Tire brand at 2 levels: \code{A} and \code{B}.}
    \item{mph}{Miles per Hour as predictor variable.}
    \item{cost}{Response variable: operating costs per mile.}
  }
}
\details{
  See Westfall et al. (1999, p. 109)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{
data(tire)
C <- c(0,1,-1,0,10,-10)
for ( x in seq(15,70,5) ) { C <- rbind( C,c(0,1,-1,0,x,-x) ) }
# numerate the contrasts
rownames(C) <- paste("C", 1:nrow(C), sep="")

# simultaneous confidence intervals for difference of regression functions
summary(simint(cost ~ make*mph, data=tire,
               cmatrix=C, eps=0.001))
}
\keyword{datasets}
