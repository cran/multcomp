\name{litter}
\alias{litter}
\non_function{}
\title{ Litter Weights Data Set }
\usage{data(litter)}
\description{
  Dose response of litter weights in rats; data set taken from Westfall et al. (1999, p. 109). Specific dose comparisons in a dose response experiment with two covariates. 
}
\format{
  This data frame contains the following variables
  \describe{
    \item{dose}{Dosages at 4 levels: \code{0}, \code{5}, \code{50}, 
                \code{500}.}
    \item{gesttime}{Gestation time as covariate.}
    \item{number}{Number of animals in litter as covariate.}
    \item{weight}{Response variable: average post-birth weights 
                  in the entire litter.}
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
data(litter)

# define the contrast matrix
C <- matrix(c(0,0,0,0,0,0,1,1,1,.75,.384,.887,-1,0,0,.25,.37,.113,
       0,-1,0,-.25,.246,-.339,0,0,-1,-.75,-1,-.661,0,0,0,
       0,0,0,0,0,0,0,0,0), ncol=7)
# numerate the contrasts
rownames(C) <- paste("C", 1:nrow(C), sep="")

# simultaneous confidence intervals for above comparisons 
# in an ANCOVA
summary(simint(weight ~ dose + gesttime + number, data=litter,
               alternative="greater", cmatrix=C))
}
\keyword{datasets}
