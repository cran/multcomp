\name{cholesterol}
\alias{cholesterol}
\non_function{}
\title{ Cholesterol Reduction Data Set }
\usage{data(cholesterol)}
\description{
 Cholesterol reduction for five treatments; data set taken from 
 Westfall et al. (1999, p. 153). All pairwise comparisons according 
 to Tukey in a balanced one-way layout.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{trt}{Treatment at 5 levels: \code{A}, \code{B}, \code{C},
               \code{D}, \code{E}}
    \item{response}{Response variable.}
  }
}
\details{
  See Westfall et al. (1999, p. 153)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{
data(cholesterol)

# adjusted p-values for all-pairwise comparisons in a one-way layout 
# tests for restricted combinations
simtest(response ~ trt, data=cholesterol, type="Tukey",
        ttype="logical")

# adjusted p-values all-pairwise comparisons in a one-way layout 
# (tests for free combinations -> p-values will be larger)
simtest(response ~ trt, data=cholesterol, type="Tukey",
        ttype="free")

# enter now the estimates as parameters
# begin with degrees of freedom
nu <- as.integer(45)
# estimates
parm <- c(10.6151, -4.8331, -1.3901, 1.7597, 4.7461, 10.3325) 
# build the covariance matrix
N <- rep(2, 5)
contrast <- contrMat(N, type="Tukey")
covm <- rep(-0.20254649, 36)
covm <- matrix(covm, ncol=6)
covm[1,2:6] <- rep(0.02893521, 5)
covm[2:6,1] <- rep(0.02893521, 5)
covm[1,1] <- 0.14467606
for (i in  2:6) { covm[i,i] <- 0.83912115 } 

# use the work-horse directly (and add zero column for the intercept)

csimint(estpar=parm, df=nu, covm=covm, cmatrix=cbind(0, contrast))    
csimtest(estpar=parm, df=nu, covm=covm, cmatrix=cbind(0, contrast),
         ttype="logical")      

}
\keyword{datasets}
