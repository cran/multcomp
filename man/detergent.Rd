\name{detergent}
\alias{detergent}
\non_function{}
\title{ Detergent Durability Data Set }
\usage{data(detergent)}
\description{
 Detergent durability in a fractional design experiment; data set 
 taken from Westfall et al. (1999, p. 189). Balanced incomplete block design.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{detergent}{detergent at 5 levels: 1 ... 5}
    \item{block}{block at 10 levels: 1 ... 10}
    \item{plates}{Response variable: number of plates washed before 
     the foam disappears.}
  }
}
\details{
  See Westfall et al. (1999, p. 189)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{

data(detergent)

N <- rep(2, 5)

# BIBD: prepare the contrast matrix = all-pair comparisons for 
# the 5 levels of detergent
C <- contrMat(N, type="Tukey")
# the additional 10 columns of are for the 10 blocks
C <- cbind( matrix(0, ncol=10, nrow=10), C )
# numerate the contrasts
colnames(C) <- NULL
rownames(C) <- paste("C", 1:nrow(C), sep="")

# adjusted p-values
summary(simtest(plates ~ block+detergent, data=detergent,
                ttype="logical", cmatrix=C))

}
\keyword{datasets}
