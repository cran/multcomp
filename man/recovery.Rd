\name{recovery}
\alias{recovery}
\non_function{}
\title{ Recovery Time Data Set }
\usage{data(recovery)}
\description{
 Recovery time data set taken from Westfall 
 et al. (1999, p. 80). Many-to-one comparisons according to 
 Dunnett in an unbalanced one-way layout. 
}
\format{
  This data frame contains the following variables
  \describe{
    \item{blanket}{Blanket type at 4 levels: \code{b0}, \code{b1},
                   \code{b3}, \code{b4}.}
    \item{minutes}{Response variable: recovery time after a surgical 
                   procedure.}
  }
}
\details{
  See Westfall et al. (1999, p. 80)
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{
data(recovery)

# one-sided simultaneous confidence intervals for Dunnett 
# in the one-way layout
simint(minutes ~ blanket, data=recovery, conf.level=0.9, 
       alternative="less",eps=0.0001)

# same results, but specifying the contrast matrix by hand
C <- c(0, 0, 0, -1, -1, -1, 1, 0, 0, 0, 1, 0, 0, 0, 1)
C <- matrix(C, ncol=5)
# numerate the contrasts
rownames(C) <- paste("C", 1:nrow(C), sep="")
test <- simint(minutes~blanket, data=recovery, conf.level=0.9, 
               alternative="less",eps=0.0001, cmatrix=C)
print(test)

# same results, but more detailed information using the summary method
summary(test)
}
\keyword{datasets}
