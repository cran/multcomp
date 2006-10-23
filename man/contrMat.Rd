\name{contrMat}
\alias{contrMat}
\title{ Contrast Matrices }
\description{
  Computes contrast matrices for several multiple comparison procedures.
}
\usage{
contrMat(n, type=c("Dunnett", "Tukey", "Sequen", "AVE", 
                   "Changepoint", "Williams", "Marcus", 
                   "McDermott"), base = 1)
}
\arguments{
  \item{n}{ a (possibly named) vector of sample sizes for each group or a
            vector of factor levels. }
  \item{type}{ type of contrast. }
  \item{base}{ an integer specifying which group is considered the baseline
          group for Dunnett contrasts.}
}
\details{
  Computes the requested matrix of contrasts for comparisons of mean levels.
}
\value{
  The matrix of contrasts with appropriate row names is returned. 
}
\references{
  Frank Bretz, Alan Genz and Ludwig A. Hothorn (2001), On the numerical 
  availability of multiple comparison procedures. \emph{Biometrical Journal},
  \bold{43}(5), 645--656.
} 
\examples{
 n <- c(10,20,30,40)
 names(n) <- paste("group", 1:4, sep="")
 contrMat(n)	# Dunnett is default
 contrMat(n, base = 2)	# use second level as baseline
 contrMat(n, type = "Tukey")
 contrMat(n, type = "Sequen")
 contrMat(n, type = "AVE")
 contrMat(n, type = "Changepoint")
 contrMat(n, type = "Williams")
 contrMat(n, type = "Marcus")
 contrMat(n, type = "McDermott")
}
\keyword{misc}
