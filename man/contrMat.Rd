\name{contrMat}
\alias{contrMat}
\title{ Contrast Matrices }
\description{
  Computes contrast matrices for several multiple comparision procedures.
}
\usage{
contrMat(n, type=c("Dunnett", "Tukey", "Sequen", "AVE", 
                   "Changepoint", "Williams", "Marcus", 
                   "McDermott", "Tetrade"), nlevel=NULL)
}
\arguments{
  \item{n}{ a (possibly named) vector of sample sizes for each group. }
  \item{type}{ type of contrast. }
  \item{nlevel}{ number of levels of two factors (for \code{type="Tetrade"}
                  only).}
}
\details{
  Computes the requested matrix of contrasts.
}
\value{
  The matrix of contrasts with appropriate rownames is returned. 
}
\references{
  Frank Bretz, Alan Genz and Ludwig A. Hothorn (2001), On the numerical 
  availability of multiple comparison procedures. \emph{Biometrical Journal},
  \bold{43}(5), 645--656.
} 
\author{ Frank Bretz <bretz@ifgb.uni-hannover.de> and  \cr
 Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{
 n <- c(10,20,30,40)
 a <- 2
 b <- 3
 names(n) <- paste("group", 1:4, sep="")
 contrMat(n)	# Dunnett is default
 contrMat(n, type="Tukey")
 contrMat(n, type="Sequen")
 contrMat(n, type="AVE")
 contrMat(n, type="Changepoint")
 contrMat(n, type="Williams")
 contrMat(n, type="Marcus")
 contrMat(n, type="McDermott")
 contrMat(n, type="Tetrade", nlevel=c(a,b))
}
\keyword{misc}
