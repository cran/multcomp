% file VR/MASS/ginv.Rd
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
% Copied from VR_6.3-2
\name{ginv}
\alias{ginv}
\title{
Generalized Inverse of a Matrix
}
\description{
Calculates the Moore-Penrose generalized inverse of a matrix
\code{X}.
}
\usage{
ginv(X, tol=sqrt(.Machine$double.eps))
}
\arguments{
\item{X}{
Matrix for which the Moore-Penrose inverse is required.
}
\item{tol}{
A relative tolerance to detect zero singular values.
}}
\value{
A MP generalized inverse matrix for \code{X}.
}
\references{
Venables & Ripley, Chapter 2. Copied from package \code{VR}.
}
\seealso{
\code{\link{solve}}, \code{\link{svd}}, \code{\link{eigen}}
}
\examples{
\dontrun{
# The function is currently defined as
function(X, tol = sqrt(.Machine$double.eps))
{
## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] \%*\% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}
}}
\keyword{algebra}
