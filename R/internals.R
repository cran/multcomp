# $Id: internals.R,v 1.3 2002/04/08 14:18:07 hothorn Exp $

ct <- function(x, contrasts=FALSE) diag(length(x))

getdigits <- function(x) {
  if (x > 0.1 || x <= 0) return(NA)
  if (1/x >= 10 && 1/x < 100) return(1)
  if (1/x >= 100 && 1/x < 1000) return(2)
  if (1/x >= 1000 && 1/x < 10000) return(3)
  if (1/x >= 10000) return(4)
}
