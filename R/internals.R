# $Id: internals.R,v 1.8 2003/05/12 15:23:20 hothorn Exp $

ct <- function(x, contrasts=FALSE) {
  a <- diag(length(x))
  colnames(a) <- x
  a
}

getdigits <- function(x) {
  if (x > 0.1 || x <= 0) return(NA)
  if (1/x >= 10 && 1/x < 100) return(1)
  if (1/x >= 100 && 1/x < 1000) return(2)
  if (1/x >= 1000 && 1/x < 10000) return(3)
  if (1/x >= 10000) return(4)
}

nicepaste <- function(x, pattern) {
  if (length(x) == 1) RET <- x
  if (length(x) == 2) RET <- paste(x[1], pattern, x[2], collapse="", sep=" ")
  if (length(x) > 2)
    RET <- paste(paste(x[1:(length(x)-1)], pattern, collapse="", sep=" "),
                       x[length(x)], collapse="", sep=" ")
  RET
}
