# $Id: zzz.R,v 1.5 2004/08/04 08:38:22 hothorn Exp $

.onLoad <- function(lib, pkg) {
    if(!require(mvtnorm, quietly = TRUE))
        stop("Could not load package mvtnorm")
}
