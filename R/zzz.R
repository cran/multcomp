# $Id: zzz.R,v 1.4 2003/06/23 13:35:22 hothorn Exp $

.onLoad <- function(lib, pkg) {
    if(!require(mvtnorm))
        warning("Could not load package mvtnorm")
}
