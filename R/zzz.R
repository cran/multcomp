# $Id: zzz.R,v 1.3 2002/04/08 14:18:52 hothorn Exp $

.First.lib <- function(lib, pkg) {
    if(!require(mvtnorm))
        warning("Could not load package mvtnorm")
}
