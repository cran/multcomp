# $Id: print.hmtest.R,v 1.6 2002/09/18 09:30:12 hothorn Exp $

print.hmtest <- function(x, digits=4, ...)
{
    digits <- min(digits, getdigits(x$eps))
    cat("\n")
    if (!is.na(x$ctype))
      type <- paste(x$ctype,"contrasts")
    else
      type <- "user-defined contrasts"
    writeLines(strwrap(paste("Simultaneous confidence intervals:", type),
                       prefix="\t"))
    cat("\n")
    if (!is.null(x$DNAME)) {
      cat("Call: \n")
      print(x$DNAME)
      cat("\n")
    }
    if (!is.null(x$estimate) && !is.null(x$conf.int)) {
        cint <- round(x$conf.int, digits=digits)
        conf.level <- attr(cint, "conf.level")
        attr(cint, "conf.level") <- NULL
        est <- round(x$estimate, digits=digits)
        if (length(est) == nrow(cint)) {
            writeLines(strwrap(paste(format(100 * conf.level),
                               "\% confidence intervals"), prefix="\t"))
            cat("\n")
            ecout <- cbind(est, cint)
            colnames(ecout) <- c("Estimate", "lower CI", "upper CI")
            print(ecout)
        }
    }
}
