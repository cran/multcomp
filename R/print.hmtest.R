# $Id: print.hmtest.R,v 1.8 2003/11/19 18:24:50 hothorn Exp $

print.hmtest <- function(x, digits=4, ...)
{
    digits <- min(digits, getdigits(x$eps))
    cat("\n")
    if (!is.na(x$ctype))
      type <- paste(x$ctype,"contrasts")
    else
      type <- "user-defined contrasts"
    if (x$asympt)
      writeLines(strwrap(paste("Asymptotic simultaneous confidence intervals:", type),
                         prefix="\t"))
    else 
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
            a <- switch(x$alternative, 
                        "two.sided" = c((1 - conf.level)/2, 
                                         1 - (1 - conf.level)/2),
                        "less" = c(NA, conf.level),
                        "greater" = c(1 - conf.level, NA))
            a[!is.na(a)] <- paste(round(a[!is.na(a)]*100, 1), "%")
            a[is.na(a)] <- "--"
            colnames(ecout) <- c("Estimate", a)
            print(ecout)
        }
    }
}
