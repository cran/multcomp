# $Id: summary.hmtest.R,v 1.6 2002/03/14 09:26:44 hothorn Exp $

summary.hmtest <- function(object, ...)
{
     class(object) <- "summary.hmtest"
     object
}

print.summary.hmtest <- function(x, digits = max(3, getOption("digits")-3),
                                 ...)
{
    digits <- min(digits, getdigits(x$eps))
    cat("\n")
    if (!is.na(x$ctype))
      type <- paste(x$ctype,"contrasts")
    else
      type <- "user-defined contrasts"
    cint <- round(x$conf.int, digits=digits)
    conf.level <- attr(cint, "conf.level")
    attr(cint, "conf.level") <- NULL
    writeLines(strwrap(paste("Simultaneous ", format(100 * conf.level),
                       "\% confidence intervals: ", type, sep=""),
                       prefix="\t"))
    cat("\n")
    cat("Data: ", x$DNAME, "\n")
    cat("\n")
    cat("Contrast matrix:")
    cat("\n")
    print(x$cmatrix)
    cat("\nAbsolute Error Tolerance: ", x$eps, "\n")
    cat(paste("\n", format(100 * conf.level), "\% quantile: ",
              round(x$calpha,digits=digits), "\n"))
    if (!is.null(x$estimate) && !is.null(x$conf.int)) {
        est <- round(x$estimate, digits=digits)
        if (length(est) == nrow(cint)) {
            cat("\nCoefficients:\n")
            stat <- round(x$statistics, digits=digits)
            sd <- round(x$sd, digits=digits)
            praw <- round(x$p.value.raw, digits=digits)
            pbon <- round(x$p.value.bon, digits=digits)
            padj <- round(x$p.value.adj, digits=digits)
            ecout <- cbind(est, cint, stat, sd, praw, pbon, padj)
            colnames(ecout) <- c("Estimate", "low CI,", "upp CI", "t value",
                                 "Std.Err.", "p raw", "p Bonf", "p adj")   
            print(ecout)
        }
    }
    invisible(x)
}




