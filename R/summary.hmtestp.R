# $Id: summary.hmtestp.R,v 1.6 2002/07/11 12:56:53 hothorn Exp $

summary.hmtestp <- function(object, ...)
{
     class(object) <- "summary.hmtestp"
     object
}

print.summary.hmtestp <- function(x, digits = max(3, getOption("digits")-3),
                                 ...)
{
    digits <- min(digits, getdigits(x$eps))
    cat("\n")
    if (!is.na(x$ctype))
      type <- paste(x$ctype,"contrasts")
    else
      type <- "user-defined contrasts"
    cat("\t Simultaneous tests:", type, "\n")
    cat("\n")
    cat("Call: \n")
    print(x$DNAME)
    cat("\n")
    cat("\t", type, "for factor",  x$FNAME$mainF)
    if (length(x$FNAME$coVar) > 0) {
      if (length(grep("[+]", x$FNAME$coVar)) > 0)
        cat(", covariables: ", x$FNAME$coVar, "\n")
      else
        cat(", covariable: ", x$FNAME$coVar, "\n")
    } else {
      cat("\n")
    }
    cat("\n")
    cat("Contrast matrix:")
    cat("\n")
    print(x$cmatrix)
    cat("\n")
    cat("\nAbsolute Error Tolerance: ", x$eps, "\n")
    if (!is.null(x$estimate)) {
        est <- round(x$estimate, digits=digits)
        cat("\nCoefficients:\n")
        stat <- round(x$statistics, digits=digits)
        sd <- round(x$sd, digits=digits)
        praw <- round(x$p.value.raw, digits=digits)
        pbon <- round(x$p.value.bon, digits=digits)
        padj <- round(x$p.value.adj, digits=digits)
        ecout <- cbind(est, stat, sd, praw, pbon, padj)
        colnames(ecout) <- c("Estimate", "t value",
                             "Std.Err.", "p raw", "p Bonf", "p adj")   
        print(ecout)
    }
    invisible(x)
}




