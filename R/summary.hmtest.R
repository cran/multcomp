# $Id: summary.hmtest.R,v 1.11 2003/11/19 18:24:50 hothorn Exp $

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
    if (x$asympt) 
      writeLines(strwrap(paste("Asymptotic simultaneous ", format(100 * conf.level),
                         "\% confidence intervals: ", type, sep=""),
                         prefix="\t"))
    else 
      writeLines(strwrap(paste("Simultaneous ", format(100 * conf.level),
                         "\% confidence intervals: ", type, sep=""),
                         prefix="\t"))
    cat("\n")
    if (!is.null(x$DNAME)) {
      cat("Call: \n")
      print(x$DNAME)
      cat("\n")
    }
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

            a <- switch(x$alternative, 
                        "two.sided" = c((1 - conf.level)/2,
                                         1 - (1 - conf.level)/2),
                        "less" = c(NA, conf.level),
                        "greater" = c(1 - conf.level, NA))
            a[!is.na(a)] <- paste(round(a[!is.na(a)]*100, 1), "%")
            a[is.na(a)] <- "--"

            colnames(ecout) <- c("Estimate", a, "t value",
                                 "Std.Err.", "p raw", "p Bonf", "p adj")   
            print(ecout)
        }
    }
    invisible(x)
}




