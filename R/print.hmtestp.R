# $Id: print.hmtestp.R,v 1.4 2002/07/05 16:35:57 hothorn Exp $

print.hmtestp <- function(x, digits=4, ...)
{
    digits <- min(digits, getdigits(x$eps))
    cat("\n")
    if (!is.na(x$ctype))
      type <- paste(x$ctype,"contrasts")
    else
      type <- "user-defined contrasts"
    writeLines(strwrap(paste("Simultaneous tests:", type),
                       prefix="\t"))
    cat("\n")
    cat("Call: \n")
    print(x$DNAME)
    cat("\n")
    cat("Contrast matrix:")
    cat("\n")
    print(x$cmatrix)
    if (!is.null(x$p.value.adj)) {
        padj <- round(x$p.value.adj, digits=digits)
        cat("\nAdjusted P-Values\n")
        cat("\n")
        colnames(padj) <- "p adj"
        print(padj)
    }
}
