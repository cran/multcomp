# $Id: print.hmtestp.R,v 1.7 2003/05/13 10:29:13 hothorn Exp $

print.hmtestp <- function(x, digits=4, ...)
{
    digits <- min(digits, getdigits(x$eps))
    cat("\n")
    if (!is.na(x$ctype))
      type <- paste(x$ctype,"contrasts")
    else
      type <- "user-defined contrasts"
    if (x$asympt)
      writeLines(strwrap(paste("Asymptotic simultaneous tests:", type),
                         prefix="\t"))
    else
      writeLines(strwrap(paste("Simultaneous tests:", type),
                         prefix="\t"))
    cat("\n")
    if (!is.null(x$DNAME)) {
      cat("Call: \n")
      print(x$DNAME)
      cat("\n")
    }
    cat("Contrast matrix:")
    cat("\n")
    print(x$cmatrix)
    if (!is.null(x$p.value.adj)) {
        padj <- round(x$p.value.adj, digits=digits)
        cat("\nAdjusted P-Values\n")
        cat("\n")
        colnames(padj) <- "p adj"
        if (is.null(rownames(x$estimate)))
           rownames(padj) <- rownames(x$estimate)
        print(padj)
    }
}
