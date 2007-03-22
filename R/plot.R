
### uhhh -- mainly copy and paste from plot.TukeyHSD
### with modifications by Richard M. Heiberger <rmh@temple.edu>
plot.confint.glht <- function(x, xlim, xlab, ...) {

    xi <- x$confint
    ### make sure one-sided intervals are drawn correctly
    xrange <- c(min(xi[,"lwr"]), max(xi[, "upr"]))
    if (!is.finite(xrange[1])) xrange[1] <- min(xi[,"Estimate"])
    if (!is.finite(xrange[2])) xrange[2] <- max(xi[,"Estimate"])
    yvals <- nrow(xi):1
    if (missing(xlim))
        xlim <- xrange
    plot(c(xi[, "lwr"], xi[, "upr"]), rep.int(yvals, 2), 
         type = "n", axes = FALSE, xlab = "", ylab = "", xlim = xlim, ...)
    axis(1, ...)
    axis(2, at = nrow(xi):1, labels = dimnames(xi)[[1]], 
         las = 1, ...)
    abline(h = yvals, lty = 1, lwd = 0, col = "lightgray")
    abline(v = 0, lty = 2, lwd = 0, ...)
    left <- xi[, "lwr"]
    left[!is.finite(left)] <- min(c(0, xrange[1] * 2))
    right <- xi[, "upr"]
    right[!is.finite(right)] <- max(c(0, xrange[2] * 2))
    segments(left, yvals, right, yvals, ...)
    points(xi[, "lwr"], yvals, pch = "(", ...)
    points(xi[, "upr"], yvals, pch = ")", ...)
    points(xi[, "Estimate"], yvals, pch = 20, ...)
    main <- list(...)$main
    if (is.null(main)) {
        if (attr(x, "type") == "adjusted") {
            main <- paste(format(100 * attr(x$confint, "conf.level"), 2), 
                          "% family-wise confidence level\n", sep = "")
        } else {
            main <- paste(format(100 * attr(x$confint, "conf.level"), 2),
                          "% confidence level\n", sep = "")
        }
    }
    if (missing(xlab))
          xlab <- "Linear Function"
    title(main = main, xlab = xlab)
    box()
}

plot.glht <- function(x, ...) plot(confint(x), ...)
