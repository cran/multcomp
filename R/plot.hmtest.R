# $Id: plot.hmtest.R,v 1.6 2004/03/08 07:35:03 hothorn Exp $

plot.hmtest <- function(x, ltycint=2, ltyzero=3, 
                        main = NULL, xlab = NULL, ...) {
  est <- x$estimate
  cint <- x$conf.int
  conf.level <- attr(cint, "conf.level")
  attr(cint, "conf.level") <- NULL
  if (!is.na(x$ctype))
    type <- paste(x$ctype,"contrasts")
  else
    type <- "user-defined contrasts"
  n <- length(est)
  crange <- range(cint)
  ONE <- FALSE
  if (is.infinite(crange[1])) {
    crange <- c(min(est), max(cint[,2]))
    cint[,1] <- 2*crange[1]
    ONE <- TRUE
  }
  if (is.infinite(crange[2])) {
    crange <- range(min(cint[,1]), max(est))
    cint[,2] <- crange[2]*2
    ONE <- TRUE
  }
  crange[1] <- ifelse(crange[1] > 0, crange[1]*0.9,crange[1]*1.1)
  crange[2] <- ifelse(crange[2] > 0, crange[2]*1.1, crange[2]*0.9)
  if (is.null(xlab)) {
    if (ONE)
      xlab <- paste(format(100 * conf.level), "\%", 
                 "one-sided", "confidence intervals")
    else
      xlab <- paste(format(100 * conf.level), "\%",
                 "two-sided", "confidence intervals")     
  }

  # <FIXME>
  # strwidth needs an open graphic device.
  # if (dev.cur() == 1) 
  plot.new()
  #  plot(1:n, type="n", ...)
  # </FIXME>

  # we need to determine cex.axis:
  args <- list(...)
  cex.axis <- args$cex.axis
  if (!is.null(cex.axis)) 
    par(cex.axis=cex.axis)
  # we need to determine the left margin depending on the size of the 
  # factor levels
  oldmai <- mymai <- par("mai")
  ywidth <- max(strwidth(rownames(est), units="inches", 
                         cex=par("cex.axis")))*1.2
  if (mymai[2] < ywidth)
   mymai[2] <- ywidth
  par(mai=mymai, new=TRUE)
  if (!is.null(main)) type = main

  pr = rbind(c(crange[1], 1), c(crange[2], n))
  pargs = c(list(x = pr[,1]), list(y=pr[,2]), 
            type="n", axes=FALSE, xlab=xlab, ylab = "", main=type, args)
  do.call("plot", pargs)
  axis(1, ...)
  axis(2, 1:n, rownames(est)[n:1], las=1, ...)
  box(...)
  for (i in 1:n)  {
    segments(cint[n-i+1,1], i, cint[n-i+1,2], i, lty=ltycint, ...)
    points(cint[n-i+1,1], i, pch="(", ...)
    points(cint[n-i+1,2], i, pch=")", ...)
    points(est[n-i+1], i, pch=19, ...)
  }
  abline(v = 0, lty = ltyzero, ...)
  par(mai=oldmai)
}
