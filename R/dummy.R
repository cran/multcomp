# $Id: contrMat.R,v 1.18 2003/05/13 11:21:20 hothorn Exp $

contrMat <- function(n, type=c("Dunnett", "Tukey", "Sequen", "AVE",
                               "Changepoint", "Williams", "Marcus",
                               "McDermott","Tetrade"), nlevel=NULL, base = 1) {

    if (length(n) < 2) stop("less than 2 groups")
    type <- match.arg(type)
    if (any(n < 2)) stop("less than 2 observations in at least one group")
    k <- length(n)
    if (base < 1 || base > k) stop("base is not between 1 and ", k)
    CM <- c()
    rnames <- c()
    if (!is.null(names(n)))
        varnames <- names(n)
    else 
        varnames <- 1:length(n)

    kindx <- 1:k
    if (base != 1 && type == "Dunnett") {
      n <- c(n[base], n[-base])
      varnames <- c(varnames[base], varnames[-base])
      kindx <- c(base, (1:k)[-base])
    }

    type <- match.arg(type)

    switch(type, "Dunnett" = {
        for(i in 2:k)
            CM <- rbind(CM, as.numeric(kindx==i)-as.numeric(kindx==1))
        rnames <- paste(varnames[2:k], "-", varnames[1], sep="")
    }, "Tukey" = {
        for (i in 1:(k-1)) {
            for(j in (i+1):k) {
                CM  <- rbind(CM, as.numeric(kindx==j)-as.numeric(kindx==i))
                rnames <- c(rnames, paste(varnames[j], "-", varnames[i],
                                          sep=""))
            }
        }
    }, "Sequen" =  {
        for (i in 2:k) {
            CM  <- rbind(CM, as.numeric(kindx==i)-as.numeric(kindx==i-1))
            rnames <- c(rnames, paste(varnames[i], "-", varnames[i-1],
                                      sep=""))
        }
    }, "AVE" = {
        help <- c(1,  -n[2:k]/sum(n[2:k]))
        CM <- rbind(CM, help)
        for (i in 2:(k-1)) {
            x <- sum(n[1:(i-1)])+sum(n[(i+1):k])
            help <- c(-n[1:(i-1)]/x, 1, -n[(i+1):k]/x)
            CM <- rbind(CM, help)
        }
        help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
        CM  <- rbind(CM, help)
        rnames <- paste("C", 1:nrow(CM))
    }, "Changepoint" = {
        for (i in 1:(k-1)) {
            help <- c(-n[1:i]/sum(n[1:i]), n[(i+1):k]/sum(n[(i+1):k]))
            CM <- rbind(CM, help)
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM), sep=""))
    }, "Williams" = {
        for (i in 1:(k-2)) {
            help <-  c(-1, rep(0, k-i-1), n[(k-i+1):k]/sum(n[(k-i+1):k]))
            CM <- rbind(CM, help)
        }
        help <- c(-1, n[2:k]/sum(n[2:k]))
        CM <- rbind(CM, help)
        rnames <- c(rnames, paste("C", 1:nrow(CM), sep=""))
    }, "Marcus" = {
        cm1 <- matrix(0, nrow=k-1, ncol=k)
        cm2 <- cm1
        for (i in 1:(k-1)) {
            cm1[i,(i+1):k] <- n[(i+1):k]/sum(n[(i+1):k])
            cm2[i,1:i] <- n[1:i]/sum(n[1:i])
        }
        row <- k*(k-1)/2
        index <- 1
        for (i in 1:(k-1)) {
            for (j in 1:i) {
                help <- cm1[i,]-cm2[j,]
                CM <- rbind(CM, help)
                index <- index+1
            }
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM), sep=""))
     }, "McDermott" = {
         for(i in 1:(k-2)) {
             help  <- c(-n[1:i]/sum(n[1:i]), 1, rep(0, k-i-1))
             CM <- rbind(CM, help)
         }
         help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
         CM  <- rbind(CM, help)
         rnames <- c(rnames, paste("C", 1:nrow(CM), sep=""))
    }, "Tetrade" = {
        if (is.null(nlevel)) stop("nlevel missing")
        if (length(nlevel) != 2) stop("only two factors allowed")
        a <- nlevel[1]
        b <- nlevel[2]
	idi <- 1:a
	idj <- 1:b
        for (i1 in 1:(a-1)) {
            for (i2 in (i1+1):a) {
	        for (j1 in 1:(b-1)) {
        	    for (j2 in (j1+1):b) {
                	CM <- rbind(CM, kronecker( ( as.numeric(idi==i1)-as.numeric(idi==i2) ),
                                                   ( as.numeric(idj==j1)-as.numeric(idj==j2) ) ) ) 
		        rnames <- c(rnames, paste( "(", i1, j1, "-", i1, j2, ")", "-", 
                                                   "(", i2, j1, "-", i2, j2, ")",  sep=""))
            	    }
        	}
	    }
        }
    },)
    rownames(CM) <- rnames
    if (type=="Tetrade")
      colnames(CM) <- NULL
    else 
      colnames(CM) <- varnames
    CM
}

contr.Dunnett <- function(n, base = 1, contrasts=TRUE) {
  if (length(n) == 1)  {
    if (contrasts) {
      mginv(contrMat(rep(10, n), base = base, type="Dunnett"))
    } else {
      diag(n)
    }
  }
  if (length(n) > 1) { 
    if (contrasts) {
      x <- rep(10, length(n))
      names(x) <- n
      mginv(contrMat(x, base = base, type="Dunnett")) 
    } else {
      diag(length(n))
    }
  } 
}

contr.Tukey <- function(n, contrasts=TRUE) {
  if (!contrasts) stop("contrasts is false")
  if (length(n) == 1) 
    mginv(contrMat(rep(10, n), type="Tukey"))
  if (length(n) > 1) {
    x <- rep(10, length(n))
    names(x) <- n
    mginv(contrMat(x, type="Tukey"))
  }
}
# $Id: internals.R,v 1.8 2003/05/12 15:23:20 hothorn Exp $

ct <- function(x, contrasts=FALSE) {
  a <- diag(length(x))
  colnames(a) <- x
  a
}

getdigits <- function(x) {
  if (x > 0.1 || x <= 0) return(NA)
  if (1/x >= 10 && 1/x < 100) return(1)
  if (1/x >= 100 && 1/x < 1000) return(2)
  if (1/x >= 1000 && 1/x < 10000) return(3)
  if (1/x >= 10000) return(4)
}

nicepaste <- function(x, pattern) {
  if (length(x) == 1) RET <- x
  if (length(x) == 2) RET <- paste(x[1], pattern, x[2], collapse="", sep=" ")
  if (length(x) > 2)
    RET <- paste(paste(x[1:(length(x)-1)], pattern, collapse="", sep=" "),
                       x[length(x)], collapse="", sep=" ")
  RET
}
# Copyright (C) 1994-9  W. N. Venables and B. D. Ripley.
# Copied from VR_6.3-2

mginv <-  function(X, tol = sqrt(.Machine$double.eps))
     {
     ## Generalized Inverse of a Matrix
       dnx <- dimnames(X)
       if(is.null(dnx)) dnx <- vector("list", 2)
       s <- svd(X)
       nz <- s$d > tol * s$d[1]
       structure(
         if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
         dimnames = dnx[2:1])
     }
# $Id: parseformula.R,v 1.16 2003/08/26 11:52:14 hothorn Exp $

parseformula <- function(formula, data=list(), subset, na.action, 
                         whichf=NULL, ...) {
    if(missing(na.action))
        na.action <- getOption("na.action")

    # get additional arguments (type may be needed here)
    cargs <- list(...)

    # the model terms
    mt <- terms(formula, data=data)
    
    # check for intercept and interactions, Tetrade or specific contrasts
    INTERCEPT <- (attr(mt, "intercept") == 1)
    INTERACTIONS <- any(attr(mt, "order") > 1)
    TETRADE <- FALSE
    if (!is.null(cargs$type))
      TETRADE <- (cargs$type == "Tetrade")
    CMATRIX <- FALSE
    if (!is.null(cargs$cmatrix))
      CMATRIX <- TRUE

    # get the model.frame
    m <- match.call(expand.dots = FALSE)
    m$whichf <- NULL
    if(missing(subset)) 
        m$subset <- NULL
    else 
        m$subset <- subset
    if(missing(na.action)) 
        m$na.action <- NULL
    else  
        m$na.action <- na.action
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())

    # all variables used
    xvars <- as.character(attr(mt, "variables"))[-1]

    # variables on the right hand side only
    if ((yvar <- attr(mt, "response")) > 0) {
      resp <-  xvars[yvar]
      xvars <- xvars[-yvar]
    } else {
      stop("no response specified")
    }

    # extract the response
    response <- attr(attr(mf, "terms"), "response")
    y <- mf[[response]]

    # handle offsets
    offset <- model.offset(mf)
    if(!is.null(offset) && length(offset) != length(y)) {
      stop("Number of offsets is ", length(offset),
           ", should equal ", length(y), " (number of observations)")
    }

    if (!is.null(offset)) {
      y <- y - offset
    }

    ncovar <- c(0,0)
    nlevel <- NULL
    nobs <- NULL

    # get the terms on the right hand side as in the call
    rhs <- attr(mt, "term.labels")

    mainF <- NULL
    coVar <- NULL

    # determine, which variables are factors
    fact <- unlist(lapply(mf[xvars], is.factor))

    # remove unused factor levels (e.g. induced by `subset')
    mf[xvars[fact]] <- lapply(mf[xvars[fact]], factor)

    # check if any matches whichf, if given
    if (!missing(whichf)) {
      if(!any(xvars[fact] == whichf)) {
        err <- paste("no variable", whichf, "given in formula")
        stop(err)
      }
    }

    # sorry, at least one factor required
    if (sum(fact) == 0) stop("at least one factor required")

    nlevel <- unlist(lapply(mf[xvars[fact]], nlevels))

    flevels <- lapply(mf[xvars[fact]], levels)

    # everything clear
    if (sum(fact) == 1) {
      mainF <- xvars[fact]
      coVar <- rhs[rhs != mainF]
    }

    # take a deeper look
    if (sum(fact) > 1) {
      # ok, the user tells us what to do
      if (!missing(whichf)) {
        fact[xvars != whichf] <- FALSE
        mainF <- xvars[fact]
        coVar <- rhs[rhs != mainF]
      } else {
        # cannot determine which one to choose
        if (!INTERACTIONS) {
          # ok, the user gets want he wants: build a design matrix using the
          # whole right hand side of formula (we check the dimensions later)
          if (CMATRIX)
            mainF <- rhs
          else 
            stop("more than one factor specified but whichf not given!")
        } else {
          # factors NOT involved in an interaction term
          fnames <- xvars[fact]
          fnames <- fnames[fnames %in% attr(mt, "term.labels")]
          # some procedure as above
          if (length(fnames) == 1) {
            mainF <- xvars[xvars == fnames]
            coVar <- rhs[rhs != mainF]
          } 
          if (length(fnames) > 1) {
            if (!missing(whichf)) {
              fact[xvars != whichf] <- FALSE
              mainF <- xvars[xvars == whichf]
              coVar <- rhs[rhs != mainF]
            } else {
              # ok, the user gets want he wants!
              if (CMATRIX)
                mainF <- rhs
              else 
                stop("more than one factor specified but whichf not given!")
            }
          }
          if (length(fnames) == 0) {
            # ok, everything is in interactions
            # this is allowed for TETRADE or CMATRIX 
            if (CMATRIX) {
              mainF <- unlist(strsplit(rhs[!rhs %in% xvars], ":"))
              coVar <- rhs[rhs != rhs[!rhs %in% xvars]]
            } else {
              if (TETRADE) {
                if(sum(fact) == 2) {
                  mainF <- unlist(strsplit(rhs[!rhs %in% xvars], ":"))
                  coVar <- rhs[rhs != rhs[!rhs %in% xvars]]
                } else {
                  stop("Tetrade assumes 2 factors in one interaction term")
                }
              } else {
                stop("use Tetrade or cmatrix for interactions")
              }
            }
          }                
        }
      } 
    } 

    ncovar <- c(0,0)


    ### <FIXME> check if contrasts to mainF have been previously defined
    ### via C() or contrasts()
    ### </FIXME>

    ### interaction contrasts
    if (length(mainF) > 1) mainF <- paste(mainF, collapse = ":")

    if (!is.null(mainF)) {
      eval(parse(text=paste("thisctrs <- list(", 
                         paste(mainF, "=\"ct\"", collapse=", "), ")")))
      x <- model.matrix(attr(mf, "terms"), mf, contrasts=thisctrs)

      fm <- attr(mt, "factors")
      ### determine which columns of the design matrix 
      ### correspond to the factor of interest.
      indx <- which(colnames(fm) == mainF)
      mainFindx <- which(attr(x, "assign") == indx)

      ### ???
      intera <- grep(":", colnames(x)[mainFindx])
      if (length(intera) > 0 & (length(intera) < length(mainFindx))) {
        mainFindx <- mainFindx[-intera]
      }

      nobs <- apply(x[,mainFindx], 2, sum)

      # get number of covariables
      ncovar[1] <- mainFindx[1] - 1
      ncovar[2] <- ncol(x) - max(mainFindx)

      # get the number of observations at each factor level
      # NOT needed for Tetrade in contrMatr
      nobs <- apply(x[,mainFindx], 2, sum)

      if (TETRADE) {
        if (INTERCEPT) {
          tx <- x[,-1]
        } else {
          tx <- x
        }
        col <- 0
        newx <- c()
        for (i in 1:ncol(tx))
          newx <- rbind(newx, tx[tx[,i] == 1,])
        tx <- newx
        for (i in 1:nlevel[2]) {
          for (j in 1:nlevel[1]) {
            col <- col + 1
            newx[, i + nlevel[2]*(j-1)] <- tx[,col]
          }
        }
        if (INTERCEPT) {
            x <- cbind(x[,1],newx)
        } else {
            x <- newx
        }
      }
      xall <- x
    } else {
      stop("Could not find main effect")
    }

    xall <- x
    x <- x[,mainFindx]

    if (CMATRIX) {
      # sanity checks here: we need to reject contrast matrices of 
      # not matching dimensions
      if (is.matrix(cargs$cmatrix)) {
        if (ncol(cargs$cmatrix) != ncol(xall) && ncol(cargs$cmatrix) != ncol(x)) {
          stop("dimensions of x and cmatrix do not match")
        } else {
          # the dimensions of cmatrix match the dimension of the complete
          # design matrix, therefore information on the number of
          # covariables obsolete

          # <FIXME> this may be completly wrong since the variables
          # may have been reordered
          if (ncol(cargs$cmatrix) == ncol(xall)) {
            colnames(cargs$cmatrix) <- colnames(xall)
            ncovar <- c(0,0)
          }
          # </FIXME>
        }
      }
      if (is.vector(cargs$cmatrix)) {
        if (length(cargs$cmatrix) != ncol(xall) 
            && length(cargs$cmatrix) != ncol(x) ) {
          stop("dimensions of x and cmatrix do not match")
        } else {
          if (length(cargs$cmatrix) == ncol(xall))
              ncovar <- c(0,0)
        }
      }
    }

    cargs$nlevel <- nlevel
    cargs$nzerocol <- ncovar
    cargs$nobs <- nobs
    if (length(coVar) > 1) coVar <- nicepaste(coVar, "+")

    fnames <- list(response = resp, mainF = mainF, coVar=coVar)

    RET <- list(y=y, x=xall, cargs=cargs, fnames=fnames)
    RET
}
# $Id: plot.hmtest.R,v 1.4 2002/08/15 07:28:05 hothorn Exp $

plot.hmtest <- function(x, ltycint=2, ltyzero=3, ...) {
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
  if (ONE)
    xlab <- paste(format(100 * conf.level), "\%", 
                 "one-sided", "confidence intervals")
  else
    xlab <- paste(format(100 * conf.level), "\%",
                 "two-sided", "confidence intervals")     


  # <FIXME>
  # strwidth needs and open graphic device.
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
  plot(rbind(c(crange[1], 1), c(crange[2], n)), type="n", axes=FALSE,
       xlab=xlab, ylab="", main=type, ...)
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
# $Id: print.hmtest.R,v 1.7 2003/05/13 10:29:13 hothorn Exp $

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
            colnames(ecout) <- c("Estimate", "lower CI", "upper CI")
            print(ecout)
        }
    }
}
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
# $Id: simint.R,v 1.45 2003/07/03 09:13:21 hothorn Exp $

simint <- function(y, ...) UseMethod("simint")

simint.default <- function(y, x=NULL, type=c("Dunnett", "Tukey",
                     "Sequen", "AVE", "Changepoint", "Williams", "Marcus",
                     "McDermott","Tetrade"), cmatrix=NULL, conf.level=0.95,
                     alternative=c("two.sided","less","greater"), asympt=FALSE,
                     eps=0.001, maxpts=1e+06, nlevel=NULL, nzerocol=c(0,0),
                     ...)
{
    ctype <- match.arg(type)

    addargs <- list(...)
    base <- 1
    if (!is.null(addargs$base))
      base <- addargs$base

    # Compute the parameter estimates and their covariance

    xpxi   <- mginv(t(x) %*% x)
    rankx  <- sum(diag((xpxi %*% (t(x) %*% x))))
    n      <- nrow(x)
    p      <- ncol(x)
    df     <- round(n-rankx)
    estpar <- xpxi %*% t(x) %*% y
    mse    <- t(y-x %*% estpar) %*% (y-x %*% estpar)/df
    covm   <- mse[1,1]*xpxi

    # compute the appropriate contrast matrix, if not given

    if (!is.null(cmatrix)) ctype <- "user-defined"

    if (is.null(cmatrix)) {
        nobs <- apply(x[,(nzerocol[1] + 1):(ncol(x) - nzerocol[2])], 2, sum)
        cm <- contrMat(n = nobs, type = ctype, nlevel = nlevel, base = base)
    } else {
        cm <- cmatrix
    }

    if (nzerocol[1] != 0)
      cm <- cbind(matrix(0, ncol=nzerocol[1], nrow=nrow(cm)), cm)
    if (nzerocol[2] != 0)
      cm <- cbind(cm, matrix(0, ncol=nzerocol[2], nrow=nrow(cm)))

    ### test if cm %*% estpar is estimable

    a <- cm %*%  xpxi %*% (t(x) %*% x)
    b <- cm
    attributes(a) <- NULL
    attributes(b) <- NULL
    areequal <- all.equal(a, b)

    if (!is.logical(areequal))
        warning("at least one contrast is not estimable")

    csimint(estpar, df, covm, cm, ctype, conf.level, alternative,
            asympt, eps, maxpts)
}

csimint <- function(estpar, df, covm, cmatrix=NULL, ctype="user-defined",
                    conf.level=0.95,
                    alternative=c("two.sided","less","greater"), asympt=FALSE,
                    eps=0.001, maxpts=1000000)
{
    if (!is.vector(estpar) & !is.matrix(estpar)) stop("estpar not a vector")
    p <- length(estpar)
    if (missing(df) & !asympt) {
      stop("df is missing")
    } else {
      if (missing(df)) df <- 0
      if (!all.equal(df - floor(df), 0)) stop("df not an integer")
    }
    if (!is.matrix(covm)) {
      if (length(covm) == 1) 
        covm <- as.matrix(covm)
      else 
        stop("covm is not a matrix")
    }
    cm <- cmatrix
    if (ctype !="user-defined") cmatrix <- NULL

    alternative <- match.arg(alternative)

    if (asympt) df <- 0                          
 
    covm  <- cm %*% covm %*% t(cm)
    d     <- 1/sqrt(diag(covm))                            
    if (length(d) > 1)
      d <- diag(d)              
    cr    <- d %*% covm %*% d

    ests  <- cm %*% estpar
    ses   <- sqrt(diag(covm))
    tvals <- ests/ses
    dim   <- ncol(cr)

    # compute the p-values

    pfct <- function(q, conf=FALSE) {
        switch(alternative, "two.sided" = {
                  low <- rep(-abs(q), dim)
                  upp <- rep( abs(q), dim)
           }, "less" = {
                  low <- rep(-abs(q), dim)
                  upp <- rep(    Inf, dim)
           }, "greater" = {
                  low <- rep(   -Inf, dim)
                  upp <- rep( abs(q), dim)
           },)
           if (conf)
               pmvt(lower=low, upper=upp, df=df, corr=cr,
                    abseps=eps/10, maxpts=maxpts)-conf.level
           else 
               pmvt(lower=low, upper=upp, df=df, corr=cr,
                    abseps=eps, maxpts=maxpts)
    }

    switch(alternative, "two.sided" = {
        if (df>0) rawp <- 2*(1-pt(abs(tvals),df))     
        else      rawp <- 2*(1-pnorm(abs(tvals)))
    }, "less" = {
        if (df>0) rawp <- pt(tvals,df) 
        else      rawp <- pnorm(tvals)
    }, "greater" = {
       if (df>0) rawp <- 1-pt(tvals,df)
       else      rawp <- 1-pnorm(tvals)
    },)

    adjp <- 1-apply(tvals, 1, pfct)
    bonp <- pmin(1,dim*rawp)

   # and the simultaneous confidence intervals

    calpha <- uniroot(pfct, lower=0, upper=5, tol=eps, conf=TRUE)$root    

    switch(alternative, "two.sided" = {  
        LowerCL <- ests - calpha*ses
        UpperCL <- ests + calpha*ses
    }, "less" = {
        LowerCL <- rep(-Inf, dim)
        UpperCL <- ests + calpha*ses
    }, "greater" = {
        LowerCL <- ests - calpha*ses
        UpperCL <- rep( Inf, dim)
    },)

    cint <- cbind(LowerCL, UpperCL)
    colnames(cint) <- c("lower", "upper")
    attr(cint, "conf.level") <- conf.level

    RET <- list(cmatrix = cm, ctype = ifelse(is.null(cmatrix), ctype, NA), 
                estimate = ests, sd = ses, statistics = tvals,
                p.value.raw = rawp, p.value.bon = bonp,
                p.value.adj = adjp, conf.int = cint, eps=eps, calpha=calpha,
                asympt = asympt)
    class(RET) <- "hmtest"
    RET
}

simint.formula <-
function(formula, data=list(), subset, na.action, whichf, ...)
{
    cl <- match.call()
    if (!missing(subset)) subset <- cl$subset
    if (!missing(na.action)) na.action <- cl$na.action

    pf <- parseformula(formula, data, subset, na.action, whichf, ...)
    x <- pf$x
    y <- pf$y
    cargs <- pf$cargs

    attr(x, "contrasts") <- NULL
    attr(x, "assign") <- NULL
    y <- do.call("simint", c(list(y=y, x=x), cargs))
    y$DNAME <- cl
    y$FNAMES <- pf$fnames
    y
}

simint.lm <- function(y, psubset=NULL, conf.level=0.95, cmatrix = NULL, 
                      alternative=c("two.sided","less","greater"), 
                      asympt=FALSE, eps=0.001, maxpts=1000000, ...) {

  p <- length(coef(y))
  pnames <- names(coef(y))
  ctype <- "model"
  if (p < 2) return(y)

  if (!is.null(psubset)) {
    switch(class(psubset), "integer" = {
      if (!all(psubset %in% 1:p))
        stop("incorrect psubset")
    }, "numeric" = { 
      if (!all(psubset %in% 1:p))
        stop("incorrect psubset")
    }, "character" = {
      if (all(psubset %in% names(coef(y)))) {
        psubset <- which(names(coef(y)) %in% psubset)
      } else {
        stop("incorrect psubset")
      }
    }, {
      stop("psubset is neither of class integer, numeric nor character")
    })
  } else {
    psubset <- 1:p
  }

  estpar <- coef(y)[psubset]
  df <- summary(y)$df[2]    
  covm = vcov(y)[psubset, psubset]
  pnames <- pnames[psubset]
  p <- length(psubset)
  
  if (!is.null(cmatrix)) {
    if (!is.matrix(cmatrix))
      stop("argument cmatrix is no a matrix")
    if (ncol(cmatrix) != p)
        stop("argument cmatrix does not have", p, " columns")
    if (is.null(rownames(cmatrix)))
      rownames(cmatrix) <- paste("C", 1:nrow(cmatrix), sep="")
  } else {
    cmatrix <- diag(p)
  }

  RET <- csimint(estpar = estpar, df = df, covm = covm, 
           cmatrix = cmatrix, conf.level = conf.level, 
           alternative = alternative,
           eps = eps, maxpts = maxpts)
  if (!is.null(pnames) & is.null(rownames(RET$estimate))) 
    rownames(RET$estimate) <- pnames
  RET$ctype <- ctype
  return(RET)
}

simint.glm <- function(y, psubset=NULL, conf.level=0.95, cmatrix = NULL, 
                      alternative=c("two.sided","less","greater"), 
                      asympt=FALSE, eps=0.001, maxpts=1000000, ...) {

  p <- length(coef(y))
  pnames <- names(coef(y))
  ctype <- "model"
  if (p < 2) return(y)

  if (asympt) {
    df <- summary(y)$df[2]
    if (y$family$family != "gaussian" || y$family$link != "identity") {
      warning("cannot compute confidence intervals based on t distribtion")
      asympt <- FALSE
    }
  }  
  if (asympt) df <- 0 else df <- summary(y)$df[2]


  if (!is.null(psubset)) {
    switch(class(psubset), "integer" = {
      if (!all(psubset %in% 1:p))
        stop("incorrect psubset")
    }, "numeric" = { 
      if (!all(psubset %in% 1:p))
        stop("incorrect psubset")
    }, "character" = {
      if (all(psubset %in% names(coef(y)))) {
        psubset <- which(names(coef(y)) %in% psubset)
      } else {
        stop("incorrect psubset")
      }
    }, {
      stop("psubset is neither of class integer, numeric nor character")
    })
  } else {
    psubset <- 1:p
  }

  estpar <- coef(y)[psubset]
  covm = vcov(y)[psubset, psubset]
  pnames <- pnames[psubset]
  p <- length(psubset)
  
  if (!is.null(cmatrix)) {
    if (!is.matrix(cmatrix))
      stop("argument cmatrix is no a matrix")
    if (ncol(cmatrix) != p)
        stop("argument cmatrix does not have", p, " columns")
    if (is.null(rownames(cmatrix)))
      rownames(cmatrix) <- paste("C", 1:nrow(cmatrix), sep="")
  } else {
    cmatrix <- diag(p)
  }

  RET <- csimint(estpar = estpar, df = df, covm = covm, 
           cmatrix = cmatrix, conf.level = conf.level, 
           alternative = alternative,
           eps = eps, maxpts = maxpts)
  if (!is.null(pnames) & is.null(rownames(RET$estimate))) 
    rownames(RET$estimate) <- pnames
  RET$ctype <- ctype
  return(RET)
}

# $Id: simtest.R,v 1.56 2003/08/13 10:14:07 bretz Exp $

simtest <- function(y, ...) UseMethod("simtest")

simtest.default <- function(y, x=NULL, type=c("Dunnett", "Tukey",
                     "Sequen", "AVE", "Changepoint", "Williams", "Marcus",
                     "McDermott","Tetrade"), cmatrix=NULL,
                     alternative=c("two.sided","less", "greater"), asympt=FALSE,
                     ttype=c("free","logical"), eps=0.001, maxpts=1e+06,
                     nlevel=NULL, nzerocol=c(0,0), ...)

{
    ctype <- match.arg(type)
    alternative <- match.arg(alternative)

    addargs <- list(...)
    base <- 1
    if (!is.null(addargs$base))
      base <- addargs$base


    xpxi   <- mginv(t(x) %*% x)
    rankx  <- sum(diag((xpxi %*% (t(x) %*% x))))
    n      <- nrow(x)
    p      <- ncol(x)
    df     <- round(n-rankx)
    estpar <- xpxi %*% t(x) %*% y
    mse    <- t(y-x %*% estpar) %*% (y-x %*% estpar)/df
    covm   <- mse[1,1]*xpxi

    # compute the appropriate contrast matrix, if not given

    if (!is.null(cmatrix)) ctype <- "user-defined"
   
    if (is.null(cmatrix)) { 
      nobs <- apply(x[,(nzerocol[1] + 1):(ncol(x) - nzerocol[2])], 2, sum) 
      cm <- contrMat(n = nobs, type = ctype, nlevel = nlevel, base = base)
      if (alternative == "greater")
        cm <- -cm
    } else {
      cm <- cmatrix
    }

    if (nzerocol[1] != 0)
      cm <- cbind(matrix(0, ncol=nzerocol[1], nrow=nrow(cm)), cm)
    if (nzerocol[2] != 0)
      cm <- cbind(cm, matrix(0, ncol=nzerocol[2], nrow=nrow(cm)))

    ### test if cm %*% estpar is estimable

    a <- cm %*%  xpxi %*% (t(x) %*% x)
    b <- cm
    attributes(a) <- NULL
    attributes(b) <- NULL
    areequal <- all.equal(a, b)

    if (!is.logical(areequal)) 
        warning("at least one contrast is not estimable")



    csimtest(estpar, df, covm, cm, ctype, ttype,
         alternative,  asympt, eps, maxpts)
}


csimtest <- function(estpar, df, covm, cmatrix=NULL, ctype="user-defined",
                    ttype=c("free","logical"), alternative=c("two.sided",
                    "less","greater"), asympt=FALSE,
                    eps=0.001, maxpts=1000000)
{
    if (!is.vector(estpar) & !is.matrix(estpar)) stop("estpar not a vector")
    p <- length(estpar)
    if (missing(df) & !asympt) {
      stop("df is missing")
    } else {
      if (missing(df)) df <- 0
      if (!all.equal(df - floor(df), 0)) stop("df not an integer")
    }
    if (!is.matrix(covm)) {
      if (length(covm) == 1) 
        covm <- as.matrix(covm)
      else 
        stop("covm is not a matrix")
    }
    cm <- cmatrix
    if (ctype !="user-defined") cmatrix <- NULL

    alternative <- match.arg(alternative)
    ttype <- match.arg(ttype)

    if (asympt) df <- 0                          
 
    covm  <- cm %*% covm %*% t(cm)                            
    d     <- 1/sqrt(diag(covm))                            
    if (length(d) > 1)
      d <- diag(d)              
    cr    <- d %*% covm %*% d

    ests  <- cm %*% estpar
    ses   <- sqrt(diag(covm))
    tvals <- ests/ses
    dim   <- ncol(cr)
    delta <- rep(0,dim)

    switch(alternative, "two.sided" = {
	tvals <- -abs(tvals)
        if (df>0) rawp <- 2*pt(tvals,df)     
        else      rawp <- 2*pnorm(tvals)
    }, "less" = {
        if (df>0) rawp <- pt(tvals,df) 
        else      rawp <- pnorm(tvals)
    }, "greater" = {
       if (df>0) rawp <- pt(tvals,df)
       else      rawp <- pnorm(tvals)
    },)


    pvals    <- as.matrix(rawp)
    tvals    <- as.matrix(tvals)    
    covcont  <- covm

    k        <- nrow(cm)
    ### <FIXME>
    ### we need to check if the p-values are unique
#    if (any(duplicated(pvals))) {
#        warning("duplicated pvalues")
#    }
    ### </FIXME>
#    r        <- t(rank(pvals))
    ordpval1 <- order(pvals)
    ordpval2 <- (1:length(ordpval1))[order(ordpval1)]
    r        <- t(ordpval2)
# end fix
    ir       <- r
    ir[,r]   <- 1:nrow(pvals)
    origord  <- t(ir)         
    cord     <- cm[ir,]
    tvalsord <- tvals[ir,]
    pvalsord <- pvals[ir,]
    ccord    <- covcont[ir,ir]
    crrccord <- solve(sqrt(diag(diag(ccord)))) %*% ccord %*% solve(sqrt(diag(diag(ccord))))
    cct      <- t(cord)
   
    if (ttype == "logical" & k > 2) {
        for (iout in 1:(k-2)) {
	    limit <- 2^(k-iout-1)
	    in2   <- rep(0, k-iout-1)
	    zero  <- as.matrix(rep(0, k))
	    in1   <- as.matrix(zero)
	    y     <- cct[,1:iout]

            for (kk in 1:limit) {
                if (kk == limit) {
                    in2 <- rep(0, k-iout-1)
		}
		else { 
		    ii <- 1
	            zz <- kk
                    while( zz %% 2 == 0) {    
		        ii<-ii+1
                        zz<-zz/2
                    }
                    if (in2[ii] == 0) 
                        in2[ii] <- 1
                    else 
                        in2[ii] <- 0 
                }

                locbin  <- as.matrix( rbind( as.matrix(rep(0, iout)), 
                                             as.matrix(rbind(1, as.matrix(in2) ) ) ) )
                loc1 <- 1
                for (jj in 1:nrow(locbin)) {
                    if (locbin[jj] == 0) {
                    }
                    else {
                        loc1 <- c(loc1, jj)
                    }
                }
                loc1 <- t(loc1)
                loc1 <- t(loc1[2:ncol(loc1)])
                x    <- as.matrix(cct[,loc1])
                res    <- y - x %*% mginv(t(x) %*% x) %*% t(x) %*% y
                ssemat <- diag(t(res) %*% res)
                if (all(ssemat > .00000001)) {
                    if (identical(in1,zero)) {
                        in1 <- locbin
                        }
                    else {
                        check <- in1 - rep(locbin, ncol(in1))
                        diff  <- apply(check,2,max) - apply(check,2,min)
                        if (min(diff) == 2)
                            in1 <- cbind(as.matrix(in1),as.matrix(locbin))
                        else {
                            mindx <- which.min(diff)
                            if (sum(check[,mindx]) == -1) 
                                in1[,mindx] <- locbin
                        }
                    }
                }
            }
            in1   <- t(in1)
            ncont <- nrow(in1)
            in1   <- cbind(as.matrix(rep(iout+1,ncont)),in1)
            if (iout == 1)
                inbig <- in1
            else
                inbig <- rbind(inbig, in1)
        }
        big <- rbind(t(rep(1,k+1)),inbig)   
    }
    

    if (ttype == "free" | k < 3)
         big <- t(rep(1,k+1))
    lastset <-  cbind( cbind( rep(k,1), t(rep(0,k-1)) ), rep(1,1) ) 
    big <- rbind(big, lastset)
    stepj <- big[,1]
    if (ttype == "free") 
        stepj <- 1:k
    SubsetK <- big[,2:ncol(big)]
    if (ttype == "free") {
        m <- t(rep(1,k))
        for (i in 2:k) {
            r <- cbind( t(rep(0,i-1)), t(rep(1,k-i+1)) )
            m <- rbind(m,r)
        }
        SubsetK <- m
    }
    nbig <- nrow(big)
    if (ttype == "logical") {
        yyy <- rnorm(nrow(big))
        aaa <- as.factor(big[,1])
    }
    else {
        yyy <- rnorm(nrow(as.matrix(stepj)))
        aaa <- as.factor(stepj)
    }   
    contrasts(aaa) <- "ct"
    ff <- (yyy ~ aaa)
    mf <- model.frame(ff)
    des <- model.matrix(ff, mf)
    des <- des[,2:ncol(des)]
    if (ttype == "logical")
        contonly <- big[,2:(k+1)]
    else 
        contonly <- SubsetK
    tcmpr <- des %*% tvalsord


    if (ttype == "free") 
        nbig <- k
    count   <- rep(0,nbig)
    countc  <- count
    countc2 <- count
    nnn <- ncol(crrccord)
    for (i in 1:nnn) {
        for (j in 1:nnn) {
            if ( (contonly[i,j] + des[i,j])  > 0 )
                contonly[i,j] <- 1
            else
                contonly[i,j] <- 0
        }
    }

    stepj   <- as.matrix(stepj)
    subsets <- cbind(contonly, stepj)
    gls     <- rep(0,nrow(contonly))
    stdgls  <- rep(0,nrow(contonly))

    for (i1 in 1:nrow(contonly)) {
        loct <- 1
        for (jj in 1:ncol(contonly)) {
            if (contonly[i1,jj] != 0)
                loct <- c(loct, jj)
        }
        loct <- t(loct)
        loct <- t(loct[2:ncol(loct)]) 
        cort  <- as.matrix(crrccord[loct,loct])
        nt    <- nrow(cort)
        low3  <-  tcmpr[i1]* rep(1,nt)
        upp3  <- -tcmpr[i1]* rep(1,nt)
        if (alternative != "two.sided")
           upp3 <- rep(Inf,nt)
#        maxpts <- 1000000
        delta  <- rep(0,nt) 
        prob       <- pmvt(lower=low3, upper=upp3, df=df, 
                           delta=delta, corr=cort, abseps=eps, maxpts=maxpts)
        gls[i1]    <- 1-prob
        stdgls[i1] <- attr(prob, "error")
    }

    glsbig <- des * ( as.matrix(gls) %*% t(rep(1,k)) )
    glsp   <- apply(glsbig,2,max)
    glsin  <- 1:ncol(glsbig)
    for (i in 1:ncol(glsbig)) {
        glsin[i]  <- which.max(glsbig[,i])
    }
    glsin  <- t(glsin)  
    stdgls <- t(stdgls)  
    stdgls  <- stdgls[glsin]

    for (i in 2:k) {
        if (glsp[i] < glsp[i-1]) {
            glsp[i]   <- glsp[i-1]
            stdgls[i] <- stdgls[i-1]
        }
    }

    seadjp  <- as.matrix(stdgls)

    # <FIXME>
    # is a warning appropriate?
    if (any(seadjp > eps)) warning("error > eps")
    # </FIXME>

    adjpgls <- as.matrix(glsp)
    adjp    <- as.matrix(adjpgls)

    if (ttype == "logical")
        totals  <-  apply(contonly,1,sum)
    else 
        totals <- k:1
    if (alternative == "two.sided") {
        if (df  == 0)
            bon <- 2*(pnorm(tcmpr)) * totals
        else 
            bon <- 2*(pt(tcmpr,df)) * totals
    }       
    else {
        if (df == 0) 
            bon <- (pnorm(tcmpr)) * totals
        else 
            bon <- (pt(tcmpr,df)) * totals
    }      

    bonbig  <- des * ( as.matrix(bon) %*% t(rep(1,k)) )
    bonp    <- apply(bonbig,2,max)
    bonmult <- bonp/pvalsord

    for (i in 2:k) {
        if (bonp[i] < bonp[i-1]) 
            bonp[i] <- bonp[i-1]
    }
    rawp     <- pvalsord
    ests <- ests[ir,]
    tvals<- tvals[ir,]
    if (alternative == "greater") {
        ests <- -ests
        cm <- -cm
        tvals <- -tvals
    }
    adjpbon  <- bonp
    adjpbon  <- pmin(1,adjpbon)
    rownames(adjp) <- rownames(cord)

    RET <- list(cmatrix = cm, ctype = ifelse(is.null(cmatrix), ctype, NA),
                estimate = as.matrix(ests), sd = ses, statistics = tvals,
                p.value.raw = rawp, p.value.bon = adjpbon,
                p.value.adj = adjp, eps=eps, asympt = asympt)

    class(RET) <- "hmtestp"
    RET
}

simtest.formula <-
function(formula, data=list(), subset, na.action, whichf, ...)
{
    cl <- match.call()
    if (!missing(subset)) subset <- cl$subset
    if (!missing(na.action)) na.action <- cl$na.action

    pf <- parseformula(formula, data, subset, na.action, whichf, ...)
    x <- pf$x
    y <- pf$y
    cargs <- pf$cargs

    attr(x, "contrasts") <- NULL
    attr(x, "assign") <- NULL
    y <- do.call("simtest", c(list(y=y, x=x), cargs))
    y$DNAME <- cl
    y$FNAMES <- pf$fnames
    y
}

simtest.lm <- function(y, psubset = NULL, cmatrix = NULL, 
                      ttype=c("free","logical"), 
                      alternative=c("two.sided","less","greater"), 
                      asympt=FALSE,
                      eps=0.001, maxpts=1000000, ...) {


  p <- length(coef(y))
  pnames <- names(coef(y))
  ctype <- "model"
  if (p < 2) return(y)

  if (!is.null(psubset)) {
    switch(class(psubset), "integer" = {
      if (!all(psubset %in% 1:p))
        stop("incorrect psubset")
    }, "numeric" = { 
      if (!all(psubset %in% 1:p))
        stop("incorrect psubset")
    }, "character" = {
      if (all(psubset %in% names(coef(y)))) {
        psubset <- which(names(coef(y)) %in% psubset)
      } else {
        stop("incorrect psubset")
      } 
    }, {
      stop("psubset is neither of class integer, numeric nor character")
    })
  } else {
    psubset <- 1:p
  }

  estpar <- coef(y)[psubset]
  df <- summary(y)$df[2]    
  covm = vcov(y)[psubset, psubset]
  pnames <- pnames[psubset]
  p <- length(psubset)
  
  if (!is.null(cmatrix)) {  
    if (!is.matrix(cmatrix))
      stop("argument cmatrix is no a matrix")
    if (ncol(cmatrix) != p)
        stop("argument cmatrix does not have", p, " columns")
    if (is.null(rownames(cmatrix)))
      rownames(cmatrix) <- paste("C", 1:nrow(cmatrix), sep="")
  } else {
    cmatrix <- diag(p)
  }

  RET <- csimtest(estpar = estpar, df = df, covm = covm,
           cmatrix = cmatrix,
           alternative = alternative, ttype = ttype,
           eps = eps, maxpts = maxpts)
  if (!is.null(pnames) & is.null(rownames(RET$estimate)))
    rownames(RET$estimate) <- pnames
  RET$ctype <- ctype
  return(RET)
}

simtest.glm <- function(y, psubset = NULL, cmatrix = NULL, 
                      ttype=c("free","logical"), 
                      alternative=c("two.sided","less","greater"), 
                      asympt=FALSE,
                      eps=0.001, maxpts=1000000, ...) {


  p <- length(coef(y))
  pnames <- names(coef(y))
  ctype <- "model"
  if (p < 2) return(y)

  if (asympt) {
    if (y$family$family != "gaussian" || y$family$link != "identity") {
      warning("cannot compute confidence intervals based on t distribtion")
      asympt <- FALSE
    }
  }  
     
  if (asympt) df <- 0 else df <- summary(y)$df[2]


  if (!is.null(psubset)) {
    switch(class(psubset), "integer" = {
      if (!all(psubset %in% 1:p))
        stop("incorrect psubset")
    }, "numeric" = { 
      if (!all(psubset %in% 1:p))
        stop("incorrect psubset")
    }, "character" = {
      if (all(psubset %in% names(coef(y)))) {
        psubset <- which(names(coef(y)) %in% psubset)
      } else {
        stop("incorrect psubset")
      } 
    }, {
      stop("psubset is neither of class integer, numeric nor character")
    })
  } else {
    psubset <- 1:p
  }

  estpar <- coef(y)[psubset]
  covm = vcov(y)[psubset, psubset]
  pnames <- pnames[psubset]
  p <- length(psubset)
  
  if (!is.null(cmatrix)) {  
    if (!is.matrix(cmatrix))
      stop("argument cmatrix is no a matrix")
    if (ncol(cmatrix) != p)
        stop("argument cmatrix does not have", p, " columns")
    if (is.null(rownames(cmatrix)))
      rownames(cmatrix) <- paste("C", 1:nrow(cmatrix), sep="")
  } else {
    cmatrix <- diag(p)
  }

  RET <- csimtest(estpar = estpar, df = df, covm = covm,
           cmatrix = cmatrix,
           alternative = alternative, ttype = ttype,
           eps = eps, maxpts = maxpts)
  if (!is.null(pnames) & is.null(rownames(RET$estimate)))
    rownames(RET$estimate) <- pnames
  RET$ctype <- ctype
  return(RET)
}


# $Id: summary.hmtest.R,v 1.10 2003/05/13 10:29:13 hothorn Exp $

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
            colnames(ecout) <- c("Estimate", "low CI,", "upp CI", "t value",
                                 "Std.Err.", "p raw", "p Bonf", "p adj")   
            print(ecout)
        }
    }
    invisible(x)
}




# $Id: summary.hmtestp.R,v 1.8 2003/05/13 10:29:13 hothorn Exp $

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
    if (x$asympt)
      cat("\t Asymptotic simultaneous tests:", type, "\n")
    else
      cat("\t Simultaneous tests:", type, "\n")
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




# $Id: zzz.R,v 1.4 2003/06/23 13:35:22 hothorn Exp $

.onLoad <- function(lib, pkg) {
    if(!require(mvtnorm))
        warning("Could not load package mvtnorm")
}
