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
    if(missing(subset)) m$subset <- NULL
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

    ncovar <- c(0,0)
    nlevel <- NULL
    nobs <- NULL

    # get the terms on the right hand side as in the call
    rhs <- attr(mt, "term.labels")

    mainF <- NULL
    coVar <- NULL

    # determine, which variables are factors
    fact <- unlist(lapply(mf[xvars], is.factor))

    # check if any matches whichf, if given
    if (!missing(whichf)) {
      if(!any(xvars[fact] == whichf)) {
        err <- paste("no variable", whichf, "given in formula")
        stop(err)
      }
    }

    # sorry, at least one factor requires
    if (sum(fact) == 0) stop("at least one factor required")

    nlevel <- unlist(lapply(mf[xvars[fact]], nlevels))

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
        # cannot determine, which one to choose
        if (!INTERACTIONS) {
          stop("more than one factor but whichf not given!")
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
              stop("more than one factor specified but whichf not given!")
            }
          }
          if (length(fnames) == 0) {
            # ok, everything is in interactions
            # this is allowed for TETRADE or CMATRIX 
            if (CMATRIX) {
              mainF <- rhs[!rhs %in% xvars]
              coVar <- rhs[rhs != mainF]
            } else {
              if (TETRADE) {
                if(sum(fact) == 2) {
                  mainF <- rhs[!rhs %in% xvars]
                  coVar <- rhs[rhs != mainF]
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

    if (!is.null(mainF)) {
      formulaEff <- as.formula(paste(deparse(formula[[2]], width=500), "~",
                                     mainF, "-1"))
      x <- model.matrix(terms(formulaEff, data=data), mf)
      if (TETRADE) {
        tx <- x
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
        x <- newx
      }
      xall <- x
    } else {
      stop("Could not find main effect")
    }

    # get the number of observations at each factor level
    # NOT needed for Tetrade in contrMatr
    nobs <- apply(x, 2, sum)

    xcov <- NULL
    ncovar <- c(0,0)
    if (!is.null(coVar)) {
      if (length(coVar) > 0) {
        coVar <- nicepaste(coVar, "+")
        formulaCov <- paste(deparse(formula[[2]], width=500), "~",
                            coVar, "-1")
        formulaCov <- as.formula(formulaCov)
        xcov <- model.matrix(terms(formulaCov, data=data), mf)
        xall <- cbind(x, xcov)
        ncovar[2] <- ncol(xcov)
      }
    }

    if (INTERCEPT) {
      xall <- cbind(1, xall)
      colnames(xall)[1] <- "(Intercept)"
      ncovar[1] <- 1
    }

    if (CMATRIX) {
      if (is.matrix(cargs$cmatrix)) {
        if (ncol(cargs$cmatrix) != ncol(xall) && ncol(cargs$cmatrix) != ncol(x)) {
          stop("dimensions of x and cmatrix do not match")
        } else {
          # the dimensions of cmatrix match the dimension of the complete
          # design matrix, therefore information on the number of
          # covariables obsolete
          if (ncol(cargs$cmatrix) == ncol(xall))
            ncovar <- c(0,0)
        }
      }
      if (is.vector(cargs$cmatrix)) {
        if (length(cargs$cmatrix) != ncol(xall) && length(cargs$cmatrix) != ncol(x) ) {
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

    fnames <- list(response = resp, mainF = mainF, coVar=coVar)

    RET <- list(y=y, x=xall, cargs=cargs, fnames=fnames)
    RET
}
