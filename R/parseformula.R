# $Id: parseformula.R,v 1.17 2003/09/11 14:41:56 hothorn Exp $

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


    if (!is.null(mainF)) {
      eval(parse(text=paste("thisctrs <- list(", 
                         paste(mainF, "=\"ct\"", collapse=", "), ")")))
      x <- model.matrix(attr(mf, "terms"), mf, contrasts=thisctrs)

      ### interaction contrasts
      if (length(mainF) > 1) mainF <- paste(mainF, collapse = ":")

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
