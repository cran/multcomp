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

