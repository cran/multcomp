# $Id: simint.R,v 1.28 2002/04/12 08:17:27 hothorn Exp $

simint <- function(y, ...) UseMethod("simint")

simint.default <- function(y, x=NULL, type=c("Dunnett", "Tukey",
                     "Sequen", "AVE", "Changepoint", "Williams", "Marcus",
                     "McDermott","Tetrade"), cmatrix=NULL, conf.level=0.95,
                     alternative=c("two.sided","less","greater"), asympt=FALSE,
                     eps=0.001, maxpts=1000000, nlevel=NULL, ...)
{
    ctype <- match.arg(type)

    # Compute the parameter estimates and their covariance

    xpxi   <- ginv(t(x) %*% x)
    rankx  <- sum(diag((xpxi %*% (t(x) %*% x))))
    n      <- nrow(x)
    p      <- ncol(x)
    df     <- n-rankx
    estpar <- xpxi %*% t(x) %*% y
    mse    <- t(y-x %*% estpar) %*% (y-x %*% estpar)/df
    covm   <- mse[1,1]*xpxi

    # compute the appropriate contrast matrix, if not given

    if (!is.null(cmatrix)) ctype <- "user-defined"

    nobs  <- apply(x, 2, sum)[2:p]	# omit Intercept
    if (is.null(cmatrix)) 
        cm <- contrMat(nobs, ctype, nlevel) 
    else cm <- cmatrix 

    csimint(estpar, as.integer(df), covm, cm, ctype, conf.level, alternative,
            asympt, eps, maxpts)
}

csimint <- function(estpar, df, covm, cmatrix=NULL, ctype="user-defined",
                    conf.level=0.95,
                    alternative=c("two.sided","less","greater"), asympt=FALSE,
                    eps=0.001, maxpts=1000000)
{
    if (!is.vector(estpar) & !is.matrix(estpar)) stop("estpar not a vector")
    p <- length(estpar)
    if (!is.integer(df)) stop("df not an integer")
    if (!is.matrix(covm)) stop("covm is not a matrix")
    cm <- cmatrix
    if (ctype !="user-defined") cmatrix <- NULL

    alternative <- match.arg(alternative)

    if (asympt) df <- 0                          
 
    covm  <- cm %*% covm %*% t(cm)                            
    d     <- diag(1/sqrt(diag(covm)))              
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

    RET <- list(cmatrix = cm[,2:p], ctype = ifelse(is.null(cmatrix), ctype, NA), 
                estimate = ests, sd = ses, statistics = tvals,
                p.value.raw = rawp, p.value.bon = bonp,
                p.value.adj = adjp, conf.int = cint, eps=eps, calpha=calpha)
    class(RET) <- "hmtest"
    RET
}

simint.formula <-
function(formula, data=list(), subset, na.action, ...)
{
    if(missing(na.action))
        na.action <- getOption("na.action")
    mt <- terms(formula, data=data)
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data) 
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    namD <- names(mf)
    dnames <- "Intercept"
    nlevel <- c()
    COVAR <- FALSE
    INTERACTIONS <- !is.null(grep(":", as.character(formula[3])))
    for(nn in namD[-1]) {
        if (is.factor(mf[[nn]])) {
            contrasts(mf[[nn]]) <- "ct"
            dnames <- c(dnames, levels(mf[[nn]]))
            nlevel <- c(nlevel, nlevels(mf[[nn]]))
        } else {
            # stop("no covariable allowed")
            COVAR <- TRUE
            dnames <- c(dnames, nn)
        }
    }
    cargs <- list(...)
    if (!is.null(nlevel)) cargs$nlevel <- nlevel
    response <- attr(attr(mf, "terms"), "response")
    y <- mf[[response]]
    x <- model.matrix(mt, mf)

#    <FIXME: Tetrade contrasts assume x-cols to be ordered in another way as
#   model.matrix returns> 
    col <- 0
    if (!is.null(cargs$type)) {
      if (cargs$type == "Tetrade") {
        tx <- x[,2:ncol(x)]
	col <- 0
        if (length(nlevel) == 2) {
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
          x <- cbind(1, newx)
        } else {
          stop("can use Tetrade contrasts with two factors only")
        }

      }
    }
#   </FIXME> 

    if (is.null(cargs$cmatrix)) {
      if (COVAR)             
        stop("cmatrix not given but covariables specified")
    } else { 
      if (is.matrix(cargs$cmatrix)) {
        if (ncol(cargs$cmatrix) != ncol(x)) 
          stop("dimensions of x and cmatrix do not match")
      }
      if (is.vector(cargs$cmatrix)) {
        if (length(cargs$cmatrix) != ncol(x))
          stop("dimensions of x and cmatrix do not match")
      }                                     
    }

    # <FIXME>: problem with interactions!!!
    if (!INTERACTIONS)
      colnames(x)[1:length(dnames)] <- dnames
    # </FIXME>
    attr(x, "contrasts") <- NULL
    attr(x, "assign") <- NULL
    y <- do.call("simint", c(list(y=y, x=x), cargs))
    if (length(namD) > 2) {
        y$DNAME <- paste(namD[1], "by", paste(namD[-1],
                 collapse = " + "))
    } else {
        y$DNAME <- paste(namD, collapse = " by ")
    }
    y
}
