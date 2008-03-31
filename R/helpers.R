

### oh dear!
### Cox models don't have any intercept ...
model.matrix.coxph <- function(object, ...) {
    mm <- model.matrix.default(object)
    at <- attributes(mm)
    mm <- mm[,-1]
    at$dim[2] <- at$dim[2] - 1
    at$dimnames[[2]] <- at$dimnames[[2]][-1]
    at$assign <- at$assign[-1]
    attributes(mm) <- at
    mm
}

model.matrix.aovlist <- function(object, ...)
    stop(sQuote("glht"), " does not support objects of class ", 
         sQuote("aovlist"))

model.matrix.lme <- function(object, ...)
    model.matrix(terms(object), data = model.frame(object), ...)

model.frame.lme <- function(object, ...)
    object$data

### some methods of (g)lmer objects
model.frame.mer <- function(object, ...) {
    x <- object@frame
    if (is.null(x))
        stop("models of class ", class(object), 
             " need to be fitted with argument ", 
             sQuote("model = TRUE"))
    x
}
model.frame.lmer <- model.frame.mer
model.frame.glmer <- model.frame.mer

model.matrix.mer <- function(object, ...) {
    return(model.matrix(terms(object), model.frame(object)))
}
model.matrix.lmer <- model.matrix.mer
model.matrix.glmer <- model.matrix.mer


### and now for (g)lmer2 (thanks for Manuel Eugster)
terms.lmer2 <- function(x, ...) return(x@terms) 
model.frame.lmer2 <- model.frame.mer
model.matrix.lmer2 <- function(object, ...) {
  return(model.matrix(terms(object), model.frame(object)))
}
terms.glmer2 <- function(x, ...) return(x@terms) 

model.frame.glmer2 <- model.frame.mer
model.matrix.glmer2 <- function(object, ...) {
  return(model.matrix(terms(object), model.frame(object)))
}


### extract coefficients, covariance matrix and 
### degrees of freedom (if available) from `model'
modelparm <- function(model, coef., vcov., df, ...) 
    UseMethod("modelparm")

modelparm.default <- function(model, coef. = coef, vcov. = vcov, 
                              df = NULL, ...) 
{

    ### extract coefficients and their covariance matrix
    beta <- try(coef.(model))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("model"), " found!")

    sigma <- try(vcov.(model))
    if (inherits(sigma, "try-error"))
        stop("no ", sQuote("vcov"), " method for ",
             sQuote("model"), " found!")       
    sigma <- as.matrix(sigma)

    if (any(length(beta) != dim(sigma))) 
        stop("dimensions of coefficients and covariance matrix don't match")

    ### determine degrees of freedom
    if (is.null(df)) {
        df <- 0
        ### check if a linear model was supplied
        if (class(model)[1] %in% c("aov", "lm")) {
            class(model) <- "lm"
            df <- summary(model)$df[2]
        }
        if (inherits(model, "parm"))
            df <- model$df
    } else {
        if (df < 0) stop(sQuote("df"), " is not positive")
    }

    ### try to identify non-estimable coefficients
    ### coef.aov removes NAs, thus touch coefficients 
    ### directly
    ocoef <- coef.(model)
    if (inherits(model, "aov")) ocoef <- model$coefficients
    estimable <- rep(TRUE, length(ocoef))
    if (any(is.na(ocoef))) {
        estimable[is.na(ocoef)] <- FALSE
        beta <- ocoef[estimable]
    }

    ### just in case...
    if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable))
        stop("could not extract coefficients and covariance matrix from ", 
             sQuote("model"))

    RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
    class(RET) <- "modelparm"
    RET
}

### mixed effects models (package `lme4')
modelparm.mer <- function(model, coef. = fixef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)
modelparm.lmer <- modelparm.mer
modelparm.glmer <- modelparm.mer
modelparm.lmer2 <- modelparm.mer
modelparm.glmer2 <- modelparm.mer


### extract coefficients, covariance matrix and 
### degrees of freedom (if available) from `model'
modelparm <- function(model, coef., vcov., df, ...) 
    UseMethod("modelparm")

modelparm.default <- function(model, coef. = coef, vcov. = vcov, 
                              df = NULL, ...) 
{

    ### extract coefficients and their covariance matrix
    beta <- try(coef.(model))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("model"), " found!")

    sigma <- try(vcov.(model))
    if (inherits(sigma, "try-error"))
        stop("no ", sQuote("vcov"), " method for ",
             sQuote("model"), " found!")       
    sigma <- as.matrix(sigma)

    if (any(length(beta) != dim(sigma))) 
        stop("dimensions of coefficients and covariance matrix don't match")

    ### determine degrees of freedom
    if (is.null(df)) {
        df <- 0
        ### check if a linear model was supplied
        if (class(model)[1] %in% c("aov", "lm")) {
            class(model) <- "lm"
            df <- summary(model)$df[2]
        }
    } else {
        if (df < 0) stop(sQuote("df"), " is not positive")
    }

    ### try to identify non-estimable coefficients
    ### coef.aov removes NAs, thus touch coefficients 
    ### directly
    ocoef <- coef.(model)
    if (inherits(model, "aov")) ocoef <- model$coefficients
    estimable <- rep(TRUE, length(ocoef))
    if (any(is.na(ocoef))) {
        estimable[is.na(ocoef)] <- FALSE
        beta <- ocoef[estimable]
    }

    ### just in case...
    if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable))
        stop("could not extract coefficients and covariance matrix from ", 
             sQuote("model"))

    RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
    class(RET) <- "modelparm"
    RET
}

### mixed effects models (package `lme4')
modelparm.mer <- function(model, coef. = fixef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)
modelparm.lmer <- modelparm.mer
modelparm.glmer <- modelparm.mer
modelparm.lmer2 <- modelparm.mer
modelparm.glmer2 <- modelparm.mer

### package `nlme'
modelparm.lme <- function(model, coef. = nlme:::fixef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### survreg models (package `survival')
vcovsurvreg <- function(object, ...) {
    sigma <- vcov(object)
    p <- length(coef(object))
    return(sigma[1:p, 1:p])
}

modelparm.survreg <- function(model, coef. = coef, vcov. = vcovsurvreg, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

modelparm.aovlist <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
    stop(sQuote("glht"), " does not support objects of class ", sQuote("aovlist"))


### modified from package MASS  
MPinv <- function (X, tol = sqrt(.Machine$double.eps))
{
    if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X)))
        stop("X must be a numeric or complex matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X))
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    if (all(Positive)) 
        RET <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))   
    else if (!any(Positive))
        RET <- array(0, dim(X)[2:1])
    else RET <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
        t(Xsvd$u[, Positive, drop = FALSE]))
    return(list(MPinv = RET, rank = sum(Positive)))
}

### meaningless ...
chkdots <- function(...) {

    lst <- list(...)
    if (length(lst) > 0) {
        warning("Argument(s) ", sQuote(names(lst)), " passed to ", sQuote("..."), 
                " are ignored", call. = TRUE)
    }
}
