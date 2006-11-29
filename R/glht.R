
### general linear hypotheses
glht <- function(model, linfct, ...) UseMethod("glht", linfct)

### K coef(model) _!alternative_ rhs
glht.matrix <- function(model, linfct, 
    alternative = c("two.sided", "less", "greater"), rhs = 0, ...) {

    ### extract coefficients and their covariance matrix, df
    mpar <- modelparm(model, ...)

    alternative <- match.arg(alternative)
    if (!is.numeric(rhs))
        stop(sQuote("rhs"), " is not a numeric vector")

    if (ncol(linfct) != length(mpar$coef))
        stop(sQuote("ncol(linfct)"), " is not equal to ", 
             sQuote("length(coef(model))"))

    if (is.null(colnames(linfct)))
        colnames(linfct) <- names(mpar$coef)

    if (is.null(rownames(linfct))) # {
        rownames(linfct) <- 1:nrow(linfct)
#    } else {
        ### alt <- switch(alternative, 
        ###    "two.sided" = "==", "less" = ">=", "greater" = "<=")
        ### rownames(linfct) <- paste(rownames(linfct), alt, rhs)
#    }

    if (length(rhs) == 1) rhs <- rep(rhs, nrow(linfct))
    if (length(rhs) != nrow(linfct))
        stop(sQuote("nrow(linfct)"), " is not equal to ",
             sQuote("length(rhs)"))

    RET <- list(model = model, linfct = linfct, rhs = rhs,
                coef = mpar$coef, vcov = mpar$vcov, 
                df = mpar$df, alternative = alternative,
                type = NULL)
    class(RET) <- "glht"
    RET
}

### symbolic description of H_0
glht.character <- function(model, linfct, ...) {
    ### extract coefficients and their covariance matrix
    beta <- try(coef(model))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("model"), " found!")

    tmp <- chrlinfct2matrix(linfct, names(beta))
    return(glht(model, linfct = tmp$K, rhs = tmp$m, 
                alternative = tmp$alternative))
}

### symbolic description of H_0
glht.expression <- function(model, linfct, ...) 
    glht(model, deparse(linfct), ...)

### multiple comparison procedures
glht.mcp <- function(model, linfct, ...) {

    ### extract factors and contrast matrices from `model'
    tmp <- mcp2matrix(model, linfct = linfct)

    args <- list(model = model, linfct = tmp$K)
    if (!is.null(tmp$alternative))
        args$alternative <- tmp$alternative
    if (any(tmp$m != 0))
        args$rhs <- tmp$m
    args <- c(args, list(...))

    ret <- do.call("glht", args)
    ret$type <- tmp$type
    ret$focus <- names(linfct)
    return(ret)
}
