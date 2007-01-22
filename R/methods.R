
### methods for `glht' objects
coef.glht <- function(object, rhs = FALSE, ...) 
{
    if (rhs) return(object$rhs)
    drop(object$linfct %*% object$coef)
}

vcov.glht <- function(object, ...) 
    object$linfct %*% tcrossprod(object$vcov, object$linfct)

summary.glht <- function(object, test = adjusted(), ...) {
    ts <- test(object)
    object$test <- ts
    class(object) <- switch(class(ts), "mtest" = "summary.glht",
                                       "gtest" = "summary.gtest")
    class(object) <- c(class(object), "glht")
    return(object)
}

confint.glht <- function(object, parm, level = 0.95, calpha = adjusted_calpha(), ...) 
{
    type <- attr(calpha, "type")
    if (is.function(calpha))
        calpha <- calpha(object, level, ...)
    if (!is.numeric(calpha) || length(calpha) != 1)
        stop(sQuote("calpha"), " is not a scalar")
    error <- attr(calpha, "error")
    attributes(calpha) <- NULL

    betahat <- coef(object)
    ses <- sqrt(diag(vcov(object)))
    switch(object$alternative, "two.sided" = {
            LowerCL <- betahat - calpha * ses
            UpperCL <- betahat + calpha * ses
        }, "less" = {
            LowerCL <- rep(-Inf, length(ses))
            UpperCL <- betahat - calpha * ses
        }, "greater" = {
            LowerCL <- betahat - calpha * ses
            UpperCL <- rep( Inf, length(ses))
    })

    ci <- cbind(LowerCL, UpperCL)
    colnames(ci) <- c("lower", "upper")
    object$confint <- cbind(betahat, ci)
    colnames(object$confint) <- c("Estimate", "lwr", "upr")
    attr(object$confint, "conf.level") <- level
    attr(object$confint, "calpha") <- calpha
    attr(object$confint, "error") <- error
    attr(object, "type") <- type
    class(object) <- c("confint.glht", "glht")
    return(object)
}
