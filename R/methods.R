
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
    class(object) <- c(class(ts), class(object))
    return(object)
}

confint.glht <- function(object, parm, level = 0.95, ...) 
{
    pq <- pqglht(object)
    ci <- pq$qfunction(conf.level = level, ...)
    object$confint <- cbind(pq$coefficients, ci)
    colnames(object$confint) <- c("Estimate", "lwr", "upr")
    attr(object$confint, "conf.level") <- level
    attr(object$confint, "calpha") <- attr(ci, "calpha")
    attr(object$confint, "error") <- attr(ci, "error")
    class(object) <- c("confint.glht", "glht")
    return(object)
}
