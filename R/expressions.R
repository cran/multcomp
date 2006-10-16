
### determine if an expression `x' can be interpreted as numeric
is_num <- function(x) {
    if (length(x) == 1) return(is.numeric(x))
    if (length(x) == 2) return(is.name(x[[1]]) && is.numeric(x[[2]]))
    return(FALSE)
}

### extract coefficients and variable names
coefs <- function(ex) {

    ### `a'
    if (length(ex) == 1 && !is_num(ex))
        return(list(coef = 1, var = as.character(ex)))

    ### `-a'
    if (length(ex) == 2 && (ex[[1]] == "-" && !is_num(ex[[2]])))
        return(list(coef = -1, var = as.character(ex[[2]])))

    if (length(ex) == 3) {

        ### `2 * a'
        if (ex[[1]] == "*" && (is_num(ex[[2]]) && !is_num(ex[[3]])))
            return(list(coef = eval(ex[[2]]), var = as.character(ex[[3]])))

        cf <- coefs(ex[[3]])
        if (ex[[1]] == "-") 
            cf$coef <- cf$coef * (-1)

        return(cf)
    }
    stop("cannot interpret expression ", sQuote(ex), " as linear function")
}

### extract left hand side of an expression
lhs <- function(ex) {

    if (length(ex) != 1)
        stop("expression is not of length 1")

    if (length(ex[[1]]) != 3)
        stop("expression ", sQuote(ex), 
             " does not contain a left and right hand side")

    return(ex[[1]][[2]])
}

### extract right hand side of an expression
rhs <- function(ex) {

    if (length(ex) != 1)
        stop("expression is not of length 1")

    if (length(ex[[1]][[3]]) == 2)
        return(-ex[[1]][[3]][[2]])

    rhs <- ex[[1]][[3]]
    if (!is_num(rhs) || length(rhs) > 1)
        stop("right hand side of expression ", sQuote(ex), 
             " is not a scalar numeric")
    return(rhs)
}

### extract direction of the _alternative_
side <- function(ex) {

    side <- as.character(ex[[1]][[1]])
    if (!(side %in% c("<=", ">=", "==", "=")))
        stop("does not contain ", sQuote("<=, >=, =="))
    alternative <- switch(side, 
        "<=" = "greater",
        ">=" = "less",
        "==" = "two.sided",
        "=" = "two.sided")
    return(alternative)
}

### extract coefficients and variable names
expression2coef <- function(ex) {

    cf <- c()
    nm <- c()

    m <- rhs(ex)
    x <- lhs(ex)
    alternative <- side(ex)

    while (TRUE) {

        tmp <- coefs(x)
        cf <- c(cf, tmp$coef)
        nm <- c(nm, tmp$var)

        ### x == "A"
        if (is.name(x)) break
        ### x == "-1"
        if (is_num(x)) break
        ### x == "3 * A"
        if (length(x) == 3 && is_num(x[[2]])) break
        ### x == "-3 * A"
        if (length(x) == 3 && is_num(x[[2]])) break
        x <- x[[2]]   
    }
    return(list(coef = cf, names = nm, 
                m = m, alternative = alternative, 
                lhs = deparse(lhs(ex))))
}

### interpret character representations of linear functions
chrlinfct2matrix <- function(ex, var) {

    if (!is.character(ex))
        stop("argument ", sQuote(ex), " is not of type character")
    if (!is.character(var))
        stop("argument ", sQuote(var), " is not of type character")

    K <- matrix(0, nrow = length(ex), ncol = length(var))
    colnames(K) <- var
    rownames(K) <- 1:length(ex)
    m <- rep(0, length(ex))

    for (i in 1:length(ex)) {

        tmp <- expression2coef(parse(text = ex[i]))

        if (!all(tmp$names %in% var))
            stop("variable(s) ", sQuote(tmp$names[!tmp$names %in% var]), 
                 " not found")

        for (n in tmp$names)
            K[i, var == n] <- tmp$coef[tmp$names == n]

        m[i] <- tmp$m

        if (i == 1) {
            alternative <- tmp$alternative
        } else {
            if (tmp$alternative != alternative)
                stop("mix of alternatives currently not implemented")
        }

        rownames(K)[i] <- tmp$lhs[1]
    }
    list(K = K, m = m, alternative = alternative)
}
