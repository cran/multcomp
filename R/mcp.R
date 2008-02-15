
### multiple comparison procedures for levels of 
### factors in AN(C)OVA models
mcp <- function(...) {

    linfct <- list(...)

    linfct <- lapply(linfct, function(x) {
        if (is.numeric(x) && !is.matrix(x)) {
            return(matrix(x, nrow = 1))
        } else {
            return(x)
        }})

    if (is.null(names(linfct)))
        stop(sQuote("linfct"), " doesn't have a ", sQuote("names"), 
             " attribute")

    classes <- sapply(linfct, function(x) inherits(x, "matrix") || 
                                          inherits(x, "character"))

    if (length(linfct) == 1 && linfct[[1]] == "Means") {
        class(linfct) <- "means"
        return(linfct)
    }

    if (all(classes)) {
        class(linfct) <- "mcp"
        return(linfct)
    }

    stop("Arguments don't consist of either matrices or characters")
}

### extract factors and contrast matrices used in `model'
factor_contrasts <- function(model) {

    ### extract model matrix, frame and terms
    mm <- try(model.matrix(model))
    if (inherits(mm, "try-error"))
        stop("no ", sQuote("model.matrix"), " method for ", 
             sQuote("model"), " found!")

    mf <- try(model.frame(model))
    if (inherits(mf, "try-error"))
        stop("no ", sQuote("model.frame"), " method for ", 
             sQuote("model"), " found!")

    tm <- try(terms(model))
    if (inherits(tm, "try-error"))
        stop("no ", sQuote("terms"), " method for ", 
             sQuote("model"), " found!")

    list(contrasts = attr(mm, "contrasts"),
         factors = attr(tm, "factors"),
         intercept = attr(tm, "intercept") != 0,
         mm = mm, 
         mf = mf)
}

### convert linear hypotheses supplied as single matrices,
### type arguments or expressions into one matrix
mcp2matrix <- function(model, linfct) {

    ### extract factors and contrasts
    fc <- factor_contrasts(model)
    contrasts <- fc$contrasts
    factors <- fc$factors
    intercept <- fc$intercept
    mf <- fc$mf
    mm <- fc$mm

    alternative <- NULL

    ### linear hypotheses
    if (!is.list(linfct) || is.null(names(linfct)))
        stop(sQuote("linfct"), "is not a named list")
    nhypo <- names(linfct)
    checknm <- nhypo %in% rownames(factors)
    if (!all(checknm)) 
        stop("Variable(s) ", sQuote(nhypo[!checknm]), " have been specified in ",
             sQuote("linfct"), " but cannot be found in ", sQuote("model"), "! ")
    if (any(checknm)) {
        checknm <- sapply(mf[nhypo[checknm]], is.factor)
        if (!all(checknm))
            stop("Variable(s) ", sQuote(paste(nhypo[!checknm], collapse = ", ")), " of class ", 
                  sQuote(paste(sapply(mf[nhypo[!checknm]], class), collapse = ", ")), 
                  " is/are not contained as a factor in ", sQuote("model"), ".")
    }
    m <- c()
    ctype <- c()
    for (nm in nhypo) {
        if (is.character(linfct[[nm]])) {

            Kchr <- function(kch) {
                ### check if kch is suitable as `type' argument to `contrMat'
                types <- eval(formals(contrMat)$type)
                pm <- pmatch(kch, types)
                ### if yes, compute K from `contrMat'
                if (!is.na(pm)) {
                    tmpK <- contrMat(table(mf[[nm]]), type = types[pm])
                    ctype <<- c(ctype, types[pm])
                } else {
                    ### if not, interpret kch as an expression
                    tmp <-  chrlinfct2matrix(kch, levels(mf[[nm]]))
                    tmpK <- tmp$K
                    m <<- c(m, tmp$m)
                    alternative <<- tmp$alternative
                }
                if (is.null(rownames(tmpK)))
                    rownames(tmpK) <- paste(kch, 1:nrow(tmpK), sep = "_")
                if (length(nhypo) > 1)
                    rownames(tmpK) <- paste(nm, rownames(tmpK), sep = ": ")
                list(K = tmpK)
            }
            
            tmp <- lapply(linfct[[nm]], Kchr)
            linfct[[nm]] <- do.call("rbind", lapply(tmp, function(x) x$K))
        }
    }

    ### transform linear hypotheses using model contrasts
    hypo <- vector(mode = "list", length = length(nhypo))
    names(hypo) <- nhypo

    for (nm in nhypo) {
        ### extract contrast matrix for each factor from model fit
        if (is.character(contrasts[[nm]])) {
            C <- do.call(contrasts[[nm]], 
                         list(n = nlevels(mf[[nm]])))
        } else {
            C <- contrasts[[nm]]
        }
        ### and transform the original linear hypotheses 
        ### K beta to K C beta^* 
        if (intercept) {
            Kstar <- linfct[[nm]] %*% C
        } else {
            ### model.matrix has `contrasts' argument even if no intercept
            ### was fitted and the contrast actually hasn't been applied
            Kstar <- linfct[[nm]]
        }
        pos <- factors[nm,] == 1
        ### average over interaction terms (if any)
        if (sum(pos) > 1) {
            Kinter <- c()
            for (i in which(pos)[-1]) {
                k <- sum(attr(mm, "assign") == i) / ncol(Kstar)
                ivar <- rownames(factors)[factors[ ,i] == 1]
                ivar <- ivar[ivar != nm]
                classes <- sapply(mf[, ivar, drop = FALSE], is.factor)
                if (all(classes)) {
                    fact <- 1 / (k + 1)
                } else {
                    fact <- 1
                    warning("covariate interactions found -- please choose appropriate contrast")
                }
                if (sum(factors[1:which(rownames(factors) == nm), i]) == 1) {
                    Kinter <- cbind(Kinter, 
                        Kstar[,rep(1:ncol(Kstar), k), drop = FALSE] * fact)
                } else {
                    Kinter <- cbind(Kinter, 
                        Kstar[,rep(1:ncol(Kstar), rep(k, ncol(Kstar))), 
                              drop = FALSE] * fact)
                }
            }
            Kstar <- cbind(Kstar, Kinter)
        }
        hypo[[nm]] <- list(K = Kstar,
            where = attr(mm, "assign") %in% which(factors[nm,] == 1))
    }

    ### combine all single matrices computed so far into
    ### one matrix of all linear hypoheses
    Ktotal <- matrix(0, nrow = sum(sapply(hypo, function(x) nrow(x$K))),
                     ncol = ncol(mm))
    colnames(Ktotal) <- colnames(mm)

    count <- 1
    for (h in hypo) {
        Ktotal[count:(count + nrow(h$K) - 1), h$where] <- h$K
        count <- count + nrow(h$K)
    }
    if (!is.matrix(Ktotal)) Ktotal <- matrix(Ktotal, nrow = 1)
    rownames(Ktotal) <- unlist(lapply(hypo, function(x) rownames(x$K)))

    if (is.null(ctype))
        ctype <- "User-defined"
    ctype <- paste(unique(ctype), collapse = ", ")
    attr(Ktotal, "type") <- ctype

    if (length(m) == 0) m <- 0
    list(K = Ktotal, m = m, alternative = alternative, type = ctype)
}

### contributed by Richard M. Heiberger <rmh@temple.edu>
meanslinfct <- function (model, focus, mmm.data = model$model, 
                         formula.in = terms(model))
{
    mmm.factor <- sapply(mmm.data, inherits, "factor")
    mmm.levels <- lapply(mmm.data[mmm.factor], levels)
    mmm.rows <- sapply(mmm.levels, length)
    n.mmm.rows <- prod(mmm.rows)
    mmm.new <- mmm.data[1:n.mmm.rows, ]
    mmm.factor.names <- names(mmm.data)[mmm.factor]
    mmm.rows.forward <- cumprod(mmm.rows)
    mmm.rows.forward.prev <- c(1, mmm.rows.forward)
    names(mmm.rows.forward.prev) <- c(names(mmm.rows.forward),
        "all")
    for (i in mmm.factor.names) mmm.new[[i]] <- gl(mmm.rows[i],
        mmm.rows.forward.prev[i], n.mmm.rows, labels = mmm.levels[[i]])
    mmm.numeric.names <- names(mmm.data)[!mmm.factor]
    for (i in mmm.numeric.names) mmm.new[[i]][] <- mean(mmm.data[[i]])
    none.data <- model.matrix(formula.in, data = mmm.new)
    none.linfct <- aggregate(none.data, by = mmm.new[focus],
        FUN = mean)[, -1]
    rownames(none.linfct) <- levels(mmm.new[[focus]])
    data.matrix(none.linfct)
}
