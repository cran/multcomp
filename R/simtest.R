# $Id: simtest.R,v 1.46 2003/02/18 08:27:56 hothorn Exp $

simtest <- function(y, ...) UseMethod("simtest")

simtest.default <- function(y, x=NULL, type=c("Dunnett", "Tukey",
                     "Sequen", "AVE", "Changepoint", "Williams", "Marcus",
                     "McDermott","Tetrade"), cmatrix=NULL,
                     alternative=c("two.sided","less", "greater"), asympt=FALSE,
                     ttype=c("free","logical"), eps=0.001, maxpts=1e+06,
                     nlevel=NULL, nzerocol=c(0,0), ...)

{
    ctype <- match.arg(type)
    alternative <- match.arg(alternative)

    xpxi   <- mginv(t(x) %*% x)
    rankx  <- sum(diag((xpxi %*% (t(x) %*% x))))
    n      <- nrow(x)
    p      <- ncol(x)
    df     <- ceiling(n-rankx)
    estpar <- xpxi %*% t(x) %*% y
    mse    <- t(y-x %*% estpar) %*% (y-x %*% estpar)/df
    covm   <- mse[1,1]*xpxi

    # compute the appropriate contrast matrix, if not given

    if (!is.null(cmatrix)) ctype <- "user-defined"
   
    if (is.null(cmatrix)) { 
      nobs <- apply(x[,(nzerocol[1] + 1):(ncol(x) - nzerocol[2])], 2, sum) 
      cm <- contrMat(nobs, ctype, nlevel)
      if (alternative == "greater")
        cm <- -cm
    } else {
      cm <- cmatrix
    }

    if (nzerocol[1] != 0)
      cm <- cbind(matrix(0, ncol=nzerocol[1], nrow=nrow(cm)), cm)
    if (nzerocol[2] != 0)
      cm <- cbind(cm, matrix(0, ncol=nzerocol[2], nrow=nrow(cm)))

    csimtest(estpar, df, covm, cm, ctype, ttype,
         alternative,  asympt, eps, maxpts)
}


csimtest <- function(estpar, df, covm, cmatrix=NULL, ctype="user-defined",
                    ttype=c("free","logical"), alternative=c("two.sided",
                    "less","greater"), asympt=FALSE,
                    eps=0.001, maxpts=1000000)
{
    if (!is.vector(estpar) & !is.matrix(estpar)) stop("estpar not a vector")
    p <- length(estpar)
    if (!all.equal(df - floor(df), 0)) stop("df not an integer")
    if (!is.matrix(covm)) stop("covm is not a matrix")
    cm <- cmatrix
    if (ctype !="user-defined") cmatrix <- NULL

    alternative <- match.arg(alternative)
    ttype <- match.arg(ttype)

    if (asympt) df <- 0                          
 
    covm  <- cm %*% covm %*% t(cm)                            
    d     <- diag(1/sqrt(diag(covm)))              
    cr    <- d %*% covm %*% d

    ests  <- cm %*% estpar
    ses   <- sqrt(diag(covm))
    tvals <- ests/ses
    dim   <- ncol(cr)
    delta <- rep(0,dim)

    switch(alternative, "two.sided" = {
	tvals <- -abs(tvals)
        if (df>0) rawp <- 2*pt(tvals,df)     
        else      rawp <- 2*pnorm(tvals)
    }, "less" = {
        if (df>0) rawp <- pt(tvals,df) 
        else      rawp <- pnorm(tvals)
    }, "greater" = {
       if (df>0) rawp <- pt(tvals,df)
       else      rawp <- pnorm(tvals)
    },)


    pvals    <- as.matrix(rawp)
    tvals    <- as.matrix(tvals)    
    covcont  <- covm

    k        <- nrow(cm)
    r        <- t(rank(pvals))
    ir       <- r
    ir[,r]   <- 1:nrow(pvals)
    origord  <- t(ir)         
    cord     <- cm[ir,]
    tvalsord <- tvals[ir,]
    pvalsord <- pvals[ir,]
    ccord    <- covcont[ir,ir]
    crrccord <- solve(sqrt(diag(diag(ccord)))) %*% ccord %*% solve(sqrt(diag(diag(ccord))))
    cct      <- t(cord)
   
    if (ttype == "logical" & k > 2) {
        for (iout in 1:(k-2)) {
	    limit <- 2^(k-iout-1)
	    in2   <- rep(0, k-iout-1)
	    zero  <- as.matrix(rep(0, k))
	    in1   <- as.matrix(zero)
	    y     <- cct[,1:iout]

            for (kk in 1:limit) {
                if (kk == limit) {
                    in2 <- rep(0, k-iout-1)
		}
		else { 
		    ii <- 1
	            zz <- kk
                    while( zz %% 2 == 0) {    
		        ii<-ii+1
                        zz<-zz/2
                    }
                    if (in2[ii] == 0) 
                        in2[ii] <- 1
                    else 
                        in2[ii] <- 0 
                }

                locbin  <- as.matrix( rbind( as.matrix(rep(0, iout)), 
                                             as.matrix(rbind(1, as.matrix(in2) ) ) ) )
                loc1 <- 1
                for (jj in 1:nrow(locbin)) {
                    if (locbin[jj] == 0) {
                    }
                    else {
                        loc1 <- c(loc1, jj)
                    }
                }
                loc1 <- t(loc1)
                loc1 <- t(loc1[2:ncol(loc1)])
                x    <- as.matrix(cct[,loc1])
                res    <- y - x %*% mginv(t(x) %*% x) %*% t(x) %*% y
                ssemat <- diag(t(res) %*% res)
                if (all(ssemat > .00000001)) {
                    if (identical(in1,zero)) {
                        in1 <- locbin
                        }
                    else {
                        check <- in1 - rep(locbin, ncol(in1))
                        diff  <- apply(check,2,max) - apply(check,2,min)
                        if (min(diff) == 2)
                            in1 <- cbind(as.matrix(in1),as.matrix(locbin))
                        else {
                            mindx <- which.min(diff)
                            if (sum(check[,mindx]) == -1) 
                                in1[,mindx] <- locbin
                        }
                    }
                }
            }
            in1   <- t(in1)
            ncont <- nrow(in1)
            in1   <- cbind(as.matrix(rep(iout+1,ncont)),in1)
            if (iout == 1)
                inbig <- in1
            else
                inbig <- rbind(inbig, in1)
        }
        big <- rbind(t(rep(1,k+1)),inbig)   
    }
    

    if (ttype == "free" | k < 3)
         big <- t(rep(1,k+1))
    lastset <-  cbind( cbind( rep(k,1), t(rep(0,k-1)) ), rep(1,1) ) 
    big <- rbind(big, lastset)
    stepj <- big[,1]
    if (ttype == "free") 
        stepj <- 1:k
    SubsetK <- big[,2:ncol(big)]
    if (ttype == "free") {
        m <- t(rep(1,k))
        for (i in 2:k) {
            r <- cbind( t(rep(0,i-1)), t(rep(1,k-i+1)) )
            m <- rbind(m,r)
        }
        SubsetK <- m
    }
    nbig <- nrow(big)
    if (ttype == "logical") {
        yyy <- rnorm(nrow(big))
        aaa <- as.factor(big[,1])
    }
    else {
        yyy <- rnorm(nrow(as.matrix(stepj)))
        aaa <- as.factor(stepj)
    }   
    contrasts(aaa) <- "ct"
    ff <- (yyy ~ aaa)
    mf <- model.frame(ff)
    des <- model.matrix(ff, mf)
    des <- des[,2:ncol(des)]
    if (ttype == "logical")
        contonly <- big[,2:(k+1)]
    else 
        contonly <- SubsetK
    tcmpr <- des %*% tvalsord


    if (ttype == "free") 
        nbig <- k
    count   <- rep(0,nbig)
    countc  <- count
    countc2 <- count
    nnn <- ncol(crrccord)
    for (i in 1:nnn) {
        for (j in 1:nnn) {
            if ( (contonly[i,j] + des[i,j])  > 0 )
                contonly[i,j] <- 1
            else
                contonly[i,j] <- 0
        }
    }

    stepj   <- as.matrix(stepj)
    subsets <- cbind(contonly, stepj)
    gls     <- rep(0,nrow(contonly))
    stdgls  <- rep(0,nrow(contonly))

    for (i1 in 1:nrow(contonly)) {
        loct <- 1
        for (jj in 1:ncol(contonly)) {
            if (contonly[i1,jj] != 0)
                loct <- c(loct, jj)
        }
        loct <- t(loct)
        loct <- t(loct[2:ncol(loct)]) 
        cort  <- as.matrix(crrccord[loct,loct])
        nt    <- nrow(cort)
        low3  <-  tcmpr[i1]* rep(1,nt)
        upp3  <- -tcmpr[i1]* rep(1,nt)
        if (alternative != "two.sided")
           upp3 <- rep(Inf,nt)
#        maxpts <- 1000000
        delta  <- rep(0,nt) 
        prob       <- pmvt(lower=low3, upper=upp3, df=df, 
                           delta=delta, corr=cort, abseps=eps, maxpts=maxpts)
        gls[i1]    <- 1-prob
        stdgls[i1] <- attr(prob, "error")
    }

    glsbig <- des * ( as.matrix(gls) %*% t(rep(1,k)) )
    glsp   <- apply(glsbig,2,max)
    glsin  <- 1:ncol(glsbig)
    for (i in 1:ncol(glsbig)) {
        glsin[i]  <- which.max(glsbig[,i])
    }
    glsin  <- t(glsin)  
    stdgls <- t(stdgls)  
    stdgls  <- stdgls[glsin]

    for (i in 2:k) {
        if (glsp[i] < glsp[i-1]) {
            glsp[i]   <- glsp[i-1]
            stdgls[i] <- stdgls[i-1]
        }
    }

    seadjp  <- as.matrix(stdgls)

    # <FIXME>
    # is a warning appropriate?
    if (any(seadjp > eps)) warning("error > eps")
    # </FIXME>

    adjpgls <- as.matrix(glsp)
    adjp    <- as.matrix(adjpgls)

    if (ttype == "logical")
        totals  <-  apply(contonly,1,sum)
    else 
        totals <- k:1
    if (alternative == "two.sided") {
        if (df  == 0)
            bon <- 2*(pnorm(tcmpr)) * totals
        else 
            bon <- 2*(pt(tcmpr,df)) * totals
    }       
    else {
        if (df == 0) 
            bon <- (pnorm(tcmpr)) * totals
        else 
            bon <- (pt(tcmpr,df)) * totals
    }      

    bonbig  <- des * ( as.matrix(bon) %*% t(rep(1,k)) )
    bonp    <- apply(bonbig,2,max)
    bonmult <- bonp/pvalsord

    for (i in 2:k) {
        if (bonp[i] < bonp[i-1]) 
            bonp[i] <- bonp[i-1]
    }
    rawp     <- pvalsord
    ests <- ests[ir,]
    tvals<- tvals[ir,]
    if (alternative == "greater") {
        ests <- -ests
        cm <- -cm
        tvals <- -tvals
    }
    adjpbon  <- bonp
    adjpbon  <- pmin(1,adjpbon)
    rownames(adjp) <- rownames(cord)

    RET <- list(cmatrix = cm, ctype = ifelse(is.null(cmatrix), ctype, NA),
                estimate = ests, sd = ses, statistics = tvals,
                p.value.raw = rawp, p.value.bon = bonp,
                p.value.adj = adjp, eps=eps)

    class(RET) <- "hmtestp"
    RET
}

simtest.formula <-
function(formula, data=list(), subset, na.action, whichf, ...)
{
    cl <- match.call()
    pf <- parseformula(formula, data, subset, na.action, whichf, ...)
    x <- pf$x
    y <- pf$y
    cargs <- pf$cargs

    attr(x, "contrasts") <- NULL
    attr(x, "assign") <- NULL
    y <- do.call("simtest", c(list(y=y, x=x), cargs))
    y$DNAME <- cl
    y$FNAMES <- pf$fnames
    y
}


