# $Id: contrMat.R,v 1.14 2002/04/10 06:46:39 hothorn Exp $

contrMat <- function(n, type=c("Dunnett", "Tukey", "Sequen", "AVE",
                               "Changepoint", "Williams", "Marcus",
                               "McDermott","Tetrade"), nlevel=NULL) {

    if (length(n) < 2) stop("less than 2 groups")
    if (any(n < 2)) stop("less than 2 observations in at least one group")
    k <- length(n)
    CM <- c()
    rnames <- c()
    if (!is.null(names(n)))
        varnames <- names(n)
    else 
        varnames <- 1:length(n)

    type <- match.arg(type)

    switch(type, "Dunnett" = {
        for(i in 2:k)
            CM <- rbind(CM, as.numeric((1:k)==i)-as.numeric((1:k)==1))
        rnames <- paste(varnames[2:k], "-", varnames[1], sep="")
    }, "Tukey" = {
        for (i in 1:(k-1)) {
            for(j in (i+1):k) {
                CM  <- rbind(CM, as.numeric((1:k)==j)-as.numeric((1:k)==i))
                rnames <- c(rnames, paste(varnames[j], "-", varnames[i],
                                          sep=""))
            }
        }
    }, "Sequen" =  {
        for (i in 2:k) {
            CM  <- rbind(CM, as.numeric((1:k)==i)-as.numeric((1:k)==i-1))
            rnames <- c(rnames, paste(varnames[i], "-", varnames[i-1],
                                      sep=""))
        }
    }, "AVE" = {
        help <- c(1,  -n[2:k]/sum(n[2:k]))
        CM <- rbind(CM, help)
        for (i in 2:(k-1)) {
            x <- sum(n[1:(i-1)])+sum(n[(i+1):k])
            help <- c(-n[1:(i-1)]/x, 1, -n[(i+1):k]/x)
            CM <- rbind(CM, help)
        }
        help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
        CM  <- rbind(CM, help)
        rnames <- paste("C", 1:nrow(CM))
    }, "Changepoint" = {
        for (i in 1:(k-1)) {
            help <- c(-n[1:i]/sum(n[1:i]), n[(i+1):k]/sum(n[(i+1):k]))
            CM <- rbind(CM, help)
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM), sep=""))
    }, "Williams" = {
        for (i in 1:(k-2)) {
            help <-  c(-1, rep(0, k-i-1), n[(k-i+1):k]/sum(n[(k-i+1):k]))
            CM <- rbind(CM, help)
        }
        help <- c(-1, n[2:k]/sum(n[2:k]))
        CM <- rbind(CM, help)
        rnames <- c(rnames, paste("C", 1:nrow(CM), sep=""))
    }, "Marcus" = {
        cm1 <- matrix(0, nrow=k-1, ncol=k)
        cm2 <- cm1
        for (i in 1:(k-1)) {
            cm1[i,(i+1):k] <- n[(i+1):k]/sum(n[(i+1):k])
            cm2[i,1:i] <- n[1:i]/sum(n[1:i])
        }
        row <- k*(k-1)/2
        index <- 1
        for (i in 1:(k-1)) {
            for (j in 1:i) {
                help <- cm1[i,]-cm2[j,]
                CM <- rbind(CM, help)
                index <- index+1
            }
        }
        rnames <- c(rnames, paste("C", 1:nrow(CM), sep=""))
#    }, "McDermott" = {
#        for(i in 1:(k-2)) {
#            help  <- c(-1, n[2:(i+1)]/sum(n[2:(i+1)]), rep(0, k-i-1))
#            CM <- rbind(CM, help)
#        }
#        help <- c(-1, n[2:k]/sum(n[2:k]) )
#        CM  <- rbind(CM, help)
#        rnames <- c(rnames, paste("C", 1:nrow(CM), sep=""))
     }, "McDermott" = {
         for(i in 1:(k-2)) {
             help  <- c(-n[1:i]/sum(n[1:i]), 1, rep(0, k-i-1))
             CM <- rbind(CM, help)
         }
         help <- c(-n[1:(k-1)]/sum(n[1:(k-1)]), 1)
         CM  <- rbind(CM, help)
         rnames <- c(rnames, paste("C", 1:nrow(CM), sep=""))
    }, "Tetrade" = {
        if (is.null(nlevel)) stop("nlevel missing")
        if (length(nlevel) != 2) stop("only two factors allowed")
        a <- nlevel[1]
        b <- nlevel[2]
	idi <- 1:a
	idj <- 1:b
        for (i1 in 1:(a-1)) {
            for (i2 in (i1+1):a) {
	        for (j1 in 1:(b-1)) {
        	    for (j2 in (j1+1):b) {
                	CM <- rbind(CM, kronecker( ( as.numeric(idi==i1)-as.numeric(idi==i2) ),
                                                   ( as.numeric(idj==j1)-as.numeric(idj==j2) ) ) ) 
		        rnames <- c(rnames, paste( "(", i1, j1, "-", i1, j2, ")", "-", 
                                                   "(", i2, j1, "-", i2, j2, ")",  sep=""))
            	    }
        	}
	    }
        }
    },)
    CM <- cbind(rep(0, nrow(CM)), CM)
    rownames(CM) <- rnames
    if (type != "Tetrade")
      colnames(CM) <- c("Intercept", varnames)
    else
      colnames(CM) <- NULL
    CM
}
