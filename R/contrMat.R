# $Id: contrMat.R,v 1.19 2003/12/19 09:34:41 hothorn Exp $

contrMat <- function(n, type=c("Dunnett", "Tukey", "Sequen", "AVE",
                               "Changepoint", "Williams", "Marcus",
                               "McDermott","Tetrade"), nlevel=NULL, base = 1) {

    if (length(n) < 2) stop("less than 2 groups")
    type <- match.arg(type)
    if (any(n < 2)) stop("less than 2 observations in at least one group")
    k <- length(n)
    if (base < 1 || base > k) stop("base is not between 1 and ", k)
    CM <- c()
    rnames <- c()
    if (!is.null(names(n)))
        varnames <- names(n)
    else 
        varnames <- 1:length(n)

    kindx <- 1:k
#    if (base != 1 && type == "Dunnett") {
#      n <- c(n[base], n[-base])
#      varnames <- c(varnames[base], varnames[-base])
#      kindx <- c(base, (1:k)[-base])
#    }

    type <- match.arg(type)

    switch(type, "Dunnett" = {
        for(i in kindx[-base])
            CM <- rbind(CM, as.numeric(kindx == i) - as.numeric( kindx == base))
        rnames <- paste(varnames[kindx[-base]], "-", varnames[base], sep="")
    }, "Tukey" = {
        for (i in 1:(k-1)) {
            for(j in (i+1):k) {
                CM  <- rbind(CM, as.numeric(kindx==j)-as.numeric(kindx==i))
                rnames <- c(rnames, paste(varnames[j], "-", varnames[i],
                                          sep=""))
            }
        }
    }, "Sequen" =  {
        for (i in 2:k) {
            CM  <- rbind(CM, as.numeric(kindx==i)-as.numeric(kindx==i-1))
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
    rownames(CM) <- rnames
    if (type=="Tetrade")
      colnames(CM) <- NULL
    else 
      colnames(CM) <- varnames
    CM
}

contr.Dunnett <- function(n, base = 1, contrasts=TRUE) {
  if (length(n) == 1)  {
    if (contrasts) {
      mginv(contrMat(rep(10, n), base = base, type="Dunnett"))
    } else {
      diag(n)
    }
  }
  if (length(n) > 1) { 
    if (contrasts) {
      x <- rep(10, length(n))
      names(x) <- n
      mginv(contrMat(x, base = base, type="Dunnett")) 
    } else {
      diag(length(n))
    }
  } 
}

contr.Tukey <- function(n, contrasts=TRUE) {
  if (!contrasts) stop("contrasts is false")
  if (length(n) == 1) 
    mginv(contrMat(rep(10, n), type="Tukey"))
  if (length(n) > 1) {
    x <- rep(10, length(n))
    names(x) <- n
    mginv(contrMat(x, type="Tukey"))
  }
}
