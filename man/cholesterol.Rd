\name{cholesterol}
\alias{cholesterol}
\non_function{}
\title{ Cholesterol Reduction Data Set }
\usage{data(cholesterol)}
\description{
 Cholesterol reduction for five treatments; data set taken from 
 Westfall et al. (1999, p. 153). All pairwise comparisons according 
 to Tukey in a balanced one-way layout.
}
\format{
  This data frame contains the following variables
  \describe{
    \item{trt}{Treatment at 5 levels, see below.}
    \item{response}{Response variable.}
  }
}
\details{
  See Westfall et al. (1999, p. 153). The treatment formulations 
are defined as: 20mg one a day (\code{1time}), 10mg twice a day
(\code{2times})
and 5mg four times a day (\code{4times}). In addition, two competing drugs
were used as control (\code{drugD}) and (\code{drugE}).
}
\source{
  P. H. Westfall, R. D. Tobias, D. Rom, R. D. Wolfinger, Y. Hochberg (1999).
  \emph{Multiple Comparisons and Multiple Tests Using the SAS System}.
  Cary, NC: SAS Institute Inc.  
}
\examples{
data(cholesterol)

# adjusted p-values for all-pairwise comparisons in a one-way layout 
# tests for restricted combinations
simtest(response ~ trt, data=cholesterol, type="Tukey",
        ttype="logical")

# adjusted p-values all-pairwise comparisons in a one-way layout 
# (tests for free combinations -> p-values will be larger)
simtest(response ~ trt, data=cholesterol, type="Tukey",
        ttype="free")


# the following lines illustrate the basic principles of
# parameter estimation used in all functions in this package
# and how the low-level functions can be used with raw parameter
# estimates.

# the full design matrix (with reduced rank!)
x <- cbind(1, 
           matrix(c(rep(c(rep(1,10), rep(0,50)), 4), 
                    rep(1, 10)), nrow = 50))
y <- cholesterol$response

xpxi   <- mginv(t(x) \%*\% x)
rankx  <- sum(diag((xpxi \%*\% (t(x) \%*\% x))))
n      <- nrow(x)
p      <- ncol(x)
df     <- round(n-rankx)

# parameter estimates and their correlation
parm   <- xpxi \%*\% t(x) \%*\% y
mse    <- t(y-x \%*\% parm) \%*\% (y-x \%*\% parm)/df
covm   <- mse[1,1]*xpxi

# the contrast matrix
contrast <- contrMat(table(cholesterol$trt), type="Tukey")

# use the work-horse directly (and add zero column for the intercept)

csimint(estpar=parm, df=df, covm=covm, cmatrix=cbind(0, contrast))    
csimtest(estpar=parm, df=df, covm=covm, cmatrix=cbind(0, contrast),
         ttype="logical")      

# only a subset of all pairwise hypotheses:
#
# * drug D versus all other formulations and
# * all pairwise comparisions for "1time", "2times" and "4times"
#
csubset = contrast[c(1,3,5,6,8,10),]
csubset
simint(response ~ trt, data=cholesterol, cmatrix = csubset)

}
\keyword{datasets}
