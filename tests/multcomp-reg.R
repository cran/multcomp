 
library(multcomp)

set.seed(290875)

y <- rnorm(40)
x <- factor(c(rep(0,10), rep(1, 10), rep(2, 10), rep(3,10)))
lmod <- lm(y ~ x, contrasts=list(x = "contr.Dunnett"))
glmod <- glm(y ~ x, contrasts=list(x = "contr.Dunnett"), family=gaussian)

simint(y ~ x, type="Dunnett")
simint(y ~ x, type="Dunnett", base = 2)
simint(lmod, psubset=2:4)
simint(lmod, psubset=names(coef(lmod))[2:4])
cmatrix <- cbind(0, diag(3))
rownames(cmatrix) <- names(coef(lmod))[2:4]
simint(lmod, cmatrix=cmatrix)
stopifnot(all(round(confint(lmod)[2,], 3) == 
              round(simint(lmod, psubset=2)$conf.int, 3)))

simtest(y ~ x, type="Dunnett")
simtest(y ~ x, type="Dunnett", base = 2)
simtest(lmod, psubset=2:4)
simtest(lmod, psubset=names(coef(lmod))[2:4])
simtest(lmod, cmatrix=cmatrix)

simint(y ~ x, type="Dunnett")
simint(glmod, psubset=2:4, asympt=FALSE)
simint(glmod, psubset=names(coef(lmod))[2:4], asympt=FALSE)
simint(glmod, psubset=2:4)
simint(glmod, cmatrix=cmatrix)
simint(glmod, psubset=names(coef(lmod))[2:4])

simtest(y ~ x, type="Dunnett", asympt = TRUE)
simtest(glmod, psubset=2:4)
simtest(glmod, psubset=names(coef(lmod))[2:4])
