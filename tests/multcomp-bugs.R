
set.seed(290875)

library(multcomp)

data(detergent)

N <- rep(2, 5)

C <- contrMat(N, type="Tukey")

# should be equal
one <- simtest(plates ~ block+detergent, data=detergent, whichf="detergent",
        ttype="logical", cmatrix=C)
two <- simtest(plates ~ block+detergent, data=detergent, whichf="detergent",
        ttype="logical", type="Tukey")
three <- simtest(plates ~ block+detergent, data=detergent, ttype="logical",
        cmatrix=cbind(matrix(0, ncol=11, nrow=10), C))
stopifnot(round(max(abs(one$p.value.adj - two$p.value.adj)), 2) == 0)
stopifnot(round(max(abs(two$p.value.adj -  three$p.value.adj)), 2) == 0)


