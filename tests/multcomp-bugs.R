
actversion <- paste(R.version$major, R.version$minor, sep=".")
thisversion <- "1.7.0"
if (compareVersion(actversion, thisversion) >= 0)
  RNGversion("1.6.2")
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

# Contrasts for 2 levels only failed, spotted by 
# Peter B. Mandeville <mandevip@uaslp.mx>
load("TukeyTestData.rda")
Dunn <- round(simint(LAC~SC,data=TukeyTestData)$conf.int, 3)
Tukey <- round(simint(LAC~SC,data=TukeyTestData,type="Tukey")$conf.int,3)
HSD <- round(TukeyHSD(aov(LAC ~ SC, data=TukeyTestData))$SC[2:3],3)
attributes(Dunn) <- NULL
attributes(Tukey) <- NULL
attributes(HSD) <- NULL
stopifnot(all.equal(Dunn, Tukey))
stopifnot(all.equal(Dunn, HSD))
