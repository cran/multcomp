
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

# tied pvalues in simtest
# spotted by Shin-ichi Hirata <001m9041@med.stud.kumamoto-u.ac.jp>
Response<-
c(0.0333333,0.0000000,0.1000000,0.0000000,0.0000000,0.1000000,0.4000000,
-0.3333333,0.1000000,0.0000000,0.2000000,0.0000000,-0.1666667,0.1000000,
-0.1333333,-0.0333333,0.0333333,0.0000000,0.0000000,0.0000000,0.0333333,
0.0000000,0.0000000,0.1000000,-0.1666667,0.0000000,-0.2333333,0.2000000,
0.0000000,0.2000000,0.0000000,-0.0666667,-0.1666667,0.1000000,0.2000000,
0.0000000,0.0000000,0.0000000,0.0000000,0.1000000,0.1000000,-0.1666667,
0.1000000,0.0000000,-0.2333333,0.0333333,0.4333333,0.1000000,0.2000000,
0.0000000,0.1000000)
Event<-factor(
c("b","b","b","b","b","b","c","c","c","c","c","c","c","d","d","d","d","d",
"d","e","e","e","e","e","e",
"f","f","f","f","f","f","f","g","g","g","g","g","g","h","h","h","h","h", 
"h","a","a","a","a","a","a","a"))
testdata<-cbind(as.data.frame(Response),Event)
simtest(Response~Event,data=testdata,type="Dunnett")

### test `subset' and `na.action' arguments
data(recovery)
simint(minutes ~ blanket, data=recovery, conf.level=0.9, 
       alternative="less", eps=0.01, subset = minutes > 7, na.action =
       na.fail)
simint(minutes ~ blanket, data=recovery, conf.level=0.9, 
       alternative="less", eps=0.01, subset = blanket != "b0", na.action =
       na.fail)

### spotted by Jamie Jarabek <jjarabek@stat.ufl.edu>
### fixed in mvtnorm_0.6-1 (univariate pmvt with df = 0)

x <- gl(3,10,30)
levels(x) <-c(" G1","G2","G3")
y <- rbinom(30,1,prob=c(rep(.8,10),rep(.2,10),rep(.5,10)))
toy.glm <- glm(y ~x, family=binomial)
csimint(estpar=coef(toy.glm)[2:3],df=toy.glm$df.residual,
        covm=vcov(toy.glm)[2:3,2:3],
        cmatrix=contrMat(c(10,10),type="Tukey"),asympt=TRUE)

