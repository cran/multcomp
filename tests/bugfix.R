
library("multcomp")
set.seed(290875)

### mcp didn't accept objects of class `contrMat'
### spotted by Yves Brostaux <brostaux.y@fsagx.ac.be>
amod <- aov(response ~ trt, data = cholesterol)
cht1 <- glht(amod, linfct = mcp(trt = "Tukey"))
K <- contrMat(table(cholesterol$trt), type = "Tukey")
cht2 <- glht(amod, linfct = mcp(trt = K))
stopifnot(all.equal(coef(cht1), coef(cht2)))

### several inconsistencies spotted by 
### Rich Heiberger <<rmh@temple.edu> 2006-11-28

### need to be identical
stopifnot(identical(cht1, print(cht1)))

### was: error
summary(cht1)$test

