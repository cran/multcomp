
library("multcomp")
set.seed(290875)

### mcp didn't accept objects of class `contrMat'
### spotted by Yves Brostaux <brostaux.y@fsagx.ac.be>
amod <- aov(response ~ trt, data = cholesterol)
cht1 <- glht(amod, linfct = mcp(trt = "Tukey"))
K <- contrMat(table(cholesterol$trt), type = "Tukey")
cht2 <- glht(amod, linfct = mcp(trt = K))
stopifnot(all.equal(coef(cht1), coef(cht2)))

