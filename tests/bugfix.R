
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
### Rich Heiberger <rmh@temple.edu> 2006-11-28

### need to be identical
stopifnot(identical(cht1, print(cht1)))

### was: error
summary(cht1)$test


### NAs in coefficients
tmp.data <- data.frame(EE=gl(2, 1, 24, letters[1:2]),
                FF=gl(3, 2, 24, LETTERS[3:5]),
                GG=gl(4, 6, 24, letters[6:9]))
tmp.data$x <- rep(12, 24)
tmp.data$y <- rep(7, 24)
tmp.data$z <- c(9, 14, 3, 4, 15, 1, 11, 13, 24, 10, 22, 18,
                20, 21, 6, 7, 16, 2, 19, 12, 17, 8, 23, 5)
tmp.data$w <- c(15, 9, 18, 21, 17, 11, 23, 12, 1, 10, 2, 14, 24, 7,
                13, 4, 5, 19, 16, 20, 3, 8, 22, 6)

tmp.aov <- aov(z ~ EE+FF*GG + x*y +x*EE + y*FF, data=tmp.data)

try(glht(tmp.aov, linfct=mcp(EE="Tukey")))
try(glht(tmp.aov, linfct=mcp(FF="Tukey")))
glht(tmp.aov, linfct=mcp(GG="Tukey"))

### covariate interactions: fire a warning
tmp.aov <- aov(z ~ w*GG , data=tmp.data)
glht(tmp.aov, linfct = mcp(GG = "Tukey"))

### stop with informative error message
amod <- aov(breaks ~ tension + Error(wool), data = warpbreaks)
try(glht(amod, linfct = mcp(tension = "Tukey")))

### print error, spotted by Rich
amod <- aov(breaks ~ wool * tension, data = warpbreaks)
wht <- glht(amod, linfct = mcp(tension = "Tukey"))
tmp <- confint(wht, calpha=2)
print(tmp)

### coef. and vcov. didn't pass through
### bug report by John Deke <jdeke73@gmail.com>
lmod <- lm(Fertility ~ ., data = swiss) 
my.model <- list(coef(lmod),vcov(lmod)) 
coef2 <- function(model) return(model[[1]]) 
vcov2 <- function(model) return(model[[2]]) 
a <- glht(model = my.model, linfct = c("Agriculture=0","Catholic=0"),
          coef. = coef2, vcov. = vcov2, df = 100) 
b <- glht(model = lmod, linfct = c("Agriculture=0","Catholic=0"), 
          df = 100)
stopifnot(all.equal(coef(a), coef(b)))
