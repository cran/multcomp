
library("flexmix")
library("multcomp")


foo <- function(model) {

cf <- parameters(model)
vc <- refit(model)@vcov

cf <- cf[grep("^coef", rownames(cf)),]
vc <- vc[indx <- grep("_coef", rownames(vc)), indx]

tmp <- as.vector(outer(colnames(cf), rownames(cf), paste, sep = "_"))
cf <- as.vector(t(cf))
names(cf) <- tmp

colnames(vc) <- gsub("model.1_", "", colnames(vc))
rownames(vc) <- colnames(vc)

vc <- vc[names(cf), names(cf)]

k <- model@k
p <- length(cf) / k
tmp <- rep(3, k)
Ktmp <- contrMat(tmp, "Tukey")
K <- kronecker(diag(p), Ktmp)
colnames(K) <- names(cf)
rownames(K) <- 1:nrow(K)
for (i in 1:nrow(K))
    rownames(K)[i] <- paste(colnames(K)[K[i,] == 1], colnames(K)[K[i,] == -1], sep = "-")

glht(parm(coef = cf, vcov = vc), linfct = K)
}


data("NregFix")
Model <- FLXMRglm(~ x2 + x1)
fittedModel <- stepFlexmix(y ~ 1, model = Model, nrep=5, k = 3,
                           data = NregFix, concomitant = FLXPmultinom(~ w))
summary(refit(fittedModel))

ci <- foo(fittedModel)
plot(confint(ci))
