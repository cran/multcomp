
set.seed(290875)
library("multcomp")

### drei Gruppen mit je 1000 obs
n <- 1000
group <- gl(3, n)
levels(group) <- c("A", "B", "C")

### die Mittelwerte von drei Zielgroessen
means <- list(A = c(0, 0, 0), B = c(0, 1, 2), C = c(0, 2, 3))

### Kovarianzmatrix der Zielgroessen
S <- diag(3)
S[1,2] <- S[2,1] <- 0.5
S[3,2] <- S[2,3] <- 0.2


### Erzeugen der normalverteilten Zielgroessen
Y <- sapply(names(means), function(g) {
    rmvnorm(table(group)[g], mean = means[[g]], sigma = S)
}
)
colnames(Y) <- paste("y", 1:ncol(Y), sep = "")

### multivariates lineares Modell anpassen
mod <- lm(Y ~ group)

### Schaetzer + Kovarianzmatrix rausholen
beta <- as.vector(coef(mod))
sbeta <- vcov(mod)
names(beta) <- colnames(sbeta)

### Kontrastmatrix fuer Dunnett aufsetzen
K1 <- cbind(0, diag(2))
K <- rbind(cbind(K1, 0, 0, 0, 0, 0, 0),
           cbind(0, 0, 0, K1, 0, 0, 0),
           cbind(0, 0, 0, 0, 0, 0, K1))

colnames(K) <- names(beta)
rownames(K) <- names(beta)[c(2, 3, 5, 6, 8, 9)]

### ins multcomp schiessen
gmod <- glht(parm(beta, sbeta), lin = K)

### und freuen :-)
summary(gmod)
confint(gmod)


