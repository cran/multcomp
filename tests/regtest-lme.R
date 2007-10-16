
### mixed effects models
### feature request by John Wilkinson <jnwilks@btinternet.com>
### and Dieter Menne <dieter.menne@menne-biomed.de>

library("multcomp")
K <- rbind(c(0,1,-1,0),c(0,1,0,-1),c(0,0,1,-1))

nlmeOK <- require("nlme")
lme4OK <- require("lme4")
if (lme4OK) {

    data("ergoStool", package = "nlme")

    stool.lmer <- lmer(effort ~ Type + (1 | Subject),
                       data = ergoStool)
    glme4 <- glht(stool.lmer,K)

    if (nlmeOK) {
        stool.lme <- lme(effort ~ Type, data = ergoStool,
                        random = ~ 1 | Subject)
        gnlme <- glht(stool.lme,K)
        stopifnot(all.equal(coef(glme4), coef(gnlme)))
    }
}

### and now for lmer2 as well
if (lme4OK) {

    data("ergoStool", package = "nlme")

    stool.lmer <- lmer2(effort ~ Type + (1 | Subject),
                        data = ergoStool)
    glme4 <- glht(stool.lmer, K)

    if (nlmeOK) {
        stopifnot(all.equal(coef(glme4), coef(gnlme)))
    }
}

