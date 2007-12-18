
### mixed effects models
### feature request by John Wilkinson <jnwilks@btinternet.com>
### and Dieter Menne <dieter.menne@menne-biomed.de>

library("multcomp")

nlmeOK <- require("nlme")
lme4OK <- require("lme4")
if (lme4OK) {

    data("ergoStool", package = "nlme")
    K <- glht(aov(effort ~ Type, data = ergoStool), mcp(Type = "Tukey"))$linfct

    stool.lmer <- lmer(effort ~ Type + (1 | Subject),
                       data = ergoStool)
    glme4 <- glht(stool.lmer, K)
    glme41 <- glht(stool.lmer, mcp(Type = "Tukey"))
    stopifnot(all.equal(coef(glme4), coef(glme41)))

    if (nlmeOK) {
        stool.lme <- lme(effort ~ Type, data = ergoStool,
                        random = ~ 1 | Subject)
        gnlme <- glht(stool.lme,K)
        stopifnot(all.equal(coef(glme4), coef(gnlme)))

        gnlme2 <- glht(stool.lme, linfct = mcp(Type = "Tukey"))
        stopifnot(all.equal(coef(glme4), coef(gnlme2)))
    }
}

### and now for lmer2 as well
if (lme4OK) {

    data("ergoStool", package = "nlme")

    stool.lmer <- lmer2(effort ~ Type + (1 | Subject),
                        data = ergoStool)
    glme4 <- glht(stool.lmer, K)
    glme41 <- glht(stool.lmer, mcp(Type = "Tukey"))
    stopifnot(all.equal(coef(glme4), coef(glme41)))

    if (nlmeOK) {
        stopifnot(all.equal(coef(glme4), coef(gnlme)))
    }
}

