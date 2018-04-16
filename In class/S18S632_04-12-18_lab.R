library(faraway)
data(pulp, package="faraway")
head(pulp, 12)
op <- options(contrasts=c("contr.sum", "contr.poly"))

lmod <- aov(bright ~ operator, pulp)
summary(lmod)
library(ggplot2)
ggplot(pulp, aes(x=operator, y=bright))+geom_point(position = position_jitter(width=0.1, height=0.0))
coef(lmod)
options(op)



library(lme4)
mmod <- lmer(bright ~ 1+(1|operator), pulp)
summary(mmod)
sumary(mmod)


## Predicting random effects and comparing with fixed effects


ranef(mmod)$operator
(cc <- model.tables(lmod))
cc[[1]]$operator/ranef(mmod)$operator




# BLUPs
fixef(mmod)+ranef(mmod)$operator


##############################
### Inference
################################

mmod <- lmer(bright ~ 1+(1|operator), pulp)
smod <- lmer(bright ~ 1+(1|operator), pulp, REML=FALSE)
summary(smod)
nullmod <- lm(bright ~ 1, pulp)

## Checking for random effects

lr.1 <- as.numeric(2*(logLik(smod)-logLik(nullmod)))
pvalue <- pchisq(lr.1,1,lower=FALSE)
data.frame(lr.1, pvalue)

## Parametric Bootstrap

y <- simulate(nullmod) #Simulate from the distribution under the null
lrstat <- numeric(1000)
set.seed(123)
for(i in 1:1000){
  y <- unlist(simulate(nullmod))
  bnull <- lm(y ~ 1)
  balt <- lmer(y ~ 1 + (1|operator), pulp, REML=FALSE)
  lrstat[i] <- as.numeric(2*(logLik(balt)-logLik(bnull)))
}

mean(lrstat < 0.00001)

phat = mean(lrstat > lr.1)
phat
sqrt(phat*(1-phat)/1000)


## Exact Test for Random Effects 

library(RLRsim)
exactLRT(smod, nullmod) #use exactLRT for models obtained with MLE
exactRLRT(mmod) #use exactRLRT for models obtained with REML

###################################################3
## Example 2
###################################################3

data(penicillin, package="faraway")
summary(penicillin)
penicillin$Blend <- gl(5,4)
ggplot(penicillin, aes(y=yield, x=treat, shape=Blend))+geom_point()+xlab("Treatment")
ggplot(penicillin, aes(y=yield, x=Blend, shape=treat))+geom_point()
op <- options(contrasts=c("contr.sum", "contr.poly"))
lmod <- aov(yield ~ blend + treat, penicillin)
summary(lmod)
coef(lmod)
mmod <- lmer(yield ~ treat + (1|blend), penicillin)
sumary(mmod)
options(op)

## Kenward-Roger Method (only to test fixed effects)

library(pbkrtest)
amod <- lmer(yield ~ treat + (1|blend), penicillin, REML=FALSE)
nmod <- lmer(yield ~ 1 + (1|blend), penicillin, REML=FALSE)
KRmodcomp(amod, nmod)
lrt.2 = as.numeric(2*(logLik(amod)-logLik(nmod)))
1-pchisq(lrt.2,3) #4 treatments, 4-1 = 3 df

## Parametric bootstrap

lrstat <- numeric(1000)
for(i in 1:1000){
  ryield <- unlist(simulate(nmod))
  nmodr <- refit(nmod, ryield)
  amodr <- refit(amod, ryield)
  lrstat[i] <- 2*(logLik(amodr)-logLik(nmodr))
}

mean(lrstat > lrt.2)

## Parametric bootstrap function (within R package pbkrtest)

pmod <- PBmodcomp(amod, nmod)
summary(pmod)

## Testing for random effect

rmod <- lmer(yield ~ treat + (1|blend), penicillin)
nlmod <- lm(yield ~ treat, penicillin)
lrt.3 = as.numeric(2*(logLik(rmod)-logLik(nlmod,REML=TRUE)))
lrstatf <- numeric(1000)
for(i in 1:1000){
  ryield <-  unlist(simulate(nlmod))
  nlmodr <- lm(ryield ~ treat, penicillin)
  rmodr <- refit(rmod, ryield)
  lrstatf[i] <- 2*(logLik(rmodr)-logLik(nlmodr,REML=TRUE))
}
mean(lrstatf < 0.00001)
mean(lrstatf > lrt.3)

library(RLRsim)
exactRLRT(rmod)
