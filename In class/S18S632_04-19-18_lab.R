options(digits = 4)

library(faraway)

# Example 5: PSID


data(psid, package="faraway")
head(psid)
library(dplyr)
psid20 = filter(psid, person <= 20)
library(ggplot2)
ggplot(psid20, aes(x=year, y=income))+geom_line()+facet_wrap(~ person)
ggplot(psid20, aes(x=year, y=income+100, group=person)) +geom_line() + facet_wrap(~ sex) + scale_y_log10()
lmod = lm(log(income) ~ I(year-78), subset=(person==1), psid)
coef(lmod)
library(lme4)
ml = lmList(log(income) ~ I(year-78) | person, psid)
intercepts = sapply(ml,coef)[1,]
slopes = sapply(ml,coef)[2,]
plot(intercepts,slopes,xlab="Intercept",ylab="Slope")
psex = psid$sex[match(1:85,psid$person)]
boxplot(split(slopes,psex))
t.test(slopes[psex=="M"],slopes[psex=="F"])
t.test(intercepts[psex=="M"],intercepts[psex=="F"])
library(lme4)
psid$cyear = psid$year-78
mmod = lmer(log(income) ~ cyear*sex +age+educ+(cyear|person),psid)
sumary(mmod, digits=3)
library(pbkrtest)
mmod = lmer(log(income) ~ cyear*sex +age+educ+(cyear|person),psid, REML=FALSE)
mmodr = lmer(log(income) ~ cyear + sex +age+educ+(cyear|person),psid, REML=FALSE)
KRmodcomp(mmod,mmodr)
confint(mmod, method="boot")
diagd = fortify(mmod)
ggplot(diagd,aes(sample=.resid))+stat_qq()+facet_grid(~sex)
diagd$edulevel = cut(psid$educ,c(0,8.5,12.5,20), labels=c("lessHS","HS","moreHS"))
ggplot(diagd, aes(x=.fitted,y=.resid)) + geom_point(alpha=0.3) + geom_hline(yintercept=0) + facet_grid(~ edulevel) + xlab("Fitted") + ylab("Residuals")


# Example 6: Vision


data(vision, package="faraway")
vision$npower = rep(1:4,14)
ggplot(vision, aes(y=acuity, x=npower, linetype=eye)) + geom_line() + facet_wrap(~ subject, ncol=4) + scale_x_continuous("Power",breaks=1:4,labels=c("6/6","6/18","6/36","6/60"))
mmod = lmer(acuity~power + (1|subject) + (1|subject:eye),vision)
sumary(mmod)
4.64^2/(4.64^2+3.21^2+4.07^2)
(4.64^2+3.21^2)/(4.64^2+3.21^2+4.07^2)
library(pbkrtest)
mmod = lmer(acuity~power+(1|subject)+(1|subject:eye),vision,REML=FALSE)
nmod = lmer(acuity~1+(1|subject)+(1|subject:eye),vision,REML=FALSE)
KRmodcomp(mmod, nmod)
mmodr = lmer(acuity~power+(1|subject)+(1|subject:eye),vision,REML=FALSE, subset=-43)
nmodr = lmer(acuity~1+(1|subject)+(1|subject:eye),vision,REML=FALSE, subset=-43)
KRmodcomp(mmodr, nmodr)
op = options(contrasts=c("contr.helmert", "contr.poly"))
mmodr = lmer(acuity~power+(1|subject)+(1|subject:eye),vision,subset=-43)
sumary(mmodr)
options(op)
contr.helmert(4)
plot(resid(mmodr) ~ fitted(mmodr),xlab="Fitted",ylab="Residuals")
abline(h=0)
qqnorm(ranef(mmodr)$"subject:eye"[[1]],main="")


# Example 7: jsp

# Part 1

data(jsp, package="faraway")
jspr = jsp[jsp$year==2,] # respons is third year score
ggplot(jspr, aes(x=raven, y=math))+xlab("Raven Score")+ylab("Math Score")+geom_point(position = position_jitter(),alpha=0.3)
ggplot(jspr, aes(x=social, y=math))+xlab("Social Class")+ylab("Math Score")+geom_boxplot()
glin = lm(math ~ raven*gender*social,jspr)
library(car)
Anova(glin)
glin = lm(math ~ raven*social,jspr)

Anova(glin) # let's use sig.level alpha = 0.01

glin = lm(math ~ raven+social,jspr)
summary(glin)
table(jspr$school)
mmod = lmer(math ~ raven*social*gender+(1|school)+(1|school:class),  data=jspr)
mmodr = lmer(math ~ raven*social+(1|school)+(1|school:class),  data=jspr)
KRmodcomp(mmod, mmodr)
all3 = lmer(math ~ raven*social*gender+(1|school)+(1|school:class), data=jspr, REML=FALSE)
all2 = update(all3, . ~ . - raven:social:gender)
notrs = update(all2, . ~ . -raven:social)
notrg = update(all2, . ~ . -raven:gender)
notsg = update(all2, . ~ . -social:gender)
onlyrs = update(all2, . ~ . -social:gender - raven:gender)
all1 =  update(all2, . ~ . -social:gender - raven:gender - social:raven)
nogen = update(all1, . ~ . -gender)
anova(all3, all2, notrs, notrg, notsg, onlyrs, all1, nogen)[,1:4]
jspr$craven = jspr$raven-mean(jspr$raven)
mmod = lmer(math ~ craven*social+(1|school)+(1|school:class),jspr)
sumary(mmod)
diagd = fortify(mmod)
ggplot(diagd,aes(sample=.resid))+stat_qq()
ggplot(diagd,aes(x=.fitted,y=.resid)) +geom_point(alpha=0.3) +geom_hline(yintercept=0) +xlab("Fitted") +ylab("Residuals")
qqnorm(ranef(mmod)$school[[1]],main="School effects")
qqnorm(ranef(mmod)$"school:class"[[1]],main="Class effects")
adjscores = ranef(mmod)$school[[1]]
rawscores = coef(lm(math ~ school-1,jspr))
rawscores = rawscores-mean(rawscores)
plot(rawscores,adjscores)
sint = c(9,14,29)
text(rawscores[sint],adjscores[sint]+0.2,c("9","15","30"))


library(RLRsim)
mmodc = lmer(math ~ craven*social+(1|school:class),jspr)
mmods = lmer(math ~ craven*social+(1|school),jspr)
exactRLRT(mmodc, mmod, mmods)
exactRLRT(mmods, mmod, mmodc)
schraven = lm(raven ~ school, jspr)$fit
mmodc = lmer(math ~ craven*social+schraven*social+(1|school)+  (1|school:class),jspr)
KRmodcomp(mmod, mmodc)

# Part 2

data(jsp, package="faraway")
jspr = jsp[jsp$year==2,]
mjspr = data.frame(rbind(jspr[,1:6],jspr[,1:6]),  subject=factor(rep(c("english","math"),c(953,953))),  score=c(jspr$english/100,jspr$math/40))
ggplot(mjspr, aes(x=raven, y=score))+geom_jitter(alpha=0.25)+facet_grid(gender ~ subject)
mjspr$craven = mjspr$raven-mean(mjspr$raven)
mmod = lmer(score ~ subject*gender + craven*subject + social + (1|school) + (1|school:class) + (1|school:class:id),mjspr)
sumary(mmod)
library(pbkrtest)
mmod = lmer(score ~ subject*gender+craven*subject+social+  (1|school)+(1|school:class)+(1|school:class:id),mjspr, REML=FALSE)
mmodr = lmer(score ~ subject*gender+craven+subject+social+(1|school)+(1|school:class)+(1|school:class:id),mjspr, REML=FALSE)
KRmodcomp(mmod, mmodr)
0.10^2/(0.10^2+0.12^2)
diagd = fortify(mmod)
ggplot(diagd, aes(x=.fitted,y=.resid)) + geom_point(alpha=0.3) + geom_hline(yintercept=0) + facet_grid(~ subject) + xlab("Fitted") + ylab("Residuals")
