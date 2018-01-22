### Lab1 S632

library(alr4)
library(leaps)

### Example 1: Roof Sales

data1 = read.csv(file="roof_sales.csv")
summary(data1)

y = data1$y

scatterplotMatrix(~ y + x1 + x2 + x3 + x4, diagonal="hist", smooth=F, reg.line=F, data1)
m1 =  lm(y ~ x1 + x2 + x3 + x4, data1)
summary(m1)
Anova(m1)

## Leaps and Bounds Algorithm

library(leaps) 
X=model.matrix(y~-1+x1+x2+x3+x4,data1)
y = data1$y
sets  =  leaps(X,y)

## Adjusted $R^2$

R2.vec  =  leaps(X,y,method="r2")$r2
adjR2.vec  =  leaps(X,y,method="adjr2")$adjr2

mod1  =  cbind(sets$which,R2.vec,rank(R2.vec),adjR2.vec,rank(adjR2.vec))
mod1

## PRESS


yhat.m1 = c(predict(m1))
e.m1 = c(resid(m1))
PRESS.1  =  (resid(m1))/(1-hatvalues(m1))

cbind(y, yhat.m1, e.m1, PRESS.1)

m23 = lm(y~x2+x3, data1)
m123 = lm(y~x1+x2+x3, data1)
PRESS.2  =  (y-predict(m23))/(1-hatvalues(m23))
PRESS.3  =  (y-predict(m123))/(1-hatvalues(m123))
cbind(resid(m23),PRESS.2,resid(m123),PRESS.3)

c(sum(PRESS.1^2),
  sum(PRESS.2^2),
  sum(PRESS.3^2),
  sum(abs(PRESS.1)),
  sum(abs(PRESS.2)),
  sum(abs(PRESS.3)))

PRESS  =  rep(0,15)
PRESS.abs  =  rep(0,15)
for(i in 1:15){
  X1  =  cbind(1,X[,sets$which[i,]])
  H  =  X1%*%solve(t(X1)%*%X1)%*%t(X1)
  y.hat  =  H%*%y
  hat.diag  =  diag(H)
  press.aux  =  (y-y.hat)/(1-hat.diag)
  PRESS[i]  =  sum(press.aux^2)
  PRESS.abs[i]  =  sum(abs(press.aux))
}
mod2  =  cbind(sets$which,adjR2.vec,PRESS,PRESS.abs)
mod2


## Cp

m23 = lm(y~x2+x3, data1)
dim(model.matrix(m23))
n = dim(model.matrix(m23))[1]
p = dim(model.matrix(m23))[2]
p.prime = dim(model.matrix(m1))[2]
s2.full = sigma(m1)^2
s2.red = sigma(m23)^2
Cp  =  p+(s2.red - s2.full)*(n-p)/s2.full
Cp

Cp.vec  =  leaps(X,y,method="Cp")$Cp
mod3  =  cbind(sets$which,adjR2.vec,PRESS,PRESS.abs,Cp.vec)
mod3



## AIC and BIC


sigma(m1)
n = dim(data1)[1]
p = length(coef(m1))
RSS.m1 = sum(resid(m1)^2)
aic = n*log(RSS.m1/n)+2*p
c(p,aic)
extractAIC(m1)

AIC(m1)

m0 = lm(y~1, data1)

extractAIC(m1)-extractAIC(m0)
AIC(m1)- AIC(m0)

bic = n*log(RSS.m1/n)+log(n)*p
c(p,bic)
extractAIC(m1,k = log(n))


aic.vec = rep(0,15)
bic.vec = rep(0,15)

for(i in 1:15){
  X1  =  cbind(X[,sets$which[i,]])
  mi = lm(y ~ X1)
  aic.vec[i] = extractAIC(mi)[2]
  bic.vec[i] = extractAIC(mi,k = log(n))[2]
}
mod4 = cbind(adjR2.vec,
             PRESS,
             Cp.vec,
             aic.vec,
             bic.vec)
round(mod4,3)

cbind(sets$which,apply(mod4,2,rank))

########################
# Stepwise methods
########################

## Forward Selection


m0 = lm(y ~ 1, data=data1)
m.fwd = step(m0, scope= ~ x1+x2+x3+x4, direction="forward")


## Backward Elimination

m1 = lm(y ~ x1+x2+x3+x4, data=data1)
m.bck = step(m1, scope = ~ 1, direction = "backward")


## Stepwise method

m.sw = step(m0, scope=list(lower=m0, upper=m1), direction = "both")


## Example 2: Mantel

### Forward Selection

m0 <- lm(Y ~ 1, data=mantel)
m.fwd = step(m0,scope= ~ X1 + X2 + X3, direction="forward")


### Backward Elimination

scatterplotMatrix(mantel)
m1 <- lm(Y ~ X1 + X2 + X3, data=mantel)
summary(m1)
m.bck = step(m1,scope=~1, direction="backward")


## Example 3: Highways

Highway$sigs1  =  with(Highway, (sigs * len + 1)/len)
f  =  ~ log(len) + shld + log(adt) + log(trks) + lane + slim +lwid +
  itg + log(sigs1) + acpt + htype

m0  =  lm(log(rate) ~ log(len), Highway) # the base model
m.fwd  =  step(m0, scope=f, direction="forward")
m1  =  update(m0, f)
m.bck  =  step(m1, scope = list(lower = ~ log(len), upper = m1), direction="backward")
m.sw.up  =  step(m0, scope=f)
list(lower=m0, upper=m1)
m.sw.down  =  step(m1, scope = list(lower = ~ log(len), upper = m1), direction = "both")
