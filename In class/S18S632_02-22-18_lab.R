# Review of Newton-Raphson method and 
# some results from summary of glm

# Let's find x value such that f(x)=0
# where f(x) = 2 - 1/exp(x)
# We'll start with x=0

f0 = function(x){
  return(4-exp(-x))
}
f1 = function(x){
  return(exp(-x))
}

x0 = 0
x1 = x0-f0(x0)/f1(x0)

x = seq(-3,0,length=1000)
plot(x,f0(x),type='l')
abline(0,0)
abline(f0(x0),f1(x0))
points(x1,0)

i.count = 1
print(c((i.count),x1))
while (abs(x1-x0) > 1e-6) {
  x0 = x1
  abline(f0(x1)-f1(x1)*x1,f1(x1),col=i.count)
  points(x1,0, col=i.count-1)
  Sys.sleep(5)
  x1 = x0-f0(x0)/f1(x0)
  i.count = i.count+1
  print(c((i.count),x1))
}


## Newton-Raphson method for many dimensions (betas)
## Let's use the example of coronary heart disease

library(faraway)
data(wcgs, package="faraway")


X = model.matrix(chd ~ height + cigs , wcgs)
y = as.numeric(wcgs$chd)-1

beta_t = rep(0,3)
phat = exp(X%*%beta_t)/(1+exp(X%*%beta_t))
beta_t1 = beta_t + solve(t(X)%*%diag(c(phat*(1-phat)))%*%X)%*%t(X)%*%(y-phat)

i.count = 1
print(c(i.count,beta_t1))
while (sum((beta_t1-beta_t)^2) > 1e-6){
  beta_t = beta_t1
  phat = exp(X%*%beta_t)/(1+exp(X%*%beta_t))
  beta_t1 = beta_t + solve(t(X)%*%diag(c(phat*(1-phat)))%*%X)%*%t(X)%*%(y-phat)
  i.count = i.count+1
  print(c(i.count,beta_t1))
}
m1 = glm(y~height+cigs, data = wcgs,family=binomial)
summary(m1)

varbeta = solve(t(X)%*%diag(c(phat*(1-phat)))%*%X) # X^t D X -1
varbeta
sqrt(diag(varbeta))
summary(m1)$coef


# Log-likelihood and Deviance 

ll.1 = sum(y*X%*%beta_t1-log(1+exp(X%*%beta_t1)))
Deviance.1 = -2*ll.1
Deviance.1

p0 = mean(y)
p0
ll.0 = sum(y*log(p0)+(1-y)*log(1-p0))
Deviance.0 = -2*ll.0
Deviance.0

sumary(m1)

# Saturated Model
y1 = y[y==1]
y0 = y[y==0]
L.S = prod(y1)*prod(1-y0)
ll.s = sum(log(y1))+sum(log(1-y0))
L.S
ll.s

####################################
####################################
# Poisson Regression
####################################


barplot(dpois(0:5,0.5),xlab="y",ylab="Probability",names=0:5,main="mean = 0.5")
barplot(dpois(0:10,2),xlab="y",ylab="Probability",names=0:10,main="mean = 2")
barplot(dpois(0:15,5),xlab="y",ylab="Probability",names=0:15,main="mean = 5")
data(gala, package="faraway")
summary(gala)

#Interested in all variables but Endemics (2nd variable)
gala = gala[,-2]


modl = lm(Species ~ . , gala)
plot(modl, 1)
modt = lm(sqrt(Species) ~ . , gala)
plot(modt, 1)

library(faraway)
sumary(modt)
modp = glm(Species ~ ., family=poisson, gala)
summary(modp)
sumary(modp)

halfnorm(residuals(modp))

# Checking to see if mean =  variance (Poisson characteristic)

plot(log(fitted(modp)),log((gala$Species-fitted(modp))^2), xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2))
abline(0,1)

# Dispersion parameter

(dp = sum(residuals(modp,type="pearson")^2)/modp$df.res)
sumary(modp,dispersion=dp)
modd = glm(Species ~ ., family=quasipoisson, gala)
sumary(modd)
drop1(modd,test="F")

## Purott and Reeder (1976)
## Experiment conducted to determine the effect of 
## gamma radiation on the numbers of chromosomal abnormalities
## ca: Number of chromosomal abnormalities. 
## cells: Number of cells in hundreds
## doseamt: The dose amount in Grays
## doserate: rate of dose in Grays/hour


data(dicentric, package="faraway")
round(xtabs(ca/cells ~ doseamt+doserate, dicentric),2)
with(dicentric,interaction.plot(doseamt, doserate, ca/cells, legend=T))
lmod = lm(ca/cells ~ log(doserate)*factor(doseamt), dicentric)

summary(lmod)
plot(residuals(lmod) ~ fitted(lmod),xlab="Fitted",ylab="Residuals")
abline(h=0)
dicentric$dosef = factor(dicentric$doseamt)

pmod = glm(ca ~ log(cells)+log(doserate)*dosef,  family=poisson,dicentric)
sumary(pmod)
rmod = glm(ca ~ offset(log(cells))+log(doserate)*dosef,  family=poisson,dicentric)
sumary(rmod)



data(solder, package="faraway")
modp = glm(skips ~ . , family=poisson, data=solder)
c(deviance(modp), df.residual(modp))
modp2  = glm(skips ~ (Opening +Solder + Mask + PadType + Panel)^2 ,  family=poisson, data=solder)
deviance(modp2)
pchisq(deviance(modp2),df.residual(modp2),lower=FALSE)
library(MASS)
modn = glm(skips ~ .,negative.binomial(1),solder)
modn
modn = glm.nb(skips ~ .,solder)
summary(modn)


library(pscl)
modp = glm(art ~ ., data=bioChemists, family=poisson)
sumary(modp)
ocount = table(bioChemists$art)[1:8]
pcount = colSums(predprob(modp)[,1:8])
plot(pcount,ocount,type="n",xlab="Predicted",ylab="Observed")
text(pcount,ocount, 0:7)
modh = hurdle(art ~ ., data=bioChemists)
summary(modh)
modz = zeroinfl(art ~ ., data=bioChemists)
summary(modz)
plot(fitted(modh), fitted(modz), xlab="Hurdle predictions", ylab="ZIP predictions")
abline(0,1)
modz2 = zeroinfl(art ~ fem+kid5+ment | ment, data=bioChemists)
summary(modz2)
(lrt = 2*(modz$loglik-modz2$loglik))
1-pchisq(6.1728,6)
exp(coef(modz2))
newman = data.frame(fem="Men",mar="Single",kid5=0,ment=6)
predict(modz2, newdata=newman, type="prob")
predict(modz2, newdata=newman, type="zero")


