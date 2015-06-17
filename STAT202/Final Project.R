# 1
# a
# i
rr=40/46/36*61
rr
# ii
or=(40*25)/(6*36)
or
# iii
se=sqrt(1/40+1/6+1/36+1/25)
se
log(or)+c(-1,1)*qnorm(.975)*se
exp(log(or)+c(-1,1)*qnorm(.975)*se)
# iv
data=matrix(c(40,36,6,25),2,2)
test=chisq.test(data,correct=F)
G.square = 2*sum(data*log(data/test$expected))
G.square
1-pchisq(G.square,1)

# b
Group=c(rep(1,46),rep(0,61))
Tumor=c(rep(1,40),rep(0,6),rep(1,36),rep(0,25))
mod=glm(Tumor~Group, family=binomial)
summary(mod)$coef
c=mod$coef
r1=exp(c[1]+c[2])/(1+exp(c[1]+c[2]))
r2=exp(c[1])/(1+exp(c[1]))
rr=r1/r2
rr
or=exp(c[2])
source('http://www.ics.uci.edu/~staceyah/111-202/R/Stats111-202-Rfunctions.R')
glmCI(mod)


# 2
frogs=read.csv('http://www.ics.uci.edu/~staceyah/111-202/data/frogs.csv')
attach(frogs)
mod1=glm(Mates~BodySize+Species, family=poisson, data=frogs)
mod1$coef
exp(mod1$coef)$[2]

mod2=glm(Mates~BodySize, family = poisson, data=frogs)
anova(mod2, mod1, test='LRT')

mod3=glm(Mates~(.)^2, family = poisson, data = frogs)
mod3$coef

anova(mod1,mod3,test='LRT')


plot(Mates~BodySize, type='n',main='Number of Mates vs. Body Size',xlab='Body Size (mm)',ylab='Number of Mates')
points(Mates[Species=='A']~BodySize[Species=='A'],pch=16,col='red')
points(Mates[Species=='B']~BodySize[Species=='B'],pch=17,col='blue')
c=mod3$coef
curve(exp(c[1]+c[2]*x),col='red',lty=2,add=T)
curve(exp(c[1]+c[3]+(c[2]+c[4])*x),col='blue',lty=3,add=T)
legend(locator(1),c('Species A','Species B'),pch=c(16,17),col=c('red','blue'),lty=c(2,3))


# 3
R code for question 3:
  outdoor=read.csv('http://www.ics.uci.edu/~staceyah/111-202/data/asthma.csv')
outdoor$actgrp=factor(outdoor$actgrp)
attach(outdoor)

########################## part (a) ##########################
######################## Data Summary #######################

# Summary of response:
table(lofev)

# Range/summary of quant. predictors:
summary(age)
sd(age)
summary(height)
sd(height)

# Stratified by normal/low bwt:
tapply(age,lofev,mean)
tapply(age,lofev,sd)

tapply(height,lofev,mean)
tapply(height,lofev,sd)


# Summary of cat. predictors:
table(male)
table(asthma)
table(actgrp)

# Mainly interested in how predictors are related to actgrp.
# Side-by-side boxplots with quantitative predictors:
par(mfrow=c(1,2))
boxplot(age~factor(lofev,levels=c(0,1),labels=c("Normal","Low")),ylab="Age of Child",
        main="Age vs. Normal or Low FEV1")
boxplot(height~factor(lofev,levels=c(0,1),labels=c("Normal","Low")),
        ylab="Height of Child",main="Height vs. Normal or Low FEV1")

table(male,lofev)
table(asthma,lofev)
table(actgrp,lofev)

table(male,lofev)/rowSums(table(male,lofev))
table(asthma,lofev)/rowSums(table(asthma,lofev))
table(actgrp,lofev)/rowSums(table(actgrp,lofev))

########################## part (b) ##########################
######################## Mpdel Building #######################

## Look at each predictor by itself:
## First look at the predictor of inerest: actgrp
mod.actgrp=glm(lofev~actgrp, family=binomial)
summary(mod.actgrp) # Coef actgrp2=-0.2377; p-value=0.0992
# Coef actgrp3=-0.6252; p-value=0.0001
anova(mod.actgrp,test='LRT') # 0.0005084 ***

# other predictors
mod.age=glm(lofev~age, family=binomial)
summary(mod.age) # Coef=-0.9765; p-valu<2e-16
anova(mod.age,test='LRT') # p-value< 2.2e-16 ***

mod.height=glm(lofev~height, family=binomial)
summary(mod.height)

mod.male=glm(lofev~male,family=binomial)
summary(mod.male) # Coef=-0.2392; p-value<2e-16
anova(mod.male,test='LRT') # p-value=0.02822 *

mod.asthma=glm(lofev~asthma,family=binomial)
summary(mod.asthma) # Coef=0.1678; p-value=0.243
anova(mod.asthma,test='LRT') # p-value=0.2466

# Among all predictors, only asthma is not maginally significant.

# Consider the model with all main effects
mod.main=glm(lofev~.,family=binomial,data=outdoor)
summary(mod.main)
anova(mod.main,test='LRT') #actgrp not significant
drop1(mod.main,test='LRT') #age and actgrp not significant, given other variables are in the model

# Model with all two-way interactions
mod.int=glm(lofev~(.)^2,family=binomial, data=outdoor)
summary(mod.int)

# add single predictors to the model with only actgrp
add1(mod.actgrp,scope=mod.main,test='LRT') #age, height, male are significant
mod1=update(mod.actgrp,.~.+age+height+male) # Add age, height and male?
anova(mod.actgrp,mod1,test='LRT') # p-value<2.2e-16; Yes
drop1(mod1,test='LRT') # actgrp and age are not significant
mod2=update(mod1,.~.-age) #drop age?
anova(mod1,mod2,test='LRT') # p-value=0.7297; Yes

add1(mod2,scope=mod.main,test='LRT') # asthma is significant
mod3=update(mod2,.~.+asthma) # add asthma?
anova(mod2,mod3, test='LRT') # p-value=0.02167; Yes

anova(mod3,mod.main,test='LRT') # p-value for adding age=0.7769,not add age

#add interactions?
add1(mod3,scope=mod.int,test='LRT') # No interaction is significant
drop1(mod3,test='LRT') # only actgrp is insignificant

# choose model 3 as final model
mod.final=mod3
summary(mod.final)
source("http://www.ics.uci.edu/~staceyah/111-202/R/Stats111-202-Rfunctions.R")
glmCI(mod.final)

########################## part (c) ##########################
###################### Residual Diagnositics #####################

## Evaluate variance assumption
presids = residuals(mod.final,type="pearson")
muhat = fitted(mod.final)

# Plot of squared Pearson resids vs. fitted to check variance assumption
par(mfrow=c(1,2))
plot(muhat, presids^2, xlab="Fitted Probability", ylab="Pearson Residual Squared",
     main="Squared Pearson Residuals")
sfit = supsmu(muhat, presids^2)
lines(sfit$x[order(sfit$x)], sfit$y[order(sfit$x)], col="red", lwd=2)

identify(muhat, presids^2) #1132, check this point
outdoor[1132,] #not marginally extrem point
#plot without this point
plot(muhat[-1132],presids[-1132]^2,xlab="Fitted Probability", 
     ylab="Pearson Residual Squared",main="Squared Pearson Residuals")
sfit=supsmu(muhat[-1132],presids[-1132]^2)
lines(sfit$x[order(sfit$x)], sfit$y[order(sfit$x)], col="red", lwd=2)

## Identify influential points
# Plot of Cook's Distance vs. leverage to identify influential points
par(mfrow=c(1,1))
plot(hatvalues(mod.final), cooks.distance(mod.final), xlab="Leverage",
     ylab="Cook's Distance")

# p:
length(coef(mod.final))  # 6
# n-p:
summary(mod.final)$df.residual # 1559
# Cutoff for high Cook's distance:
cd.cut=qf(0.5,6,1559)  # 0.8917; no point with a Cook's distance higher than this
# Any points with high leverage?
abline(v=3*mean(hatvalues(mod.final)),col='red',lwd=2)
identify(hatvalues(mod.final), cooks.distance(mod.final)) 
# 24  370  448  567  640  667 1017 1137 1140 1179 1340 1407 1445
infl.points=c(24,370,448,567,640,667,1017,1137,1140,1179,1340,1407,1445)
# Should look at these observations
# What R thinks are influential observations
infl = influence.measures(mod.final) 
infl.points.r=which(infl$is.inf[,10]==T)
infl.points.r
# 5 more points other than we identified before: 235, 422, 872, 1252,1493

#remove these points
outdoor=outdoor[-infl.points, ]
attach(outdoor)

#run previous code with the new data, the model selection process remains the same

########################## part (d) ##########################
####################### Goodness of Fit Test ######################

# Hosmer-Lemshow GOF:
binary.gof(mod.final) 
# Expected counts with default number of groups too low (4 below 5)
binary.gof(mod.final,ngrp=5)  # Now only one expected count below 5
# p-value = 0.9252, no evidence of lack of fit
