library(nlme)

contact = read.csv('d:study/203/FinalVA.csv')
# convert to long form
con.l = reshape(contact, idvar='id',varying=list(5:8), v.names = 'y',direction='long')
# redefine as grouped data
con.grp = groupedData(y~time | id, data=con.l)
# add age group variable
con.grp$age.grp = cut(con.grp$age,breaks=c(0,9,13,37,53))

###
##### descriptive statisitcs & plots
###
summary(contact)

plot(con.grp)

par(mfrow=c(1,3))
means1 = tapply(con.grp$y, list(con.grp$time, con.grp$gender),mean)
matplot(1:4,means1,col=c('red','blue'),lty=c(1,3),type='b',pch=c(15,16),
        xlab='Time',ylab='Mean Contact Time (seconds)',
        main='Mean Contact Time by Gender')
legend(1,4.5,c('Female','Male'), col=c('red','blue'),lty=c(1,3),pch=c(15,16))

means2 = tapply(con.grp$y, list(con.grp$time, con.grp$shape),mean)
matplot(1:4,means2,col=c('red','blue'),lty=c(1,3),type='b',pch=c(15,16),
        xlab='Time',ylab='Mean Contact Time (seconds)',
        main='Mean Contact Time by Shape')
legend(1,4.5,c('Box','Circle'), col=c('red','blue'),lty=c(1,3),pch=c(15,16))

means3 = tapply(con.grp$y, list(con.grp$time, con.grp$age.grp),mean)
matplot(1:4,means3,col=c('red','blue','green','black'),lty=1:4,type='b',pch=15:18,
        xlab='Time',ylab='Mean Contact Time (seconds)',
        main='Mean Contact Time by Age')
legend(1,6.5,c('0<Age<=9','9<Age<=13','13<Age<=37','37<Age<=52'),
       col=c('red','blue','green','black'),lty=1:4,pch=15:18)
means1
means2
means3

# plot contact time versus age in each trial
par(mfrow=c(2,2))
plot(contact$trial1~contact$age,xlab='age',ylab='contact time',main='trial1')
plot(contact$trial2~contact$age,xlab='age',ylab='contact time',main='trial2')
plot(contact$trial3~contact$age,xlab='age',ylab='contact time',main='trial3')
plot(contact$trial4~contact$age,xlab='age',ylab='contact time',main='trial4')

###
##### model building
###

# test individual variables on their own
mod.age = lme(y ~ age, random = ~1|id, method = 'ML', data=con.grp)
summary(mod.age) # p-value=0, significant

mod.gender = lme(y ~ gender, random = ~1|id, method = 'ML', data=con.grp)
summary(mod.gender) # p-value=0.0149, significant

mod.shape = lme(y ~ shape, random = ~1|id, method = 'ML', data=con.grp)
summary(mod.shape) # p-value=0.1154, not significant

mod.time = lme(y ~ time, random = ~ 1| id, method = 'ML', data=con.grp)
summary(mod.time) # p-value=0, significant

# model with age, gender and time
mod1 = lme(y ~ age+gender+time, random = ~1|id, method = 'ML', data=con.grp)
summary(mod1) # all significant

# add shape?
mod2 = update(mod1, .~.+shape)
summary(mod2) # significant
anova(mod1, mod2) # p-value=8e-4, yes

# two-way interaction
mod.2int = update(mod2, .~.^2)
summary(mod.2int) # only gender*shape, time*shape not significant, delete?
mod.2int1 = update(mod.2int, .~.-gender*shape-time*shape)
anova(mod.2int, mod.2int1) # p-value = 0.1926, delete
summary(mod.2int1) # all interactions significant
# predictors in interaction should be included
mod3 = update(mod.2int1,.~.+gender+time+shape)

# add one three-way interaction at one time
mod.3int1 = update(mod3, .~.+age*gender*shape)
anova(mod3,mod.3int1) # p-value=0.8743, not significant
mod.3int2 = update(mod3, .~.+age*gender*time)
anova(mod3,mod.3int2) # p-value=9e-4, significant
mod.3int3 = update(mod3, .~.+age*shape*time)
anova(mod3,mod.3int3) # p-value = 0.319, not significant
mod.3int4 = update(mod3, .~.+gender*shape*time)
anova(mod3,mod.3int4) # p-value = 0.5579, not significant
# age*gender*time should be included --> mod.3int2

# add four-way interaction?
mod.4int = update(mod.3int2, .~.+age*gender*shape*time)
anova(mod.3int2, mod.4int) # p-value = 0.6842, not significant

# treat time as categorical?
con.grp$time.f = factor(con.grp$time)
mod.time.f = lme(y ~ time.f, random = ~1|id, method = 'ML', data=con.grp)
summary(mod.time.f) # all p-values=0, marginally significant
mod.f1 = lme(y ~ age+gender+time.f+shape+age*gender+age*time.f+age*shape+gender*time.f+age*gender*time.f, random = ~1|id, method = 'ML', data=con.grp)
summary(mod.f1)$tTable # only age*genderM*time.f4 significant
AIC(mod.3int2) #1269.172
AIC(mod.f1) # 1274.751

# treat age as categorical?
mod.age.f = lme(y ~ age.grp, random = ~1|id, method = 'ML', data=con.grp)
summary(mod.age.f) # all p-value2=0, marginally significant
mod.f2 = lme(y ~ age.grp+gender+time+shape+age.grp*gender+age.grp*time+age.grp*shape+gender*time+age.grp*gender*time, random = ~1|id, method = 'ML', data=con.grp)
summary(mod.f2)$tTable
AIC(mod.f2) # 1215.654

# add age^n?
mod.age2 = update(mod.3int2, .~.+I(age^2))
anova(mod.3int2,mod.age2) # p-value<0.0001 yes

mod.age3 = update(mod.age2, .~.+I(age^3))
anova(mod.age2,mod.age3) # p-value=0.0187 yes

mod.age4 = update(mod.age3, .~.+I(age^4))
anova(mod.age3,mod.age4) # p-value=0.0023 yes

mod.age5 = update(mod.age4, .~.+I(age^5))
anova(mod.age4,mod.age5) # p-value=0.4425 no --> mod.age4

# add random slope?
mod.reml = update(mod.age4,method='REML')
mod.rs = update(mod.reml, random = ~1+time|id)
anova(mod.reml,mod.rs) # p-value=0.1673, no

# choose it as the final model
mod = mod.age4
summary(mod)
1.08299^2 # 1.172867 --> std of random intercept

# save the summary table
a=summary(mod)$tTable
finalmod = cbind(a[,1],a[,1]-1.96*a[,2],a[,1]+1.96*a[,2],a[,5])
write.csv(finalmod,file='D:/study/203/finaltable.csv')


### 
##### Residual Diagnostics
###

plot(mod, id=0.05, adj=-0.4,main='Figure A1: (Standardized) residuals vs. fitted values')
plot(mod, resid(.) ~ time, id=0.05, adj=-0.4, main='Figure A2: (Standardized) residuals vs. time')
plot(mod, id ~ resid(.), abline=0,main='Figure A3: Boxplots of residuals by subject')
qqnorm(mod, id=0.05, adj=-0.4,main='Figure A4: Normal quantile-quantile plot of the residuals')
qqnorm(mod, ~ranef(.), id=0.05, adj=-0.2,main='Figure A5: Normal quantile-quantile plots of the random effects')

