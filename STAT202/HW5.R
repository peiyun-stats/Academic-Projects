source("http://www.ics.uci.edu/~staceyah/111-202/R/Stats111-202-Rfunctions.R")
bwt = read.csv("http://www.ics.uci.edu/~staceyah/111-202/data/bwt.csv",header=TRUE)
attach(bwt)
summary(bwt)
plot(bwt)

cor(cbind(low,age,lwt,as.numeric(race),smoke,ptd,ht,ui,as.numeric(ftv)))
#not highly correlated

#1. Which factors contribute to low birth weight?
#2. Does the mother's smoking status have an efect on the probability of low birth weightfor her child?

chisq.test(table(low,smoke),correct=F) #No


huge.mod=glm(low~.^2,data=bwt,family=binomial)
summary(huge.mod)

anova(huge.mod, test="LRT")
drop1(huge.mod, test="LRT")

mod1=glm(low~.,data=bwt,family=binomial)
summary(mod1)
anova(mod1,huge.mod, test="LRT")

mod2=glm(low~ptd*ht,family=binomial)
summary(mod2)
anova(mod2,huge.mod,test="LRT")

mod3=glm(low~ptd+ht,family=binomial)
anova(mod3,mod2,test='LRT') #interaction is not necessary
add1(mod3,~.+age,test='LRT')#Yes
mod4=glm(low~age+ptd+ht,family=binomial)
add1(mod4,~.+lwt,test='LRT')#Yes
mod5=glm(low~age+lwt+ptd+ht,family=binomial)
add1(mod5,~.+race,test='LRT')#No
add1(mod5,~.+smoke,test='LRT')#No
add1(mod5,~.+ui,test='LRT')#no
add1(mod5,~.+ftv,test='LRT')#no

summary(glm(low~age*lwt,family=binomial)) #No
summary(glm(low~age*ptd,family=binomial)) #No
summary(glm(low~age*ht,family=binomial)) #No

summary(glm(low~lwt*ptd,family=binomial))
summary(glm(low~lwt*ht,family=binomial))

summary(mod5)
drop1(mod5,test='LRT')
mod6=glm(low~lwt+ptd+ht,family=binomial)
anova(mod6,mod5,test='LRT')

mod7=glm(low~log(lwt)+ptd+ht,family=binomial) #not much difference in AIC
AIC(mod6)
AIC(mod7) #not much difference in AIC

l=lwt
l2=lwt^2
l3=lwt^3

mod8=glm(low~l+l2+l3+ptd+ht,family=binomial)
anova(mod8,test='LRT')

binary.gof(mod6)   #GOF

presids=residuals(mod6,type='pearson')
muhat=fitted(mod6)
plot(muhat,presids^2,xlab="Fitted probability",ylab="Pearson Residual Squared")

## Add smoother:
sfit <- supsmu( muhat, presids^2 )
lines( sfit$x[ order( sfit$x ) ], sfit$y[ order( sfit$x ) ],, col="red", lwd=2 )



library(ROCR)
probs=fitted(mod6)
pred.obj = prediction(probs,low)
# Plot ROC curve
plot(performance(pred.obj,"tpr","fpr"))
# Area under the ROC curve
performance(pred.obj,"auc")

### Influential observations?
lev = hatvalues(mod6)
numb=c(1:length(low))
plot(numb,lev)
abline(h=3*mean(lev))
text(numb-2,lev,numb)

cook = cooks.distance(mod6)
plot(lev,cook,xlab='leverage',ylab="Cooks Distance",ylim=c(0,.3))
abline(h=qf(.05,length(mod6$coef),mod6$df.residual),col="blue",lwd=2)
abline(v=3*mean(lev),col='red',lwd=2)
text(lev-.02,cook,numb)
dfb = dfbeta(mod6)

par(mfrow=c(2,2))
for( i in 1:4){
  plot( numb, dfbeta(mod6)[,i], xlab="Patient ID",ylab="Delta Beta", main=dimnames(dfbeta(mod6))[[2]][i] )
  text( numb-2, dfbeta(mod6)[,i], numb )
}
