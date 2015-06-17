library(MASS)
library(e1071)
library(mi)
library(pROC)
library(boot)
library(randomForest)
library(rpart)
library(party)

data=read.table('D:/study/235/crx.txt',sep=',')

# first dealing with missing values
data[data=='?']=NA
#data.new=na.omit(data) abandon missing values

mi.info(data)
data.new= random.imp(data, imp.method = c( "bootstrap", "pca" ) )

#turn V2,V14 into continuous
data.new$V2=as.numeric(data.new$V2)
data.new$V14=as.numeric(data.new$V14)

n = dim(data.new)[1]
p = dim(data.new)[2]
colnames(data.new)[p]<-'y'

a=c(2,3,8,11,14,15,16) # 2,3,8,11,14,15 continuous
data.co=data.new[,a] # data with only continuous predictors

####
#### descriptive statistics
####

attach(data.new)
table(y)

# quant. predictors stratified by normal/low bwt:
x=rbind(tapply(V2,y,mean),tapply(V3,y,mean),tapply(V8,y,mean),tapply(V11,y,mean),tapply(V14,y,mean),tapply(V15,y,mean))
y=rbind(tapply(V2,y,sd),tapply(V3,y,sd),tapply(V8,y,sd),tapply(V11,y,sd),tapply(V14,y,sd),tapply(V15,y,sd))
z=cbind(x,y)
write.csv(z,file='D:/study/finalsummary.csv')
tapply(V3,y,mean),
tapply(V3,y,sd),
tapply(V8,y,mean),
tapply(V8,y,sd)
tapply(V11,y,mean)
tapply(V11,y,sd)
tapply(V14,y,mean)
tapply(V14,y,sd)
tapply(V15,y,mean)
tapply(V15,y,sd)

# Summary of cat. predictors:
t1=table(V1,y)
cbind(t1[,2],t1[,2]/rowSums(t1))
t4=table(V4,y)
cbind(t4[,2],t4[,2]/rowSums(t4))
t5=table(V5,y)
cbind(t5[,2],t5[,2]/rowSums(t5))
t6=table(V6,y)
cbind(t6[,2],t6[,2]/rowSums(t6))
t7=table(V7,y)
cbind(t7[,2],t7[,2]/rowSums(t7))
t9=table(V9,y)
cbind(t9[,2],t9[,2]/rowSums(t9))
t10=table(V10,y)
cbind(t10[,2],t10[,2]/rowSums(t10))
t12=table(V12,y)
cbind(t12[,2],t12[,2]/rowSums(t12))
t13=table(V13,y)
cbind(t13[,2],t13[,2]/rowSums(t13))

# Mainly interested in how predictors are related to actgrp.
# Side-by-side boxplots with quantitative predictors:
par(mfrow=c(2,3))
boxplot(V2~factor(y,levels=c('-','+'),labels=c("Denied","Approved")),ylab="V2",
        main="V2 vs. Credit Card Approval")
boxplot(V3~factor(y,levels=c('-','+'),labels=c("Denied","Approved")),ylab="V3",
       main="V3 vs. Credit Card Approval")
boxplot(V8~factor(y,levels=c('-','+'),labels=c("Denied","Approved")),ylab="V8",
        main="V8 vs. Credit Card Approval")
boxplot(V11~factor(y,levels=c('-','+'),labels=c("Denied","Approved")),ylab="V11",
        main="V11 vs. Credit Card Approval")
boxplot(V14~factor(y,levels=c('-','+'),labels=c("Denied","Approved")),ylab="V14",
        main="V14 vs. Credit Card Approval")
boxplot(V15~factor(y,levels=c('-','+'),labels=c("Denied","Approved")),ylab="V15",
        main="V15 vs. Credit Card Approval")


####
#### classification
####

nIter=100

mean.log.co=rep(0,nIter)
mean.log=rep(0,nIter)
mean.lda=rep(0,nIter)
mean.qda=rep(0,nIter)
mean.nb=rep(0,nIter)
mean.svm.co.l=rep(0,nIter)
mean.svm.co.r=rep(0,nIter)
mean.svm.l=rep(0,nIter)
mean.svm.r=rep(0,nIter)
acc=rep(0,nIter)
acc.co=rep(0,nIter)

for(i in 1:nIter){
  
  ### Divid the data into training (2/3) and test (1/3). 
  ind.tr <- sample(n, floor(2*n/3))
  ind.te <- setdiff(seq(1, n), ind.tr)

####
#### Logistic Regression with continuous predictors
####
  logist.co=glm(y~ .,family = binomial("logit"),data = data.co[ind.tr,])
  pred.class <- predict(logist, newdata=data.co[ind.te, ],type='response')
  pred.class[pred.class>0.5]=1
  pred.class[pred.class<0.5]=0
  pred.class=factor(pred.class,levels=c(0,1),labels=c('-','+'))
  act.class <- data.co$y[ind.te]
  mean.log.co[i]=mean(act.class == pred.class)
  
  
####
#### LDA, QDA, NB with continuous predictors ####
####



lda.m <- lda(y~., data=data.co[ind.tr, ])
pred.class <- predict(lda.m, data.co[ind.te, ])$class
act.class <- data.co$y[ind.te]
mean.lda[i]=mean(act.class == pred.class)

qda.m <- qda(y~., data=data.co[ind.tr, ])
pred.class <- predict(qda.m, data.co[ind.te, ])$class
act.class <- data.co$y[ind.te]
mean.qda[i]=mean(act.class == pred.class)

nb.m <- naiveBayes(y~., data=data.co, subset=ind.tr)
pred.class <- predict(nb.m, newdata=data.co[ind.te, -p])
act.class <- data.co$y[ind.te]
mean.nb[i]=mean(act.class == pred.class)

####
#### svm with continuous predictors
####
svm.co.l <- svm(y~., kernel='linear', data=data.co[ind.tr, ])
pred.class <- predict(svm.m, newdata=data.co[ind.te, ])
act.class <- data.co$y[ind.te]
mean.svm.co.l[i]=mean(act.class == pred.class)

svm.co.r <- svm(y~., kernel='radial', data=data.co[ind.tr, ])
pred.class <- predict(svm.m, newdata=data.co[ind.te, ])
act.class <- data.co$y[ind.te]
mean.svm.co.r[i]=mean(act.class == pred.class)

####
#### svm with all predictors
####

svm.l <- svm(y~., kernel='linear', data=data.new[ind.tr, ])
pred.class <- predict(svm.m, newdata=data.new[ind.te, ])
act.class <- data.new$y[ind.te]
mean.svm.l[i]=mean(act.class == pred.class)

svm.r <- svm(y~., kernel='radial', data=data.new[ind.tr, ])
pred.class <- predict(svm.m, newdata=data.new[ind.te, ])
act.class <- data.new$y[ind.te]
mean.svm.r[i]=mean(act.class == pred.class)

####
#### Random Forest
####

fit.rf.co <- randomForest(y~., data=data.co[ind.tr, ], ntree=1000)
y.hat <- predict(fit.rf, newdata=data.co[ind.te, ], type="response")
tbl<-table(data.co$y[ind.te], y.hat)
acc.co[i] <- sum(diag(tbl))/sum(tbl)

fit.rf <- randomForest(y~., data=data.new[ind.tr, ], ntree=1000)
y.hat <- predict(fit.rf, newdata=data.new[ind.te, ], type="response")
tbl<-table(data.new$y[ind.te], y.hat)
acc[i] <- sum(diag(tbl))/sum(tbl)

}

colMeans(cbind(mean.log.co,mean.lda, mean.qda, mean.nb,mean.svm.co.l,mean.svm.co.r,mean.svm.l,mean.svm.r,acc.co,acc))
#0.7614783     0.7416087     0.7200000     0.7154783     0.7614348 0.7617826     0.8494783     0.8574783     0.8654348 

####
#### Logistic Regression with all predictors
####
logist=glm(y~ .,family = binomial("logit"),data = data.new[ind.tr,])
pred.class <- predict(logist, newdata=data.new[ind.te, ],type='response')
pred.class[pred.class>0.5]=1
pred.class[pred.class<0.5]=0
pred.class=factor(pred.class,levels=c(0,1),labels=c('-','+'))
act.class <- data.new$y[ind.te]
mean=mean(act.class == pred.class)

plot(svm.l, data.new[ind.te, ], formula= V2~V3)


###### RPART

fit.part <- rpart(y~., data=data.new)

plot(fit.part, uniform=TRUE)
text(fit.part, use.n=TRUE, cex=.6)

pfit<- prune(fit.part, cp=fit.part$cptable[which.min(fit.part$cptable[,"xerror"]),"CP"])
plot(pfit, uniform=TRUE)
text(pfit, use.n=TRUE, all=TRUE, cex=.55)



varImpPlot(fit.rf)

p <- predict(fit.rf, newdata=data.new[ind.te, ], type="prob")[, 2]

roc(data.new$y[ind.te], p, plot=TRUE)

auc <- roc(data$V16[ind.te], p, plot=FALSE)$auc







####
#### ROC curves & auc
####
par(mfrow=c(3,2))
library(ROCR)

y.hat <-predict(logist.co, newdata=data.co[ind.te, ], type="response")
y.hat[y.hat>0.5]=1
y.hat[y.hat<0.5]=0
p <- predict(logist.co, newdata=data.co[ind.te, ], type="response")
roc(y[ind.te], p, main="logistic model with continuous predictors",plot=TRUE)
#0.7861

y.hat <-predict(logist, newdata=data.new[ind.te, ], type="response")
y.hat[y.hat>0.5]=1
y.hat[y.hat<0.5]=0
p <- predict(logist, newdata=data.new[ind.te, ], type="response")
roc(y[ind.te], p, main="logistic model with all predictors",plot=TRUE)
#0.8815

p=predict(lda.m,data.co[ind.te,],type='prob')$posterior[,2]
roc(y[ind.te], p, main="Linear Discriminant Analysis model",plot=TRUE)
#0.7841

p=predict(qda.m,data.co[ind.te,],type='prob')$posterior[,2]
roc(y[ind.te], p, main="Quadratic Discriminant Analysis model",plot=TRUE)
# 0.7808


p=predict(nb.m,data.co[ind.te,],type='raw')[,2]
roc(y[ind.te], p, main="Naive Bayes model",plot=TRUE)
# 0.7317

y.hat=predict(svm.co.l, data.co[ind.te, ],decision.values = TRUE) 
p=attributes(y.hat)$decision.values
roc(y[ind.te], p, main="SVM model with linear kernel with all predictors",plot=TRUE)
# 0.79

y.hat=predict(svm.co.r, data.co[ind.te, ],decision.values = TRUE) 
p=attributes(y.hat)$decision.values
roc(y[ind.te], p, main="SVM model with radial kernel with continuous predictors",plot=TRUE)
#0.8481

y.hat=predict(svm.l, data.new[ind.te, ],decision.values = TRUE) 
p=attributes(y.hat)$decision.values
roc(y[ind.te], p, main="SVM model with linear kernel with all predictors",plot=TRUE)
#0.9402

y.hat=predict(svm.r, data.new[ind.te, ],decision.values = TRUE) 
p=attributes(y.hat)$decision.values
roc(y[ind.te], p, main="SVM model with radial kernel with all predictors",plot=TRUE)
#0.9468

p=predict(fit.rf.co, newdata=data.co[ind.te, ], type="prob")[, 2]
roc(y[ind.te], p, main='Random Forest model with continuous predictors',plot=TRUE)
#0.9613

p=predict(fit.rf, newdata=data.new[ind.te, ], type="prob")[, 2]
roc(y[ind.te], p, main='Random Forest model all predictors',plot=TRUE)
# 0.987


varImpPlot(fit.rf, main='variable importance of random forest model')
