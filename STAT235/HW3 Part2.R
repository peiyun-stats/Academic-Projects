library(splines)
data = as.matrix(read.table("D:/study/235/splineExample.txt",header=TRUE))

############## a #################

x1=data[,1]
x2=data[,2]
y=data[,3]

##cross validation
k_fold=function(n=100,K=5){
  rs=runif(n)
  id=order(rs)
  k=as.integer(n*seq(1,K-1)/K)
  k=matrix(c(0,rep(k,each=2),n),ncol=2,byrow=TRUE)
  k[,1]=k[,1]+1
  l=lapply(seq.int(K),function(x,k,d) 
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}
n=length(y)
cv=k_fold(n)

#identify knots
knots1=quantile(x1,probs=seq(0,1,length.out=4))[c(2:3)]
knots2=quantile(x1,probs=seq(0,1,length.out=5))[c(2:4)]
knots3=quantile(x1,probs=seq(0,1,length.out=6))[c(2:5)]

#2 knots
sum1=0
for(i in 1:5){
  basis.fn=bs(x1[cv[[i]]$train], knots=knots1, degree=1,Boundary.knots=range(x1))
  B.model=lm(y[cv[[i]]$train]~basis.fn)
  x1.new=x1[cv[[i]]$test]
  basis.fn.x1=predict(basis.fn,x1.new)
  y.hat=cbind(1,basis.fn.x1)%*%coef(B.model)
  meanEach=mean((y[cv[[i]]$test]-y.hat)^2)
  sum1=sum1+meanEach
}
meanError1=sum1/5

#3 knots
sum2=0
for(i in 1:5){
  basis.fn=bs(x1[cv[[i]]$train], knots=knots2, degree=1,Boundary.knots=range(x1))
  B.model=lm(y[cv[[i]]$train]~basis.fn)
  x1.new=x1[cv[[i]]$test]
  basis.fn.x1=predict(basis.fn,x1.new)
  y.hat=cbind(1,basis.fn.x1)%*%coef(B.model)
  meanEach=mean((y[cv[[i]]$test]-y.hat)^2)
  sum2=sum2+meanEach
}
meanError2=sum2/5

#4 knots
sum3=0
for(i in 1:5){
  basis.fn=bs(x1[cv[[i]]$train], knots=knots3, degree=1,Boundary.knots=range(x1))
  B.model=lm(y[cv[[i]]$train]~basis.fn)
  x1.new=x1[cv[[i]]$test]
  basis.fn.x1=predict(basis.fn,x1.new)
  y.hat=cbind(1,basis.fn.x1)%*%coef(B.model)
  meanEach=mean((y[cv[[i]]$test]-y.hat)^2)
  sum3=sum3+meanEach
}
meanError3=sum3/5

#result
cbind(meanError1,meanError2,meanError3)
#1.706581  0.4488537  0.8905409
#we should choose 3 knots.

## plot of 3 knots
plot(y~x1)
abline(v=knots2[1], col='red')
abline(v=knots2[2], col='red')
abline(v=knots2[3],col='red')
basis.fn=bs(x1, knots=knots2, degree=1,Boundary.knots=range(x1))
B.model=lm(y~basis.fn)
beta = coef(B.model)

t=seq(min(x1), max(x1), .01)
basis.fn.t=predict(basis.fn, t)
y.hat=cbind(1, basis.fn.t)%*%beta
lines(t, y.hat, lwd=2)





####### b #############
plot(y~x2)
knots0=quantile(x2,probs=seq(0,1,length.out=5))[c(2:4)]
basis.fn=ns(x2, knots=knots0)
beta=coef(lm(y~basis.fn))

t=seq(min(x2), max(x2), .1)
basis.fn.t = predict(basis.fn, t)
y.hat=cbind(1, basis.fn.t)%*%beta
lines(t, y.hat, lwd=2)

########### c #########

sum=rep(0,50)
for(j in 1:50){
  for(i in 1:5){
    ss.model=smooth.spline(x1[cv[[i]]$train], y[cv[[i]]$train], df=j+1)
    x1.new=x1[cv[[i]]$test]
    y.hat=predict(ss.model,x1.new)$y
    meanEach=mean((y[cv[[i]]$test]-y.hat)^2)
    sum[j]=sum[j]+meanEach
  }
}
meanError=sum/5
which.min(meanError) #11, df=12

plot(y~x1)
ss.model=smooth.spline(x1,y,df=12)
lines(ss.model, lwd=2)

########### d ##########
basis.fn.x1=ns(x1, knots=knots2)
basis.fn.x2=ns(x2, knots=knots0)
beta = coef(lm(y~basis.fn.x1+basis.fn.x2))

t1=seq(min(x1), max(x1), length.out=500)
basis.fn.t1 = predict(basis.fn.x1, t1)
t2=seq(min(x2),max(x2), length.out=500)
basis.fn.t2 = predict(basis.fn.x2, t2)
y.hat=cbind(1, basis.fn.t1,basis.fn.t2)%*%beta

############ e ############
library(gam)
gam1=gam(y~s(x1)+s(x2))

par(ask=T)
par(mfrow=c(1,2))
plot(gam1, se=TRUE)

summary(gam1)

########### f ##############
grp.train=data[1:60,]
grp.validation=data[61:80,]
grp.test=data[81:100,]

x1.train=grp.train[,1]
x2.train=grp.train[,2]
y.train=grp.train[,3]

x1.validation=grp.validation[,1]
x2.validation=grp.validation[,2]
y.validation=grp.validation[,3]

x1.test=grp.test[,1]
x2.test=grp.test[,2]
y.test=grp.test[,3]

########## M1 ###########
basis.fn.x1=ns(x1.train, knots=knots2)
basis.fn.x2=ns(x2.train, knots=knots0)
beta = coef(lm(y.train~basis.fn.x1+basis.fn.x2))

basis.fn.x1.validation=predict(basis.fn.x1,x1.validation)
basis.fn.x2.validation=predict(basis.fn.x2,x2.validation)
y.hat.validation=cbind(1,basis.fn.x1.validation,basis.fn.x2.validation)%*%beta

basis.fn.x1.test=predict(basis.fn.x1,x1.test)
basis.fn.x2.test=predict(basis.fn.x2,x2.test)
y.hat.test=cbind(1,basis.fn.x1.test,basis.fn.x2.test)%*%beta


MSE1.valid=mean((y.hat.validation-y.validation)^2)
MSE1.test=mean((y.hat.test-y.test)^2)

######### M2 ###########
gam.train=gam(y.train~s(x1.train)+s(x2.train))
y.hat.validation=cbind(1,x1.validation,x2.validation)%*%coef(gam.train)
y.hat.test=cbind(1,x1.test,x2.test)%*%coef(gam.train)
MSE2.valid=mean((y.hat.validation-y.validation)^2)
MSE2.test=mean((y.hat.test-y.test)^2)
