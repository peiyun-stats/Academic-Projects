#3
library(MASS)
library(lasso2)
library(pls)
library(glmnet)

fire=read.csv('D:/study/235/forestfires.csv')
head(fire)
fire$month=factor(fire$month)
fire$day=factor(fire$day)
attach(fire)
dim(fire)
fire[, 1:12]=scale(fire[, 1:12])
x=as.matrix(fire[, 1:12])
y=fire[,13]
n=length(y)

#(a) M1: Multiple linear regression model that includes all the predictors

M1=lm(area~X+Y+month+day+FFMC+DMC+DC+ISI+temp+RH+wind+rain)
summary(M1)
m0=lm(area~X+Y+month+day+FFMC+DMC+DC+ISI+temp+RH+wind+rain,data=fire)

#(b) M2: Principal component regression
dat=cbind(X,Y,month,day,FFMC,DMC,DC,ISI,temp,RH,wind,rain)

pc=princomp(dat)
z=pc$scores[,1:5]
M2=lm(area ~ z)
x=as.matrix(fire[,1:12])
y=fire[,13]
M02=pcr(y~x,ncomp=5,scale=T,validation="CV",segments=5)


#(c) M3: Partial least squares
M3=plsr(y~x,ncomp=5, data=fire, scale=TRUE)
pcr.mod <- pcr(y ~ x, ncomp=5, scale=TRUE, validation = "CV", segments = 5)




pls.mod <- plsr(y ~ x, ncomp=5, data=Prostate, scale=TRUE)

#(d) M4: Ridge regression
M4=lm.ridge(y~x,lambda=seq(0,10,0.1))
ridge.res <- lm.ridge(y~x, lambda = seq(0, 10, 0.1))
plot(ridge.res)
select(ridge.res)

# Lasso with glmnet

lasso.fit=glmnet(x,y)
plot(lasso.fit)

y.hat <- predict(lasso.fit, x, s=0.01, type='response')

#(e) M5: Lasso
M5=glmnet(x,y)




# We could have the data-splitting method to choose the number of principal components



# Here, we use 5-fold CV to evaluate the performance of the PCR with 5 components

# To be even more precise, we should put this within the loop too
pc <- princomp(x)
z<- pc$scores[, 1:5] # Choosing the first 5 principal components

ind <- sample(5, n, replace=TRUE)


mse = rep(0, 5)
for(i in 1:5){

	pcr.mod <- lm(y[ind!=i] ~ z[ind!=i, ])

	y.hat <- as.matrix(cbind(1, z[ind==i, ])) %*% matrix(coef(pcr.mod))

	mse[i] = mean((y[ind==i] - y.hat)^2)


}

mean(mse)



# we can use the pcr function in the pls package




# We can use the "validation" for fine-tuning parameters. See the help document for plsr.

# Ridge regression

ridge.res <- lm.ridge(y~x, lambda = seq(0, 10, 0.1))
plot(ridge.res)
select(ridge.res)

# Lasso with glmnet

lasso.fit=glmnet(x,y)
plot(lasso.fit)

y.hat <- predict(lasso.fit, x, s=0.01, type='response')
