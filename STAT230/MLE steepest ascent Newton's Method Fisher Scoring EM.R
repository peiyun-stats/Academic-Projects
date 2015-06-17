1.  (c) 
Y=c(24,25,31,31,22,21,26,20,16,22)
a=c(734,516,754,877,814,362,764,809,223,1066)/10000
b=c(0.19,0.12,0.15,0.16,0.14,0.06,0.13,0.13,0.03,0.15)
X=a/b

# compute MLE of (alpha, beta)
# i. steepest ascent
SteepestAscent = function (X, Y, theta, delta, max.iter=10000, small=1e-8){
  for (iter in 1:max.iter){
    theta.old = theta
    alpha = theta.old[1]
    beta = theta.old[2]
    l1=c(sum(Y/(alpha+beta*X)-1),sum(X*Y/(alpha+beta*X)-X))
    theta=theta.old+1/delta*l1
    
    if (sum(abs(theta-theta.old))<small)
      break
  }
  out = list (theta, iter)
  names(out) = c("MLE","iteration.count")
  out
}

# ii. Newton's Method
NewtonsMethod = function (X, Y, theta, max.iter=10000, small=1e-8){
  for (iter in 1:max.iter){
    theta.old = theta
    alpha = theta.old[1]
    beta = theta.old[2]
    l1=c(sum(Y/(alpha+beta*X)-1),sum(X*Y/(alpha+beta*X)-X))
    J=matrix(c(sum(1/(alpha+beta*X)),sum(X/(alpha+beta*X)),sum(X/(alpha+beta*X)),sum((X^2)/(alpha+beta*X))),2,2)
    theta=theta.old+solve(J)%*%l1
    
    if (sum(abs(theta-theta.old))<small)
      break
  }
  out = list (theta, iter)
  names(out) = c("MLE","iteration.count")
  out
}

# iii. Fisher scoring
FisherScoring = function (X, Y, theta, max.iter=10000, small=1e-8){
  for (iter in 1:max.iter){
    theta.old = theta
    alpha = theta.old[1]
    beta = theta.old[2]
    l1=c(sum(Y/(alpha+beta*X)-1),sum(X*Y/(alpha+beta*X)-X))
    l2=-matrix(c(sum(Y/(alpha+beta*X)^2),sum(X*Y/(alpha+beta*X)^2),sum(X*Y/(alpha+beta*X)^2),sum((X^2)*Y/(alpha+beta*X)^2)),2,2)
    theta=theta.old+solve(-l2)%*%l1
    
    if (sum(abs(theta-theta.old))<small)
      break
  }
  out = list (theta, iter)
  names(out) = c("MLE","iteration.count")
  out
}

(d)
# starting value
theta0=lm(Y~X,weights=1/Y)$coef
theta0
# (Intercept)           X 
# 37.83374   -25.44075


(e) 
# compare
SteepestAscent(X=X, Y=Y, theta = theta, delta=.5)
$MLE
# (Intercept)           X 
# 38.98629   -26.56865 
$iteration.count
# [1] 2180

NewtonsMethod(X=X,Y=Y,theta=theta0)
$MLE
#[,1]
#[1,]  38.98630
#[2,] -26.56865
$iteration.count
#[1] 9

FisherScoring(X=X,Y=Y,theta=theta0)
$MLE
#[,1]
#[1,]  38.98630
#[2,] -26.56865
$iteration.count
#[1] 4


#Using starting value in part (d), the three methods all converge to (38.986, 0-26.569). For Steepest Ascent algorithm, 0.25, 0.5, 1, 2, 4 are used as steep size and 0.5 works the best.
#For the starting value we chose, Fisher scoring and Newton’s method converge very fast, and steepest ascent converges really slowly, even with a “good” step size.


2.
EM = function(delta,y) { 
  # E step
  a = delta[1]*dnorm(y, delta[2],sqrt(delta[4]))
  b = (1-delta[1])*dnorm(y, delta[3], sqrt(delta[5]))
  P1=a/(a+b)
  P0=b/(a+b)
  
  # M-step
  theta = mean(P1)
  mu1 = sum(y*P1)/sum(P1)
  mu2 = sum(y*P0)/sum(P0)
  sigma1 = sum(P1*(y-sum(y*P1)/sum(P1))^2)/sum(P1)
  sigma2 = sum(P0*(y-sum(y*P0)/sum(P0))^2)/sum(P0)
  c(theta, mu1, mu2, sigma1, sigma2)
}

Y = c(8.1,8.2,8.1,8.2,8.2,7.4,7.3,7.4,8.1,8.1,7.9,7.8,8.2,7.9,7.9,8.1,
      8.1)
delta = c(.5, 5, 10, 1, 2)
conv = rep (NA, 1000)
for (iter in 1:1000){
  delta.old = delta
  delta = EM (delta.old,Y)
  conv[iter] = delta[1]
}
delta
#[1] 0.176470153 7.366666588 8.064285362 0.002222225 0.016581858

plot(conv[1:100]~c(1:100), xlab = "iteration count", ylab = "theta")

# This EM algorithm doesn’t always converge. There are more than one mode.
