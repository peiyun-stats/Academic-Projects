hpd.gam = function (alpha, beta=1, p=.95, max.iter=100, small=1e-8){
  
  if (alpha <= 1)
    stop("alpha must be greater than one")
  
  # initialize
  tail=1-p
  l1=0
  l2=qgamma(tail, alpha, beta)
  for (iter in 1:max.iter){
    l=0.5*(l1+l2)
    lprob=pgamma(l,alpha,beta)
    ldens=dgamma(l,alpha,beta)
    uprob=lprob+p
    u=qgamma(uprob,alpha,beta)
    udens=dgamma(u,alpha,beta)
    if ((udens-ldens) < 0){
      l2=l
    }
    else {
      l1=l      
    }
    if ( l2 - l1 < small)
      break
  }
  intv = c(l,u)
  out = list (intv, iter)
  names(out) = c("HPD interval", "iteration.count")
  out
}


# (b)
x=seq(2,10,.01)
y=rep(NA,length(x))
z=rep(NA,length(x))
for (i in 1 : length(x)){
  y[i] = hpd.gam(x[i])$HPD[1]
  z[i] = hpd.gam(x[i])$HPD[2]
}
plot(z~x, ylim = c(0,20),xlab = 'Alpha', ylab = 'HPD Interval', type = 'l', lty = 2)
lines(y~x, lty = 1)
legend(2, 20, c('CI.low','CI.high'),lty=c(1,2))

# (c) modified for alpha < 1
hpd.gam2 = function (alpha, beta=1, p=.95, max.iter=100, small=1e-8){
  
  if (alpha <= 0)
    stop("alpha must be greater than zero")
  
  if (alpha > 0 & alpha <= 1){
    l=0
    u=qgamma(p,alpha,beta)
    iter=1
  }
  
  else{
  
  # initialize
  tail=1-p
  l1=0
  l2=qgamma(tail, alpha, beta)
  for (iter in 1:max.iter){
    l=0.5*(l1+l2)
    lprob=pgamma(l,alpha,beta)
    ldens=dgamma(l,alpha,beta)
    uprob=lprob+p
    u=qgamma(uprob,alpha,beta)
    udens=dgamma(u,alpha,beta)
    if ((udens-ldens) < 0){
      l2=l
    }
    else {
      l1=l      
    }
    if ( l2 - l1 < small)
      break
  }
  
  }
  intv = c(l,u)
  out = list (intv, iter)
  names(out) = c("HPD interval", "iteration.count")
  out
}



### 3 ###
## (a) ##
# first derivative of log-likelihood function
p = function(theta,a,b) { 
  sum((a-b)/((a-b)*theta+b))
}

# derivative of p
deri = function(theta,a,b) {  
  sum(-((a-b)/((a-b)*theta+b))^2)
}

MLE.newton = function(theta=0, y, mu1, mu2, sigma1, sigma2,may.iter=1000, small=1e-5) { 
  
  a = dnorm(y, mu1, sigma1)
  b = dnorm(y, mu2, sigma2)
  
  for (iter in 1:may.iter) {
    theta.old = theta
    theta = theta.old - p(theta.old,a,b) / deri(theta.old,a,b) 
    if (abs(theta - theta.old) < small) 
      break
  }
  out <- list(theta, iter, iter != may.iter)
  names(out) <- c("MLE", "iteration.count", "converge")
  out
}

# (b)
y = c(8.1, 8.2, 8.1, 8.2, 8.2, 7.4, 7.3, 7.4, 8.1, 8.1, 7.9, 7.8, 8.2, 7.9, 7.9, 8.1, 8.1)
x=seq(0,1,.001) # start points
z=rep(NA,length(x)) # MLE
conv=rep(NA,length(x))  # convergence
for(i in 1:length(x)){
  z[i]=MLE.newton(theta=x[i],y=y,mu1=7.5,mu2=8.0,sigma1=sqrt(.1),sigma2=sqrt(.1))$MLE
  conv[i]=MLE.newton(theta=x[i],y=y,mu1=7.5,mu2=8.0,sigma1=sqrt(.1),sigma2=sqrt(.1))$conv
}
sum(conv-1) # zero
z[1] # MLE = 0.09625715

# (c)
# first derivative of log-likelihood function for delta
p2 = function(delta,a,b,n) { 
  sum(a*exp(delta)/(a*exp(delta)+b))-n*exp(delta)/(1+exp(delta))
}

# derivative of p
deri2 = function(delta,a,b,n) {  
  sum(a*b*exp(delta)/((a*exp(delta)+b)^2)-n*exp(delta)/((1+exp(delta))^2))
}

MLE.newton2 = function(delta=0, y, mu1, mu2, sigma1, sigma2,may.iter=1000, small=1e-5) { 
  
  
  a = dnorm(y, mu1, sigma1)
  b = dnorm(y, mu2, sigma2)
  n = length(y)
  
  for (iter in 1:may.iter) {
    delta.old = delta
    delta = delta.old - p2(delta.old,a,b,n) / deri2(delta.old,a,b,n) 
    if (abs(delta - delta.old) < small) 
      break
  }
  out <- list(delta, iter, iter != may.iter)
  names(out) <- c("MLE", "iteration.count", "converge")
  out
}

# (d)
z2=rep(NA,length(x)) # MLE
conv2=rep(NA,length(x))  # convergence
for(i in 1:length(x)){
  z2[i]=MLE.newton2(delta=x[i],y=y,mu1=7.5,mu2=8.0,sigma1=sqrt(.1),sigma2=sqrt(.1))$MLE
  conv2[i]=MLE.newton2(delta=x[i],y=y,mu1=7.5,mu2=8.0,sigma1=sqrt(.1),sigma2=sqrt(.1))$conv
}
sum(conv2-1) # zero
z2[1] # MLE = -2.238935
