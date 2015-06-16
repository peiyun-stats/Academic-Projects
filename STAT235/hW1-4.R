data=read.csv('D:/study/235/235-roaches.csv')
x=data$senior
y=data$count
n=dim(data)[1]
p=dim(data)[2]

nIter= 20000
burnIn= 2000

beta=matrix(0, nIter, p)
mu0=0
sigma02=10


for(i in 2:nIter) {
  beta0.prime=rnorm(1,beta[i-1,1],1)
  postThetaPrime = getPosterior(c(beta0.prime, beta[i-1, 2]), x, y, beta0.prime, mu0, sigma02)
  postThetaCurrent = getPosterior(c(beta[i-1, 1], beta[i-1, 2]), x, y, beta[i-1, 1], mu0, sigma02)
  acceptanceRatio= min(1, exp(postThetaPrime - postThetaCurrent))

  u = runif(1)
  if (u<=acceptanceRatio) {
		beta[i,1]=beta0.prime}
  else{beta[i,1]=beta[i-1, 1]}
  beta1.prime=rnorm(1, beta[i-1, 2], 1)
      
	postThetaPrime = getPosterior(c(beta[i, 1], beta1.prime), x, y, beta1.prime, mu0, sigma02)
	postThetaCurrent = getPosterior(c(beta[i, 1], beta[i-1, 2]), x, y, beta[i-1, 2], mu0, sigma02)
	acceptanceRatio= min(1, exp(postThetaPrime - postThetaCurrent))

  u = runif(1)
  if (u<= acceptanceRatio) {
  beta[i, 2] = beta1.prime}
  else {beta[i, 2] = beta[i-1, 2]}
  
}

getPosterior=function(beta,x,y,currentBeta,mu0,sigma02){
  theta=cbind(1,x)%*%beta
  logLikelihood=sum(y*theta-exp(theta))
  logPrior=-((currentBeta-mu0)^2)/(2*sigma02)
  logPosterior=logLikelihood+logPrior
  return(logPosterior)
}

beta[nIter,]

mod=glm(y ~x, family='poisson')
mod$coe
