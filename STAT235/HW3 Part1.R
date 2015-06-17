library(coda)
data=read.csv('D:/study/235/forestfires.csv')

levels(data$month) = factor(c('0','1','2','3','4','5','6','7','8','9','10','11'))
data$month = as.numeric(as.character(data$month))
levels(data$day) = factor(c('0','1','2','3','4','5','6'))
data$day = as.numeric(as.character(data$day))

data=as.matrix(data)

nIterations = 10000

x = data[, 1:12]
y = data[,13]

n = dim(x)[1]
p = dim(x)[2] 
d = p+1

# standardizing the predictors
x=scale(x)

# Adding a column of 1's for the intercept
x = cbind(rep(1, n), x)

nu0 = mean(y)
sigma02 = var(y)

# The parameters of normal prior N(mu0, tau02) for beta's.
mu0 = rep(0, d)
tau02 = rep(100, d)

regSampProper = function(y, x){
  # The following two matrices hold the posterior samples for sigma2 and beta.
  sigma2 = rep(1, nIterations);
  beta = matrix(rep(0, nIterations*(p+1)), nrow = (p+1), ncol=nIterations);
  
  # These are y*, x* as described in the course notes.
  yStar = c(y, mu0);
  xStar = rbind(x, diag(d));
  
  for(i in 2:nIterations){
    # A contains the diagonal elements of SigmaStar^(-1/2)
    A = 1/sqrt(c(sigma2[i-1]*rep(1, n), tau02));
    # newX is SigmaStar^(-1/2)*xStar
    newX = matrix(rep(A, d), nrow=(n+d), ncol = d)*xStar;
    # newY is SigmaStar^(-1/2)*yStar
    newY = A*yStar;
    
    # Now we can get mu_n and Lambda_n as described in the course note based
    # on newX and newY
    L2 = chol(t(newX)%*%newX);
    invL2 = solve(L2);
    Lambda_n = invL2%*%t(invL2);
    mu_n = Lambda_n%*%t(newX)%*%newY;
    
    # We can now sample from the posterior distribution given sigma^2. As
    # discussed in the class, we use the Cholesky decomposition for this.
    L3 = chol(Lambda_n);
    u = rnorm(d, 0, 1);
    beta[, i] = u%*%L3 + t(mu_n);
    
    # Now, given the current beta's, we sample from the posterior distribution of
    #sigma2, which is Inv-chi2(nu_n, sigma02_n). 
    eta = x%*%beta[, i];
    eps = y - eta;
    nu = sum(eps^2)/(n);
    nu_n = nu0 + n;
    sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
    
    # I sample a new value for sigma2 from Inv-chi2(nu_n, sigma02_n), to do this, I 	
    # sample from z ~ chi2(nu_n) and then my sample from Inv-chi2 would be 	
    # nu_n*sigma02_n/z.
    z = rchisq(1, nu_n);
    sigma2[i] = nu_n*sigma02_n/z;
    
  }
  return(list(beta = beta, sigma2 = sigma2))
}



regSampImproper = function(y, x){
  # The following two matrices hold the posterior samples for sigma2 and beta.
  sigma2 = rep(1, nIterations);
  beta = matrix(rep(0, nIterations*(p+1)), nrow = (p+1), ncol=nIterations);
  
  #sigma=1
  
  for(i in 2:nIterations){
    
    # A contains the diagonal elements of SigmaStar^(-1/2)
    A = 1/sqrt(c(sigma2[i-1]*rep(1, n), tau02));
    
    # Now we can get mu_n and Lambda_n as described in the course note based
    # on newX and newY
    L2 = chol(t(x)%*%x);
    invL2 = solve(L2);
    Lambda_n = invL2%*%t(invL2);
    mu_n = Lambda_n%*%t(x)%*%y;
    
    # We can now sample from the posterior distribution given sigma^2. As
    # discussed in the class, we use the Cholesky decomposition for this.
    L3 = chol(Lambda_n);
    u = rnorm(d, 0, 1);
    beta[, i] = u%*%L3 + t(mu_n);
    
    # Now, given the current beta's, we sample from the posterior distribution of
    #sigma2, which is Inv-chi2(nu_n, sigma02_n). 
    eta = x%*%beta[, i];
    eps = y - eta;
    nu = sum(eps^2)/(n-p-1);
    nu_n = n-p-1;
    sigma02_n = sum(eps^2)/nu_n;
    
    # I sample a new value for sigma2 from Inv-chi2(nu_n, sigma02_n), to do this, I 	
    # sample from z ~ chi2(nu_n) and then my sample from Inv-chi2 would be 	
    # nu_n*sigma02_n/z.
    z = rchisq(1, nu_n);
    sigma2[i] = nu_n*sigma02_n/z;
  }
  return(list(beta = beta, sigma2 = sigma2))
}

# Running the program

### a) Proper priors
regRes=regSampProper(y,x)


beta=mcmc(t(regRes$beta))
beta.proper=colMeans(beta)
plot(beta)

sigma2=mcmc(regRes$sigma2)
plot(sigma2) 


### b) Improper priors

regRes=regSampImproper(y, x)

beta=mcmc(t(regRes$beta))
beta.improper=colMeans(beta)
plot(beta)

sigma2=mcmc(regRes$sigma2)
plot(sigma2)


