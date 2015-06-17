#1
data=read.csv('D:/study/235/235-roaches.csv')
x=data$count
hist(x)

N=length(x)

K=2
pi=rep(1/K, K)
mu=matrix(c(runif(1,min(x),max(x)),0),K,1)
d=matrix(NA,N,K)
p=matrix(NA,N,K)

# I run the algorithm for 20 iterations; in general, we should stop when the assignments don't change anymore.

for(i in 1:20){
  # Calculating the distance of each data point from the K centroids.
  for (k in 1:K){
    d[, k] <- pi[k]*dpois(x, mu[k])
  }

  total.d <- rowSums(d)
  for (k in 1:K){
    p[, k] <- pi[k]*dpois(x, mu[k])/total.d
  }


  # M-step
  pi <- colSums(p)/N

  sum.p <- colSums(p)
  for (k in 1:K){
    mu[k] <- sum(p[, k]*x)/sum.p[k]
  }
}

mu

# Hard clustering: assigning each data point to the component with the highest probability
clus <- apply(p, 1, which.max)
data <- data.frame(x, clus = factor(clus))

# Simulating data from the fitted mixture model for comparison.
Y <- NULL
for(k in 1:K){

  Y <-  rpois(round(pi[k]*50000), mu[k])

}

mix <- data.frame(y1 = Y[, 1], y2 = Y[, 2])


library(MASS)
library(RColorBrewer)
k <- 10
my.cols <- rev(brewer.pal(k, "RdYlBu"))

z <- kde2d(Y, n=100)

#2
data=read.csv('D:/study/235/235-clusterinData.csv')
N=dim(data)[1]


K = 3

# Initializing the parameters
pi <- rep(1/K, K)
mu <- cbind(runif(K, min(X[, 1]), max(X[, 1])), runif(K, min(X[, 2]), max(X[, 2])) )
sigma <- array(rbind(c(1, 0), c(0, 1)), c(2, 2, K))

d <- matrix(NA, N, K)
p <- matrix(NA, N, K)
# I run the algorithm for 20 iterations; in general, we should stop when the assignments don't change anymore.

for(i in 1:20){
  # Calculating the distance of each data point from the K centroids.
  for (k in 1:K){
    d[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])
  }

  total.d <- rowSums(d)
  for (k in 1:K){
    p[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])/total.d
  }


  # M-step
  pi <- colSums(p)/N

  sum.p <- colSums(p)
  for (k in 1:K){
    mu[k, ] <- colSums(p[, k]*X)/sum.p[k]
  }

  for (k in 1:K){
    temp <- 0
    for(j in 1:N){
      temp <- temp + p[j, k]*(X[j, ] - mu[k, ])%*%t(X[j, ] - mu[k, ])
    }
    sigma[, , k] <- temp/sum.p[k]
  }
}



K = 4

# Initializing the parameters
pi <- rep(1/K, K)
mu <- cbind(runif(K, min(X[, 1]), max(X[, 1])), runif(K, min(X[, 2]), max(X[, 2])) )
sigma <- array(rbind(c(1, 0), c(0, 1)), c(2, 2, K))

d <- matrix(NA, N, K)
p <- matrix(NA, N, K)
# I run the algorithm for 20 iterations; in general, we should stop when the assignments don't change anymore.

for(i in 1:20){
  # Calculating the distance of each data point from the K centroids.
  for (k in 1:K){
    d[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])
  }

  total.d <- rowSums(d)
  for (k in 1:K){
    p[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])/total.d
  }


  # M-step
  pi <- colSums(p)/N

  sum.p <- colSums(p)
  for (k in 1:K){
    mu[k, ] <- colSums(p[, k]*X)/sum.p[k]
  }

  for (k in 1:K){
    temp <- 0
    for(j in 1:N){
      temp <- temp + p[j, k]*(X[j, ] - mu[k, ])%*%t(X[j, ] - mu[k, ])
    }
    sigma[, , k] <- temp/sum.p[k]
  }



}


K = 5

# Initializing the parameters
pi <- rep(1/K, K)
mu <- cbind(runif(K, min(X[, 1]), max(X[, 1])), runif(K, min(X[, 2]), max(X[, 2])) )
sigma <- array(rbind(c(1, 0), c(0, 1)), c(2, 2, K))

d <- matrix(NA, N, K)
p <- matrix(NA, N, K)
# I run the algorithm for 20 iterations; in general, we should stop when the assignments don't change anymore.

for(i in 1:20){
  # Calculating the distance of each data point from the K centroids.
  for (k in 1:K){
    d[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])
  }

  total.d <- rowSums(d)
  for (k in 1:K){
    p[, k] <- pi[k]*dmvnorm(X, mu[k, ], sigma[, , k])/total.d
  }


  # M-step
  pi <- colSums(p)/N

  sum.p <- colSums(p)
  for (k in 1:K){
    mu[k, ] <- colSums(p[, k]*X)/sum.p[k]
  }

  for (k in 1:K){
    temp <- 0
    for(j in 1:N){
      temp <- temp + p[j, k]*(X[j, ] - mu[k, ])%*%t(X[j, ] - mu[k, ])
    }
    sigma[, , k] <- temp/sum.p[k]
  }



}



