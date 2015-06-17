#2
library(mvtnorm)
library(mclust)

data=read.csv('D:/study/235/235-clusterinData.csv')
X=as.matrix(data)

#a
fit1=Mclust(X,G=3) #K=3
fit2=Mclust(X,G=4) #K=4
fit3=Mclust(X,G=5) #K=5

AIC=c(AIC(fit1),AIC(fit2),AIC(fit3)) #c(1961.719, 1730.239, 1736.260)
BIC=c(BIC(fit1),BIC(fit2),BIC(fit3)) #c(2023.510, 1813.839, 1841.669)
#Both AIC and BIC are smallest among the three model when K=4

summary(fit2)
para=fit2$para
mu=para$mean
sigma=para$variance$sigma

#X is a mixture of 40/280 multi-N(mu[,1],sigma[,,1])

#b
K=4

# Initializing the parameters
pi <- rep(1/K, K)
mu <- cbind(runif(K, min(X[, 1]), max(X[, 1])), runif(K, min(X[, 2]), max(X[, 2])) )
sigma <- array(rbind(c(1, 0), c(0, 1)), c(2, 2, K))

d <- matrix(NA, N, K)
p <- matrix(NA, N, K)

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
# Hard clustering: assigning each data point to the component with the highest probability
clus <- apply(p, 1, which.max)
data<- data.frame(x1 = X[, 1], x2 = X[, 2], clus = factor(clus))

#c
fit=kmeans(X, 4)
mydata=data.frame(X, fit$cluster)
d =dist(mydata, method = "euclidean")
fit=hclust(d, method="ward")
# plot
plot(fit,main=)
rect.hclust(fit, k=4, border="red")
