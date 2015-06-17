# Problem 1
n = 100

# when alpha = 0.3
alpha1 = .3
sigmax1 = 1-alpha1^2
X1 = rnorm(n, 0, sqrt(sigmax1))

# AR random process Y
Y1 = rep(NA, n)
Y1[1] = X1[1]
for (i in 2:n){
  Y1[i] = alpha1*Y1[i-1]+X1[i]
}

# Autocorrelation Ry
Ry1 = acf(Y1,lag = 'n-1')$acf
Ry1minus = rep(NA,n-1)
for (i in 1:n-1){
  Ry1minus[i] = Ry1[n+1-i]
}
Ry1.new = c(Ry1minus, Ry1)

# power density Sy(f)
Sy1 = periodogram(Ry1)$spec
Freq1 = periodogram(Ry1)$freq
n1 = length(Sy1)
Sy1minus = rep(NA,n1)
for (i in 1:n1){
  Sy1minus[i] = Sy1[n1+1-i]
}
Sy1 = c(Sy1minus,Sy1)
Freq1minus = rep(NA,n1)
for (i in 1:n1){
  Freq1minus[i] = -Freq1[n1+1-i]
}
Freq1 = c(Freq1minus, Freq1)

# when alpha = 0.95
alpha2 = .95
sigmax2 = 1-alpha2^2
X2 = rnorm(n, 0, sqrt(sigmax2))
# AR random process Y
Y2 = rep(NA, n)
Y2[1] = X2[1]
for (i in 2:n){
  Y2[i] = alpha2*Y2[i-1]+X2[i]
}

# Autocorrelation Ry
Ry2 = acf(Y2,lag = 'n-1')$acf
Ry2minus = rep(NA,n-1)
for (i in 1:n-1){
  Ry2minus[i] = Ry2[n+1-i]
}
Ry2.new = c(Ry2minus, Ry2)

# power density Sy(f)
Sy2 = periodogram(Ry2)$spec
Freq2 = periodogram(Ry2)$freq
n2 = length(Sy2)
Sy2minus = rep(NA,n2)
for (i in 1:n2){
  Sy2minus[i] = Sy2[n2+1-i]
}
Freq2minus = rep(NA,n2)
for (i in 1:n2){
  Freq2minus[i] = -Freq2[n2+1-i]
}
Sy2 = c(Sy2minus,Sy2)
Freq2 = c(Freq2minus, Freq2)

# plots
par(mfrow = c(1,2))
plot(Y1 ~ c(1:n), ylim = c(-4,4), xlab = 'n', ylab = 'Yn',
     main = 'AR randomprocess for alpha = 0.3')
abline(h = 0,col = 'blue')
plot(Y2 ~ c(1:n), ylim = c(-4,4), xlab = 'n', ylab = 'Yn',
     main = 'AR randomprocess for alpha = 0.95')
abline(h = 0,col = 'blue')


plot(Ry1.new ~ c((-n+1):(n-1)), xlab = 'n', ylab = 'Ry',
     main = 'Autocorrelation for alpha = 0.3',type = 'h')
plot(Ry2.new ~ c((-n+1):(n-1)), xlab = 'n', ylab = 'Ry',
     main = 'Autocorrelation for alpha = 0.95',type = 'h')

plot(Sy1 ~ Freq1,type='l', xlab = 'frequency', ylab = 'power spectral density',
     main = 'Power Spectral Density for alpha = 0.3')
plot(Sy2 ~ Freq2,type='l',  xlab = 'frequency', ylab = 'power spectral density',
     main = 'Power Spectral Density for alpha = 0.95')





# Problem 2
N = 1000
Xn = matrix(rep(NA, 100*N),100,N)
for (i in 1:100){
  for (n in 1:N)
  Xn[i,n] = cos(0.2*pi*n+runif(1,-pi,pi))}

plot(Xn[1,]~c(1:N), ylim = c(-1,1), xlab = 'n', ylab = 'Xn',pch='.')
for (i in 1:100){
  points(Xn[i,]~c(1:N),pch='.')
}


### Problem 3
n = 150
r = .5
sigmax = sqrt(1-r^4)
X = rnorm(n,0,sigmax)
h2= -r^2
Y = rep(0, n)
Y_hat = rep(0,n)
for (i in 3:n){
  Y[i] = -r^2*Y[i-2]+X[i]
  Y_hat[i] = h2*Y[i-2]
}
par(mfrow = c(1,2))
plot(Y[101:n] ~ c(1:50),type = 'h', lty = 1, xlab = 'n', ylab = 'Y',
     main = '(i) r = 0.5')
points(Y[101:n] ~ c(1:50),col='red', pch = 0)
points(Y_hat[101:n] ~ c(1:50),type = 'h',lty = 3)
points(Y_hat[101:n] ~ c(1:50), col='blue', pch = 16)
abline(h = 0)
legend('topright',c('Yn','Predicted Yn'), lty = c(1,3), col = c('red','blue'), pch = c(0,16),)

n = 150
r = .95
sigmax = sqrt(1-r^4)
X = rnorm(n,0,sigmax)
h2= -r^2
Y = rep(0, n)
Y_hat = rep(0,n)
for (i in 3:n){
  Y[i] = -r^2*Y[i-2]+X[i]
  Y_hat[i] = h2*Y[i-2]
}
plot(Y[101:n] ~ c(1:50),type = 'h', lty = 1, xlab = 'n', ylab = 'Y',
     main = '(ii) r = 0.95')
points(Y[101:n] ~ c(1:50),col='red', pch = 0)
points(Y_hat[101:n] ~ c(1:50),type = 'h',lty = 3)
points(Y_hat[101:n] ~ c(1:50), col='blue', pch = 16)
abline(h = 0)
legend('topright',c('Yn','Predicted Yn'), lty = c(1,3), col = c('red','blue'), pch = c(0,16),)

