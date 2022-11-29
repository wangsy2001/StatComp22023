## -----------------------------------------------------------------------------
library(StatComp22023)

## ---- warning=FALSE-----------------------------------------------------------
library(knitr)
data(data)

## -----------------------------------------------------------------------------
head(WorldPhones, 4)

## -----------------------------------------------------------------------------
kable(head(WorldPhones, 4))

## -----------------------------------------------------------------------------
parity<-function(x){if(x%%2==1){print('odd')} else{print('even')}}
parity(5)

## ---- tidy=TRUE---------------------------------------------------------------
parity<-function(x){if(x%%2==1){print('odd')} else{print('even')}}
parity(5)

## ---- warning=FALSE-----------------------------------------------------------
library(ggplot2)

## -----------------------------------------------------------------------------
ggplot(data=beaver1, aes(x=day, y=temp, group=day)) + 
  geom_boxplot() + 
  scale_x_continuous(breaks = c(346, 347))

## ---- fig.width=3, fig.height=4-----------------------------------------------
ggplot(data=beaver1, aes(x=day, y=temp, group=day)) + 
  geom_boxplot() + 
  scale_x_continuous(breaks = c(346, 347))

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
set.seed(1)
N <- 5000
U <- runif(N)

## -----------------------------------------------------------------------------
F_inverse = function(u, a, b){
  b*(1-u)^(-1/a)
}
x <- F_inverse(U, 2, 2)

## -----------------------------------------------------------------------------
pareto <- function(x, a, b){
  a*(b/x)^(a-1)*b/x^2
}

## -----------------------------------------------------------------------------
hist(x, freq=F, breaks=100,xlim=c(0,20))
#draw density function
lines(seq(2,20,0.05), pareto(seq(2,20,0.05), 2, 2), lwd=2)

## -----------------------------------------------------------------------------
n <- 10000 #number of samples
x <- numeric(0) #samples
f = function(x) 1/beta(2,3)*x*(1-x)^2
g = function(x) I(x>0)*I(x<1) 
c <- 2

## -----------------------------------------------------------------------------
set.seed(1)
i=k=0
while(i < n){
  y <- runif(1) #y ~ g
  u <- runif(1) #u ~ u(0,1)
  if(u <= f(y)/c/g(y)){ # 接收y的条件
    x[k] <- y
    k <- k+1
  }
  i <- i+1 
}

## -----------------------------------------------------------------------------
hist(x, freq=F, ylim=c(0,2))
lines(seq(0,1,0.01), dbeta(seq(0,1,0.01),2,3), lwd=2)

## -----------------------------------------------------------------------------
n <- 1000
lambda <- rgamma(n,4,2)

## -----------------------------------------------------------------------------
x <- rexp(n, lambda)

## -----------------------------------------------------------------------------
hist(x, breaks=30, xlim=c(0,6))

## -----------------------------------------------------------------------------
pareto <- function(y, r, beta){
  r*(beta/(y+beta))^(r-1)*beta/(y+beta)^2
}

## -----------------------------------------------------------------------------
hist(x, breaks=30, freq=F, xlim=c(0,6))
lines(seq(0.05,6,0.01), pareto(seq(0.05,6,0.01), 4, 2), lwd=2)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}#??????
}

## -----------------------------------------------------------------------------
time_all = c()
num <- c(1e4, 2e4, 4e4, 6e4, 8e4)
for(k in 1:5){
  n <- num[k]
  time_each <- 0
  for (i in 1:100){
    test<-sample(1:n)
    time_each <- time_each + system.time(quick_sort(test))[1]
  }
  time_all[k] = time_each
}

## -----------------------------------------------------------------------------
x = num*log(num)
relation <- lm(time_all ~ x)
plot(x, time_all, main = "Regression",pch = 16,xlab = "n log(n)",ylab = "time")
abline(relation,cex = 1.3)

## -----------------------------------------------------------------------------
1-(exp(2)-1 - 4*(exp(1)-1)^2 + 2*exp(1))/4/((exp(2)-1)/2-(exp(1)-1)^2)

## -----------------------------------------------------------------------------
set.seed(1)
m <- 5000
u <- runif(m)

## -----------------------------------------------------------------------------
sample_MC = exp(u)
cat('the MC simulation and the variance is', mean(sample_MC), var(sample_MC))

## -----------------------------------------------------------------------------
sample_anti = (exp(u)+exp(1-u))/2
cat('the antithetic simulation and the variance is', mean(sample_anti), var(sample_anti))

## -----------------------------------------------------------------------------
cat('the percent reduction in variance is ', 1-var(sample_anti)/var(sample_MC))

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
g <- function(x){
  x^2/sqrt(2*pi)*exp(-x^2/2)
}

## -----------------------------------------------------------------------------
f1 <- function(x)  dnorm(x, mean=1.5, sd=1.5) + dnorm(x, mean=0.5, sd=1.5)

## -----------------------------------------------------------------------------
f2 <- function(x)  dnorm(x, mean=4, sd=1.5) + dnorm(x, mean=-2, sd=1.5)

## -----------------------------------------------------------------------------
set.seed(1)
x1 <- abs(rnorm(1e5, 1.5, 1.5)-1)+1 #f1
x2 <- abs(rnorm(1e5, 4, 1.5)-1)+1 #f2

## -----------------------------------------------------------------------------
x <- seq(1, 5, 0.05)
plot(x, g(x), type='l', ylim=c(0,0.6))
lines(x, f1(x), col='red')
lines(x, f2(x), col='blue')
legend("topright",legend=c('g','f1','f2'),col=c('black','red','blue'),lty=c(1,1,1))

## -----------------------------------------------------------------------------
integrate(f1, 1, Inf)
integrate(f2, 1, Inf)

## -----------------------------------------------------------------------------
f1_eva = mean(g(x1)/f1(x1)) #f1
f2_eva = mean(g(x2)/f2(x2)) #f2
cat(f1_eva, f2_eva)

## -----------------------------------------------------------------------------
f1_eva = var(g(x1)/f1(x1)) #f1
f2_eva = var(g(x2)/f2(x2)) #f2
cat(f1_eva, f2_eva)

## -----------------------------------------------------------------------------
Fx <- c(0,0.2, 0.4, 0.6, 0.8,1)
I <- -log(1-Fx*(1-exp(-1)))
I

## -----------------------------------------------------------------------------
f <- function(x){
 exp(-x)/(1-exp(-1))
}

F_inv <- function(u, Ij){
  -log(exp(-Ij) - u/5*(1-exp(-1)))
}

g <- function(x) {
  exp(-x)/(1+x^2)
}

## -----------------------------------------------------------------------------
set.seed(1)
result = numeric(5)
se = numeric(5)
n <- 2000
for(j in 1:5){
  u <- runif(n,0,1) 
  x <- F_inv(u, I[j])
  fg <- g(x) / f(x)
  result[j] <- mean(fg)
  se[j] <- sd(fg)
}
mean(result)
mean(se)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
print(c(-2/sqrt(30)*qnorm(0.975), 2/sqrt(30)*qnorm(0.975)))

## -----------------------------------------------------------------------------
n <- 30
m <- 1e5
mu <- numeric(100)
set.seed(1)
for(i in 1:m){
  x <- rlnorm(n,0,2)
  y <- log(x)
  mu[i] = mean(y)
}
quantile(mu, probs=c(0.025, 0.975))

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
set.seed(1)
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

# generate samples under H1 to estimate power
sigma1 <- 1
sigma2 <- 1.5
m <- 1000
n <- 200
power <- mean(replicate(m, expr={
                        x <- rnorm(n, 0, sigma1)
                        y <- rnorm(n, 0, sigma2)
                        count5test(x, y)
                                 }
                        )
             )
print(power)

## -----------------------------------------------------------------------------
set.seed(1)
p.value <- replicate(m, expr={
                        x <- rnorm(n, 0, sigma1)
                        y <- rnorm(n, 0, sigma2)
                        var.test(x, y, alternative = "two.sided")$p.value
                            }
                   )
         
print(sum(p.value < 0.055)/m)

## -----------------------------------------------------------------------------
set.seed(1)
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

# generate samples under H1 to estimate power
sigma1 <- 1
sigma2 <- 1.5
m <- 1000
n <- 50
power <- mean(replicate(m, expr={
                        x <- rnorm(n, 0, sigma1)
                        y <- rnorm(n, 0, sigma2)
                        count5test(x, y)
                                 }
                        )
             )
print(power)

## -----------------------------------------------------------------------------
set.seed(1)
p.value <- replicate(m, expr={
                        x <- rnorm(n, 0, sigma1)
                        y <- rnorm(n, 0, sigma2)
                        var.test(x, y, alternative = "two.sided")$p.value
                            }
                   )
         
print(sum(p.value < 0.055)/m)

## -----------------------------------------------------------------------------
set.seed(1)
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

# generate samples under H1 to estimate power
sigma1 <- 1
sigma2 <- 1.5
m <- 1000
n <- 20
power <- mean(replicate(m, expr={
                        x <- rnorm(n, 0, sigma1)
                        y <- rnorm(n, 0, sigma2)
                        count5test(x, y)
                                 }
                        )
             )
print(power)

## -----------------------------------------------------------------------------
set.seed(1)
p.value <- replicate(m, expr={
                        x <- rnorm(n, 0, sigma1)
                        y <- rnorm(n, 0, sigma2)
                        var.test(x, y, alternative = "two.sided")$p.value
                            }
                   )
         
print(sum(p.value < 0.055)/m)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
x1 = c(rep(1, 651), rep(0, 349))
x2 = c(rep(1, 676), rep(0, 324))
t.test(x1, x2)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
library(boot)
n = dim(aircondit)[1]
lambda = n/sum(aircondit$hours)
lambda # the MLE of the hazard rate

## -----------------------------------------------------------------------------
set.seed(2022)
data <- aircondit$hours
B <- 1e4
lambdastar <- numeric(B)
for(b in 1:B){
  datastar <- sample(data,replace=TRUE)
  lambdastar[b] <- n/sum(datastar)
}
round(c(bias=mean(lambdastar)-lambda, 
        se.boot=sd(lambdastar)), 4)

## -----------------------------------------------------------------------------
se.boot=sd(lambdastar)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
boot.mean <- function(x,i){
  mean(x[i])
  }

redata <- boot(data=aircondit$hours, statistic=boot.mean, R = 1e3)
result <- boot.ci(redata, type=c("norm","basic","perc","bca"))
cbind(norm=result$norm[2:3], basic=result$basic[4:5],
      percentile=result$percent[4:5],bca=result$bca[4:5])

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
n <- 20
M <- 1e4
set.seed(1)
boot.mean <- function(x,i) mean(x[i])
result.norm <- result.basic <- result.perc <- matrix(NA ,M, 2)

for(i in 1:M){
  data = rnorm(n, 0, 1)
  de <- boot(data=data,statistic=boot.mean, R = 1e3)
  result <- boot.ci(de,type=c("norm","perc","basic"))
  result.norm[i,]<-result$norm[2:3]
  result.perc[i,]<-result$percent[4:5]
  result.basic[i,]<-result$basic[4:5]
}
cat('norm =',mean(result.norm[,1]<=0 & result.norm[,2]>=0),
    'perc =',mean(result.perc[,1]<=0 & result.perc[,2]>=0),
    'basic =',mean(result.basic[,1]<=0 & result.basic[,2]>=0))

## -----------------------------------------------------------------------------
cat('miss on the left: norm =',mean(result.norm[,1]>0), 'perc =',mean(result.perc[,1]>0), 'basic =',mean(result.basic[,1]>0))
cat('\nmiss on the right: norm =',mean(result.norm[,2]<0), 'perc =',mean(result.perc[,2]<0), 'basic =',mean(result.basic[,2]<0))

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
data = bootstrap::scor
n=88

## -----------------------------------------------------------------------------
FPC <- function(x, i){
  x = x[i,]
  Sigma <- cor(x)
  eigen(Sigma)$values[1]/sum(eigen(Sigma)$values)
}
theta.hat <- FPC(data, 1:n)

## -----------------------------------------------------------------------------
theta.jack <- numeric(n)
for(i in 1:n){
  theta.jack[i] <- FPC(data,(1:n)[-i])
}
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
round(c(est=theta.hat,bias=bias.jack,
        se=se.jack),4)

## -----------------------------------------------------------------------------
rm(list=ls()) 

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)

## -----------------------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
  J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  J4 <- lm(log(y) ~ log(x))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
  yhat4 <- exp(logyhat4)
  e4[k] <- magnetic[k] - yhat4
}

## -----------------------------------------------------------------------------
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
cross <- combn(n,2)
e1 <- e2 <- e3 <- e4 <- numeric(dim(cross)[2])
for (group in 1:dim(cross)[2]) {
  k <- cross[,group]
  y <- magnetic[-k]
  x <- chemical[-k]
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[group] <- mean((magnetic[k] - yhat1)^2)
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
  J2$coef[3] * chemical[k]^2
  e2[group] <- mean((magnetic[k] - yhat2)^2)
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[group] <- mean((magnetic[k] - yhat3)^2)
  J4 <- lm(log(y) ~ log(x))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
  yhat4 <- exp(logyhat4)
  e4[group] <- mean((magnetic[k] - yhat4)^2)
}

## -----------------------------------------------------------------------------
c(mean(e1), mean(e2), mean(e3), mean(e4))

## -----------------------------------------------------------------------------
lm(magnetic ~ chemical + I(chemical^2))

## -----------------------------------------------------------------------------
rm(list=ls()) 

## -----------------------------------------------------------------------------
data(data) 
iris <- iris[1:50,1:4] #the species setosa

## -----------------------------------------------------------------------------
boot.rho <- function(data, i){
  x <- as.matrix(data[ , 1:2])
  y <- as.matrix(data[i, 3:4])
  cor.test(x, y, method = "spearman", exact = FALSE)$estimate
}
rho.hat <- boot.rho(iris, 1:50)
rho.hat

## -----------------------------------------------------------------------------
library(boot)
set.seed(1)
boot.obj <- boot(data = iris, statistic = boot.rho, R = 9999, sim = "permutation")
p.perm <- mean(c(boot.obj$t, boot.obj$t0) > rho.hat)
p.perm

## -----------------------------------------------------------------------------
test <- cor.test(as.matrix(iris[,1:2]), as.matrix(iris[,3:4]), method = "spearman", exact = FALSE)
test$p.value

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
dl <- function(x){
  exp(-abs(x))/2
}

lower <- uniroot(function(a) integrate(dl, -Inf, a)$value-0.025, lower=-6, upper=6)$root
upper <- uniroot(function(a) integrate(dl, -Inf, a)$value-0.975, lower=-6, upper=6)$root
c(lower, upper)

## -----------------------------------------------------------------------------
set.seed(1)
Gelman.Rubin <- function(chain) {
  chain <- as.matrix(chain)
  l <- ncol(chain)
  
  B <- l * var(rowMeans(chain))       
  chain.w <- apply(chain, 1, "var") 
  W <- mean(chain.w)              
  v.hat <- W*l/(l-1) + (B/l)   
  r.hat <- v.hat / W           
  return(r.hat)
}

## ---- eval=FALSE--------------------------------------------------------------
#  Metropolis <- function(sigma, x0=c(-9,-3,3,9), N_max=5000) { #generate 4 chains to calculate R
#      x <- matrix(0, nrow=4, ncol=N_max)
#      x[,1] <- x0
#      acc = 0
#      Rhat <- NULL
#      set.seed(1)
#      for (i in 2:N_max) {
#          y <- c(rnorm(1, x[1,i-1], sigma),rnorm(1, x[2,i-1], sigma),rnorm(1, x[3,i-1], sigma),rnorm(1, x[4,i-1], sigma))
#          u <- runif(1)
#          for(k in 1:4){
#            if (u <= (dl(y[k]) / dl(x[k,i-1]))){
#              x[k,i] <- y[k]
#              acc <- acc + 1
#            }
#            else x[k,i] <- x[k,i-1]
#          }
#          # calculate Rhat
#          psi <- t(apply(x[,1:i], 1, cumsum))#compute diagnostic statistics
#          for (k in 1:nrow(psi))
#              psi[k,] <- psi[k,] / (1:ncol(psi))
#          Rhat <- c(Rhat, Gelman.Rubin(psi))
#          if((i >= 3000) & (Rhat[i-1]<1.2)) break
#  
#  
#      }
#      return(list(x=x[,1:i], acc=acc/4/i, Rhat=Rhat))
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  sigma <- c(0.05, 0.5, 2, 9)
#  chain1 <- Metropolis(sigma[1])
#  chain2 <- Metropolis(sigma[2])
#  chain3 <- Metropolis(sigma[3])
#  chain4 <- Metropolis(sigma[4])

## ---- eval=FALSE--------------------------------------------------------------
#  #number of candidate points rejected
#  acc.r <- data.frame(sigma=sigma,
#                      acceptance.rates=c(chain1$acc, chain2$acc, chain3$acc, chain4$acc),
#                      convergent.length=c(sum(chain1$R>1.2)+1, sum(chain2$R>1.2)+1, sum(chain3$R>1.2)+1, sum(chain4$R>1.2)+1))
#  
#  acc.r

## ---- eval=FALSE--------------------------------------------------------------
#  par(mfrow = c(2, 2))
#  plot(chain1$x[1,], type="l", xlab='sigma=0.05')
#  abline(h=c(lower, upper))
#  plot(chain2$x[1,], type="l", xlab='sigma=0.5')
#  abline(h=c(lower, upper))
#  plot(chain3$x[1,], type="l", xlab='sigma=2')
#  abline(h=c(lower, upper))
#  plot(chain4$x[1,], type="l", xlab='sigma=9')
#  abline(h=c(lower, upper))

## ---- eval=FALSE--------------------------------------------------------------
#  par(mfrow = c(2, 2))
#  plot(chain1$Rhat[-c(1:500)], type="l", xlab=bquote(sigma == .(round(sigma[1],3))), ylab="R")
#  abline(h=1.2, lty=2)
#  plot(chain2$Rhat[-c(1:500)], type="l", xlab=bquote(sigma == .(round(sigma[2],3))), ylab="R")
#  abline(h=1.2, lty=2)
#  plot(chain3$Rhat[-c(1:300)], type="l", xlab=bquote(sigma == .(round(sigma[3],3))), ylab="R")
#  abline(h=1.2, lty=2)
#  plot(chain4$Rhat[-c(1:300)], type="l", xlab=bquote(sigma == .(round(sigma[4],3))), ylab="R")
#  abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
N_max <- 15000 #max length of chain
k <- 4 #number of chains
burn <- 1000 #burn-in length
X1 <- matrix(0, nrow=4, ncol=N_max) #the 1st dimension of the bivariate normal chain
X2 <- matrix(0, nrow=4, ncol=N_max) #the 2rd dimension of the bivariate normal chain
rho <- 0.9 #correlation
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2
Rhat <- NULL


Gelman.Rubin <- function(chain) {
  chain <- as.matrix(chain)
  l <- ncol(chain)
  
  B <- l * var(rowMeans(chain))       
  chain.w <- apply(chain, 1, "var") 
  W <- mean(chain.w)              
  v.hat <- W*l/(l-1) + (B/l)   
  r.hat <- v.hat / W           
  return(r.hat)
}


#generate the chain
set.seed(1)
X1[,1] <- c(-5,-2,2,5) #initialize
X2[,1] <- c(-5,-2,2,5) #initialize
for (i in 2:N_max) { #calculate R with Gelman-Rubin method
  for (k in 1:4){
    x2 <- X2[k,i-1]
    m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X1[k,i] <- rnorm(1, m1, s1)
    x1 <- X1[k,i]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X2[k,i] <- rnorm(1, m2, s2)
  }
  
  # calculate Rhat
  # compute diagnostic statistics
  if(i>burn){
    psi1 <- t(apply(X1[,1:i], 1, cumsum))
    psi2 <- t(apply(X2[,1:i], 1, cumsum))
    psi <- (psi1+psi2)/2
    for (k in 1:nrow(psi))
        psi[k,] <- psi[k,] / (1:ncol(psi))
    Rhat <- c(Rhat, Gelman.Rubin(psi))
    if((i >= 3000) & (Rhat[i-burn]<1.2)) break
  }
}

X1 <- X1[1, (burn+1):i]
X2 <- X2[1, (burn+1):i]

## -----------------------------------------------------------------------------
sum(Rhat>1.2)

## -----------------------------------------------------------------------------
plot(Rhat, type="l", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
fm <- lm(X2 ~ X1)
coef <- summary(fm)$coef
coef[,1]

## -----------------------------------------------------------------------------
res <- X2 - coef[1,1] - coef[2,1]*X1
plot(X1, res, main='residual')

## -----------------------------------------------------------------------------
qqnorm(res)
ks.test(res,"pnorm",mean=mean(res),sd=sqrt(var(res)))

## -----------------------------------------------------------------------------
require(graphics)
group <- X1
group[group < -2] = 1
group[group < -1] = 2
group[group < 0] = 3
group[group < 1] = 4
group[group < 2] = 5
group[X1 > 2] = 6
bartlett.test(res ~ group)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
generate_sample <- function(alpha, beta, gamma, n=500, aM=1, aY=1){
  eM = rnorm(n, 0, 1)
  eY = rnorm(n, 0, 1)
  X <- rnorm(n, 2, 2)
  M <- aM + alpha*X + eM
  Y <- aY + beta*M + gamma*X + eY
  return(data.frame(X=X, M=M, Y=Y))
}
set.seed(1)
sam1 = generate_sample(0,0,1)
sam2 = generate_sample(0,1,1)
sam3 = generate_sample(1,0,1)

## -----------------------------------------------------------------------------
library(boot)
boot.perm <- function(data, i, condition){
  x <- data$X
  y <- data$Y
  m <- data$M
  if(condition == 'X') x <- x[i] #permuting
  if(condition == 'Y') y <- y[i] #permuting
  if(condition == 'M') m <- m[i] #permuting
  fy <- lm(y ~ m + x)
  fm <- lm(m ~ x)
  b = summary(fy)$coefficients[2,1]
  b.se = summary(fy)$coefficients[2,2]
  a = summary(fm)$coefficients[2,1]
  a.se = summary(fm)$coefficients[2,2]
  ab.se = sqrt(b^2*a.se^2+a^2*b.se^2)
  t = a*b/ab.se
}

## -----------------------------------------------------------------------------
set.seed(1)
boot.med <- boot(data = sam1, statistic = boot.perm, R = 1999, sim = "permutation", condition='X')
p.value.a <- mean(abs(c(boot.med$t0, boot.med$t)) > abs(boot.med$t0))

## -----------------------------------------------------------------------------
set.seed(1)
boot.med <- boot(data = sam1, statistic = boot.perm, R = 1999, sim = "permutation", condition='Y')
p.value.b <- mean(abs(c(boot.med$t0, boot.med$t)) > abs(boot.med$t0))

## -----------------------------------------------------------------------------
set.seed(1)
boot.med <- boot(data = sam1, statistic = boot.perm, R = 1999, sim = "permutation", condition='M')
p.value.m <- mean(abs(c(boot.med$t0, boot.med$t)) > abs(boot.med$t0))

## -----------------------------------------------------------------------------
c('a=0'=p.value.a, 'b=0'=p.value.b, 'ab=0'=p.value.m)

## -----------------------------------------------------------------------------
set.seed(1)
boot.med <- boot(data = sam2, statistic = boot.perm, R = 1999, sim = "permutation", condition='X')
p.value.a <- mean(abs(c(boot.med$t0, boot.med$t)) > abs(boot.med$t0))

## -----------------------------------------------------------------------------
set.seed(1)
boot.med <- boot(data = sam2, statistic = boot.perm, R = 1999, sim = "permutation", condition='Y')
p.value.b <- mean(abs(c(boot.med$t0, boot.med$t)) > abs(boot.med$t0))

## -----------------------------------------------------------------------------
set.seed(1)
boot.med <- boot(data = sam2, statistic = boot.perm, R = 1999, sim = "permutation", condition='M')
p.value.m <- mean(abs(c(boot.med$t0, boot.med$t)) > abs(boot.med$t0))

## -----------------------------------------------------------------------------
c('a=0'=p.value.a, 'b=0'=p.value.b, 'ab=0'=p.value.m)

## -----------------------------------------------------------------------------
set.seed(1)
boot.med <- boot(data = sam3, statistic = boot.perm, R = 1999, sim = "permutation", condition='X')
p.value.a <- mean(abs(c(boot.med$t0, boot.med$t)) > abs(boot.med$t0))

## -----------------------------------------------------------------------------
set.seed(1)
boot.med <- boot(data = sam3, statistic = boot.perm, R = 1999, sim = "permutation", condition='Y')
p.value.b <- mean(abs(c(boot.med$t0, boot.med$t)) > abs(boot.med$t0))

## -----------------------------------------------------------------------------
set.seed(1)
boot.med <- boot(data = sam3, statistic = boot.perm, R = 1999, sim = "permutation", condition='M')
p.value.m <- mean(abs(c(boot.med$t0, boot.med$t)) > abs(boot.med$t0))

## -----------------------------------------------------------------------------
c('a=0'=p.value.a, 'b=0'=p.value.b, 'ab=0'=p.value.m)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
g <- function(alpha, b1, b2, b3, f0, x1, x2, x3){
  tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
  p <- 1/(1+tmp)
  mean(p) - f0
}

alpha_cal <- function(N, b1, b2, b3, f0){
  x1 <- rpois(N, 1)
  x2 <- rexp(N, rate = 1)
  x3 <- rbinom(N, size=1, prob=0.5)
  solution <- uniroot(function(a) g(a, b1, b2, b3, f0, x1, x2, x3), c(-30,30))
  solution$root
}

## -----------------------------------------------------------------------------
alpha = numeric(0)
for(f in c(0.1, 0.01, 0.001, 0.0001)){
  alpha = c(alpha, alpha_cal(1e6, 0, 1, -1, f))
}

## -----------------------------------------------------------------------------
plot(-log(c(0.1, 0.01, 0.001, 0.0001)), alpha, xlab = '-log f0')

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
u = c(11,8,27,13,16,0,23,10,24,2)
v = c(12,9,28,14,17,1,24,11,25,3)

## -----------------------------------------------------------------------------
log_like_dao <- function(lambda){
  sum((-u*exp(-lambda*u)+v*exp(-lambda*v))/(exp(-lambda*u)-exp(-lambda*v)))
}
uniroot(log_like_dao, c(0.01,1))$root

## -----------------------------------------------------------------------------
epo = 500
lambda0=1
n = length(u)
for(i in 1:epo){
  A = sum(((v+1/lambda0)*exp(-lambda0*v) - (u+1/lambda0)*exp(-lambda0*u)) / (exp(-lambda0*u)-exp(-lambda0*v)))
  lambda1 = -n/A
  if(abs(lambda0-lambda1)<1e-3){
    break
    }
  lambda0 = lambda1
}
i
lambda1

## -----------------------------------------------------------------------------
unlist(list(1,2,3))
is.vector(list(1,2,3))
as.vector(list(1,2,3))

## -----------------------------------------------------------------------------
x <- data.frame()
x

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
data(data)
data = iris[1:4]
data_every_column = apply(data, 2, scale01)
head(data_every_column)

## -----------------------------------------------------------------------------
data = iris
data_numeric_column = data.frame(lapply(data, function(x) if (is.numeric(x)) scale01(x) else x))
head(data_numeric_column)

## -----------------------------------------------------------------------------
rm(list = ls())

## -----------------------------------------------------------------------------
data(data)
data = iris[1:4]
vapply(data, sd, 0)

## -----------------------------------------------------------------------------
data = iris
vapply(data[vapply(data, is.numeric, logical(1))], sd, 0)

## -----------------------------------------------------------------------------
rm(list = ls())

## ---- echo=FALSE--------------------------------------------------------------
library(StatComp)

## -----------------------------------------------------------------------------
gibbsR <- function(N, thin) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- y <- 0
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rnorm(1, 0.9*y, sqrt(1-0.9^2))
      y <- rnorm(1, 0.9*x, sqrt(1-0.9^2))
    }
    mat[i, ] <- c(x, y)
  }
  mat
}

## ---- eval=FALSE--------------------------------------------------------------
#  gibbR=gibbsR(1000,20)
#  gibbC=gibbsC(1000,20)

## ---- eval=FALSE--------------------------------------------------------------
#  par(mfrow = c(1,2))
#  qqplot(gibbR[,1], gibbC[,1], main='Xt')
#  abline(a=0,b=1, col='red')
#  qqplot(gibbR[,2], gibbC[,2], main='Yt')
#  abline(a=0,b=1, col='red')

## ---- eval=FALSE--------------------------------------------------------------
#  library(microbenchmark)
#  compare <- microbenchmark(gibbR=gibbsR(1000,20), gibbC=gibbsC(1000,20))
#  summary(compare)

