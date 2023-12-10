## -----------------------------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
summary(lm.D9)$coef


## ----echo=FALSE---------------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
summary(lm.D9)$coef

## ----eval=FALSE---------------------------------------------------------------
#  ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#  trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#  group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#  weight <- c(ctl, trt)
#  lm.D9 <- lm(weight ~ group)
#  summary(lm.D9)$coef

## ----echo=FALSE,eval=FALSE----------------------------------------------------
#  ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#  trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#  group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#  weight <- c(ctl, trt)
#  lm.D9 <- lm(weight ~ group)
#  summary(lm.D9)$coef

## -----------------------------------------------------------------------------
xtable::xtable(head(iris))
a<-xtable::xtable(head(iris))
xtable::print.xtable(a,type="latex",file="D:/Datafiles/1.tex" ) #将生成的Latex代码保存
knitr::kable(head(iris)) #也可以直接用kable()生成表格

## -----------------------------------------------------------------------------
par(mfrow=c(2,2))
plot(lm.D9)

## -----------------------------------------------------------------------------
my.sample<-function(x,size,prob=1/length(x))
{u<-runif(size)
cp <- cumsum(prob)
a<-x[findInterval(u,cp)+1]
a
}

## -----------------------------------------------------------------------------
my.sample(0:1,10)
my.sample(1:4,50,c(0.1,0.2,0.3,0.4))

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
x<-numeric((1000))
x[u>=1/2]<--log(2*(1-u)) # F(x)=(e^x)/2,x<0
x[u<1/2]<-log(2*u) # F(x)=1-(e^-x)/2,x>=0      
x
hist(x, prob = TRUE,ylim =c(0, 0.5))
y <- seq(-5, 5, .01)
lines(y,exp(-abs(y))/2 )


## -----------------------------------------------------------------------------
betasample<-function(n,a,b) #这里a,b的值应该大于等于1
{j<-k<-0;y <- numeric(n)
while (k < n) {
u <- runif(1)
j <- j + 1
x <- runif(1) 
if (x^(a-1) * (1-x)^(b-1) > u) {
k <- k + 1
y[k] <- x
}
}
y
} 
z<-betasample(1000,3,2)
z
hist(z, prob = TRUE,ylim =c(0, 2))
x <- seq(0, 1, .01)
lines(x,dbeta(x,3,2))

## -----------------------------------------------------------------------------
fesample<-function(n){
x<-numeric(n)
i<-1
for(i in 1:n){
u<-runif(3,-1,1)
if(abs(u[3])>=abs(u[2])&abs(u[3])>=abs(u[1]))
{x[i]<-u[2]}
else{x[i]<-u[3]}
i<-i+1
}
x
}
y<-fesample(1000)
y
hist(y, prob = TRUE,ylim =c(0, 1))
x <- seq(-1, 1, .01)
lines(x,3*(1-x^2)/4 )

## -----------------------------------------------------------------------------
rho1<-1
rho2<-1/2
rho3<-1/3
m <- 1e6
pihat1<-numeric(100) #rho=1时的模拟结果
pihat2<-numeric(100) #rho=1/2时的模拟结果
pihat3<-numeric(100) #rho=1/3时的模拟结果
for(i in 1:100)
{
X <- runif(m,0,1/2)
Y <- runif(m,0,pi/2)
pihat1[i]<- 2*rho1/mean(rho1/2*sin(Y)>X)
pihat2[i] <- 2*rho2/mean(rho2/2*sin(Y)>X)
pihat3[i] <- 2*rho3/mean(rho3/2*sin(Y)>X)
}
var(pihat1)
var(pihat2)
var(pihat3)

## -----------------------------------------------------------------------------
n<-1e4
var1<-1/n*(2*exp(1)-1/2*exp(2)-3/2)
var2<-1/2/n*(10*exp(1)-3*exp(2)-5)
a<-(var1-var2)/var1

## -----------------------------------------------------------------------------
n<-1e4
x<-runif(n)
theta.hat<-mean(exp(x)) #simple Monte Carlo method
theta.hat1<-mean(exp(x[1:n/2])+exp(1-x[1:n/2]))/2 #antithetic variate approach
var1<-1/n*var(exp(x))
var2<-var(exp(x[1:n/2])+exp(1-x[1:n/2]))/2/n
a<-(var1-var2)/var1

## -----------------------------------------------------------------------------
f<-function(x){x^2*exp(-(x^2)/2)/sqrt(2*pi)}
f1<-function(x){sqrt(exp(1))*x*exp(-(x^2)/2)}
n<-1e5
u<-runif(n)
x<-sqrt(1-2*log(1-u))#用 Inverse transform method生成密度f1的随机数
est1<-mean(x/sqrt(2*pi*exp(1)))
var1<-1/n*var(x/sqrt(2*pi*exp(1)))
y<-rnorm(ceiling(n/0.1589*10))
y <- y[y > 1]
y <- y[1:n]#生成标准正态随机数并保留大于1的
est2<-mean(0.1589*y^2)
var2<-1/n*var(0.1589*y^2)

## -----------------------------------------------------------------------------
M <- 10000; k <- 5
T<-numeric(k)
a<-numeric(k)
var<-numeric(k)
g<-function(x){exp(-x)/(1+x^2)*(x>0)*(x<1)}
for(j in 1:k){
  u<-runif(M/k,(j-1)/k,j/k)
  a[j]<-(exp((1-j)/k)-exp(-j/k))/(1-exp(-1))
  x<--log(exp((1-j)/k)-(1-exp(-1))*u*a[j])
  fg<-g(x)/(exp(-x)/(1-exp(-1)))
  T[j]<-mean(fg)
  var[j]<-var(fg)
}
est1<-a%*%T
sd1<-sqrt(a^2%*%var*k/M)#stratified importance sampling estimate and its sd

u <- runif(M)
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
est2<- mean(fg)
sd2<- sd(fg)/sqrt(M)#importance sampling estimate and its sd

## -----------------------------------------------------------------------------
m<-1000
n<-20
L<-numeric(m)
U<-numeric(m)
UCL<-numeric(m)
for(i in 1:m){
x<-rchisq(n,2)
L[i]<-mean(x)-sd(x)/sqrt(n)*qt(0.975,n-1)
U[i]<-mean(x)+sd(x)/sqrt(n)*qt(0.975,n-1)
}
mean(L<=2 & 2<=U)
for(i in 1:m){
x<-rchisq(n,2)
UCL[i] <- (n-1) * var(x) / qchisq(0.05, df=n-1)
}
mean(4<=UCL)

## -----------------------------------------------------------------------------
mu0<-1
alpha<-0.05
n<-20
m<-1e4
p<-numeric(m)
for(i in 1:m){
x<-rchisq(n,1)
p[i]<-t.test(x,mu=mu0)$p.value
}
mean(p<=alpha)

for(i in 1:m){
x<-runif(n,0,2)
p[i]<-t.test(x,mu=mu0)$p.value
}
mean(p<=alpha)

for(i in 1:m){
x<-rexp(n,1)
p[i]<-t.test(x,mu=mu0)$p.value
}
mean(p<=alpha)

## -----------------------------------------------------------------------------
alpha<-0.1
m<-1000
M<-1000
p<-numeric(m)
p.adj1<-matrix(0,m,M)
p.adj2<-matrix(0,m,M)
FWER<-matrix(0,M,2)
FDR<-matrix(0,M,2)
TPR<-matrix(0,M,2)

for(i in 1:M){
 p1<-runif(.95*m) #negative
 p2<-rbeta(.05*m,0.1,1) #positive
 p<-sort(c(p1,p2))
 p.adj1[i,] <- p.adjust(p,method='fdr')
 p.adj2[i,] <- p.adjust(p,method='bonferroni')
 FWER[i,1]<-sum(p.adj1[i,]<=alpha&p%in%p1)!=0
 FWER[i,2]<-sum(p.adj2[i,]<=alpha&p%in%p1)!=0
 FDR[i,1]<-sum(p.adj1[i,]<=alpha&p%in%p1)/sum(p.adj1[i,]<=alpha)
 FDR[i,2]<-sum(p.adj2[i,]<=alpha&p%in%p1)/sum(p.adj2[i,]<=alpha)
 TPR[i,1]<-(sum(p.adj1[i,]<=alpha)-sum(p.adj1[i,]<=alpha&p%in%p1))/.05/m
 TPR[i,2]<-(sum(p.adj2[i,]<=alpha)-sum(p.adj2[i,]<=alpha&p%in%p1))/.05/m
}
a<-data.frame(FWER=round(apply(FWER,2,mean),3),FDR=round(apply(FDR,2,mean),3),TPR=round(apply(TPR,2,mean),3),row.names = c("B-H","Bonferroni"))

knitr::kable (a)



## -----------------------------------------------------------------------------
lambda<-2
n<-c(5,10,20) 
B<-1000 
m<- 1000  
lambda.hat<-matrix(0,m,length(n))
bias<-matrix(0,m,length(n))
se<-matrix(0,m,length(n))
MLE<-function(x,i){
 return(1/mean(x[i]))
}
for(j in 1:length(n)){ 
 for(i in 1:m){
  x<-rexp(n[j],lambda)
  obj<-boot::boot(data=x,statistic=MLE,R=B)
  lambda.hat[i,j]<-obj$t0
  bias[i,j]<-mean(obj$t)-obj$t0
  se[i,j]<-sd(obj$t)
 }
}
a<-data.frame(bias_theoretical=round(lambda/(n-1),3),bias_bootstrapmean=round(apply(bias,2,mean),3),se_theoretical=round(n*lambda/(n-1)/sqrt(n-2),3),se_bootstrap=round(apply(se,2,mean),3),row.names = c("n=5","n=10","n=20"))
knitr::kable (a)

## -----------------------------------------------------------------------------
library(bootstrap)
B <- 500
n <- nrow(law)
R <- 100
alpha<-.05
cor<-numeric(B)
cor2<-numeric(R)
t<-numeric(B)
corhat<-cor(law$LSAT, law$GPA)
for (b in 1:B) {
 x <- sample(1:n, size = n, replace = TRUE)
 LSAT <- law$LSAT[x] 
 GPA <- law$GPA[x]
 cor[b] <- cor(LSAT, GPA)
 for(i in 1:R){
   y<- sample(x, size = n, replace = TRUE)
   LSAT <- law$LSAT[y] 
   GPA <- law$GPA[y]
   cor2[i] <- cor(LSAT, GPA)
 }
 t[b]<-(cor[b]-corhat)/sd(cor2) 
}
t1<-corhat-sort(t)[ceiling(B*(1-alpha/2))]*sd(cor)
t2<-corhat-sort(t)[floor(B*(alpha/2))]*sd(cor)

## -----------------------------------------------------------------------------
library(boot)
set.seed(123)
R<-1000
b.mean<-function(x,i){
  mean(x[i])
}
obj<-boot(aircondit$hours,b.mean,R)
ci<-boot.ci(obj,type=c("norm","basic","perc","bca"))
ci

## -----------------------------------------------------------------------------
library(bootstrap)
x<-as.matrix(scor)
n<-nrow(scor)
m<-ncol(scor)
b.cov<-function(x,i){
  cov(x[i,])
}
COV<-matrix(0,m,m)
theta.hat <-eigen(cov(x))$values[1]/sum(eigen(cov(x))$values)
theta.jack <- numeric(n)
for(i in 1:n){
COV <- b.cov(x,(1:n)[-i])
theta.jack[i] <- eigen(COV)$values[1]/sum(eigen(COV)$values)
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,
se.jack=se.jack),3)


## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- matrix(0,n,n)
for (i in 2:n) {
  for(j in 1:(i-1)){
  
  y <- magnetic[-c(i,j)]
  x <- chemical[-c(i,j)]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[i]
  yhat2 <- J1$coef[1] + J1$coef[2] * chemical[j]
  e1[i,j] <- magnetic[i] - yhat1
  e1[j,i] <- magnetic[j] - yhat2
  
  
  J2 <- lm(y ~ x + I(x^2))
  yhat3 <- J2$coef[1] + J2$coef[2] * chemical[i] +J2$coef[3] * chemical[i]^2
  yhat4 <- J2$coef[1] + J2$coef[2] * chemical[j] +J2$coef[3] * chemical[j]^2
  e2[i,j] <- magnetic[i] - yhat3
  e2[j,i] <- magnetic[j] - yhat4

  J3 <- lm(log(y) ~ x)
  logyhat5 <- J3$coef[1] + J3$coef[2] * chemical[i]
  logyhat6 <- J3$coef[1] + J3$coef[2] * chemical[j]
  yhat5 <- exp(logyhat5)
  yhat6 <- exp(logyhat6)
  e3[i,j] <- magnetic[i] - yhat5
  e3[j,i] <- magnetic[j] - yhat6
  
  J4 <- lm(log(y) ~ log(x))
  logyhat7 <- J4$coef[1] + J4$coef[2] * log(chemical[i])
  logyhat8 <- J4$coef[1] + J4$coef[2] * log(chemical[j])
  yhat7 <- exp(logyhat7)
  yhat8 <- exp(logyhat8)
  e4[i,j] <- magnetic[i] - yhat7
  e4[j,i] <- magnetic[j] - yhat8
  }
}
 c(sum(diag(e1%*%t(e1))),sum(diag(e2%*%t(e2))), sum(diag(e3%*%t(e3))),sum(diag(e4%*%t(e4))))/n/(n-1)


## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
R <- 999 
z <- c(x, y) 
K <- 1:26
W <- numeric(R)
cvm<-function(x,y){
  n<-length(x)
  m<-length(y)
  a<-numeric(n)
  b<-numeric(m)
  for(i in 1:n) a[i]<-(sum(x<=x[i])/n-sum(y<=x[i])/m)^2
  for(i in 1:m) b[i]<-(sum(x<=y[i])/n-sum(y<=y[i])/m)^2
  return(n*m*(sum(a)+sum(b))/(m+n)^2)
}
W0 <- cvm(x, y)
for (i in 1:R) {
k <- sample(K, size = 14, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] 
W[i] <- cvm(x1, y1)
}
p <- mean(c(W0, W) >= W0)
p


## -----------------------------------------------------------------------------
maxout<- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy)))
}
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 10000 #模拟实验重复次数
R<-999 #重排法的重排次数
K<-1:50
M<-numeric(R)

alphahat <- mean(replicate(m, expr={
  x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) 
y <- y - mean(y)
z<-c(x,y)
M0<-maxout(x,y)
for (i in 1:R) {
k <- sample(K, size = n1, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] 
M[i] <- maxout(x1, y1)
}
p <- mean(c(M0, M) >= M0)
p<=.05
}))
alphahat

## -----------------------------------------------------------------------------
N <- 1e6; b1 <-0; b2 <- 1; b3<--1;f0 <- c(.1,.01,.001,.0001)
a<-numeric(length(f0))
g <- function(N,b1,b2,b3,f0){
  x1 <- rpois(N,1)
  x2<-rexp(N,1)
  x3 <- rbinom(N,1,0.5)
  f<-function(a){
    tmp <- exp(-a-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
  }
  solution <- uniroot(f,c(-100,100))
  return(round(unlist(solution)[1],5))
}
for(i in 1:length(f0)){
  a[i]<-g(N,b1,b2,b3,f0[i])
}
a
plot(-log(f0),a)

## -----------------------------------------------------------------------------
 rw.Metropolis <- function( sigma, x0, N) {
       # sigma:  standard variance of proposal distribution N(xt,sigma)
       # x0: initial value
       # N: size of random numbers required.
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
                if (u[i] <= (exp(-abs(y))) / (exp(-abs(x[i-1])))){
                    x[i] <- y 
                    k <- k + 1
                }
                else 
                    x[i] <- x[i-1]
                  
                
            }
        return(list(x=x, k=k))
        }

    N <- 1e4
    sigma <- c(.05, .5, 2,  16)

    x0 <- 10 
    rw1 <- rw.Metropolis(sigma[1], x0, N)
    rw2 <- rw.Metropolis(sigma[2], x0, N)
    rw3 <- rw.Metropolis(sigma[3], x0, N)
    rw4 <- rw.Metropolis(sigma[4], x0, N)
    print(c(rw1$k, rw2$k, rw3$k, rw4$k)/N) #acceptance rates
    
    par(mfrow=c(2,2))
    rw<-cbind(rw1$x,rw2$x,rw3$x,rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",xlab=bquote(sigma == .(round(sigma[j],3))),ylab="X", ylim=range(rw[,j]))
    }
    par(mfrow=c(1,1))


## -----------------------------------------------------------------------------
  #initialize constants and parameters
    N <- 1e4               #length of chain
    burn <- 2000            #burn-in length
    X <- numeric(N)        #the chain
    Y <- numeric(N)
    rho <-.9             #correlation
    mu1 <- 0
    mu2 <- 0
    sigma1 <- 1
    sigma2 <- 1
    s1 <- sqrt(1-rho^2)*sigma1
    s2 <- sqrt(1-rho^2)*sigma2

    ###### generate the chain #####

    X[1] <- mu1            #initialize
    Y[1] <- mu2
    for (i in 2:N) {
        x2 <- Y[i-1]
        m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
        X[i] <- rnorm(1, m1, s1)
        x1 <- X[i]
        m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
        Y[i] <- rnorm(1, m2, s2)
    }
    b <- burn + 1
    X<- X[b:N]
    Y<- Y[b:N]
    plot(X,Y, main="", cex=.5, xlab=bquote(X),
         ylab=bquote(Y), ylim=range(Y),xlim=range(X))
    
par(mfrow=c(2,2))    
fit<-lm(Y~X,as.data.frame(cbind(Y,X)))
summary(fit)
plot(fit)

## -----------------------------------------------------------------------------
set.seed(12138)
Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + B/n+(B/(n*k))     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
}

f <- function(x, sigma) {
        if (any(x < 0)) return (0)
        stopifnot(sigma > 0)
        return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

Rayleigh.chain<-function(R,sigma,k){
    x <- matrix(0, nrow=1, ncol=k)
    R.hat<-numeric(1)
    for (i in 1:k) {
      x[1,i] <- rchisq(1, df=1)
    }
    R.hat[1]<-10
    
                                                                                           
      t=2
      while (R.hat[t-1]>=1.2) {
         x<-rbind(x,0)
         R.hat<-rbind(R.hat,0)
        for(i in 1:k){
          xt <- x[t-1,i]
          y <- rchisq(1, df = xt)
          num <- f(y, sigma) * dchisq(xt, df = y)
          den <- f(xt, sigma) * dchisq(y, df = xt)
          u <- runif(k)
          if (u[i] <= num/den) x[t,i] <- y 
          else  x[t,i] <- xt
        }
        psi <- t(apply(x, 2, cumsum))
        for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))
        R.hat[t]<-Gelman.Rubin(psi)
        t<-t+1
      }

    return(list(GR.statistic=R.hat[t-1],chain=x,length=t-1) )
 
  
}
    sigma <- 4     #parameter of proposal distribution
    k <- 4          #number of chains to generate
    R <- 1.2     #target value of Gelman-Rubin statistic 
    

a<-Rayleigh.chain(R,sigma,k)
a$GR.statistic


###### coda  package #####
library(coda)
b<-mcmc.list(chain1=mcmc(a$chain[,1]),chain2=mcmc(a$chain[,2]),chain3=mcmc(a$chain[,3]),chain4=mcmc(a$chain[,4]))
gelman.diag(b)
gelman.plot(b)

## -----------------------------------------------------------------------------
u<-c(11,8,27,13,16,0,23,10,24,2)
v<-c(12,9,28,14,17,1,24,11,25,3)
n<-10
f<-function(lambda){
  a<--exp(-lambda*u)/(exp(-lambda*u)-exp(-lambda*v))
  b<-exp(-lambda*v)/(exp(-lambda*u)-exp(-lambda*v))
  t(a)%*%u+t(b)%*%v
}
solution <- uniroot(f,c(.01,10))
MLE <- solution$root

em<-function(lambda,m){
  a<-numeric(m)
  a[1]<-lambda
  for(i in 2:m){
    b<-sum((u*exp(-a[i-1]*u)-v*exp(-a[i-1]*v))/(exp(-a[i-1]*u)-exp(-a[i-1]*v))+1/a[i-1])
    a[i]<-n/b
  }
  return(a[m])
}
EM<-em(MLE,1e5)

round(c(MLE,EM),5)

## -----------------------------------------------------------------------------
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
library(boot) 
B<-A+2

m<-nrow(B)
n<-ncol(B)
a <- c(rep(0, m), 1) 
A1 <- -cbind(t(B), rep(-1, n)) 
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) 
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE)
sx$soln

## -----------------------------------------------------------------------------
a<-data.frame("x"=NULL)
a

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}
a<-data.frame(x=c(1,2,3),y=c("a","b","c"))
b<-lapply(a,is.numeric)
scale01(a[which(unlist(b))])

## -----------------------------------------------------------------------------
a<-data.frame(x=c(1,2,3),y=c(4,5,6))
vapply(a,sd,numeric(1))

b<-data.frame(x=c(1,2,3),y=c("a","b","c"))
c<-vapply(b,is.numeric,logical(1))
vapply(b[which(c)],sd,numeric(1))

## -----------------------------------------------------------------------------
N <- 5000 
a<-1
b<-1
n<-5  
Gibbs<-function(N,a,b,n){
  X <- matrix(0,N,2)
  X[1, ] <- c(0, 0) 
  for (i in 2:N) {
    y <- X[i-1, 2]
    X[i, 1] <- rbinom(1, n, y)
    x <- X[i, 1]
    X[i, 2] <- rbeta(1, x+a, n-x+b)
  }
return(X)
}


## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(gibbR=Gibbs(N,a,b,n),gibbC=GibbsCpp(N,a,b,n))
summary(ts)[,c(1,3,5,6)]

