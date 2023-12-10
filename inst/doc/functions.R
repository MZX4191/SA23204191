## -----------------------------------------------------------------------------
library(Rcpp)
sourceCpp(file="D:/STATCOMP/SA23204191/src/rmnGibbsc.cpp")
Y <- rmnGibbsc(50, c(1, 2,3 ),c(0, 0,0),matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),3,3))
Y

## -----------------------------------------------------------------------------
rmnGibbsR<-function(m,x0,mu,sigma){
  if(1-isSymmetric.matrix(sigma)){
    sigma<-(sigma+t(sigma))/2
  }
  n<-length(mu)
  x<-matrix(0,n,m)
  A<-array(0,dim=c(n,n,n))
  A[n,1,1]<-1;A[1:(n-1),2:n,1]<-diag(n-1)
  A[,,n]<-diag(n)
  if(n>=3){
    for(i in 2:(n-1)){
    A[1:(i-1),1:(i-1),i]<-diag(i-1)
    A[n,i,i]<-1
    A[i:(n-1),(i+1):n,i]<-diag(n-i)
    }
  }
  
  x[,1]<-x0
  for(i in 2:m){
    x[,i]<-x[,i-1]
    for(j in 1:n){
      a<-(A[,,j]%*%mu)[n]+(A[,,j]%*%sigma%*%t(A[,,j]))[n,1:(n-1)]%*%solve((A[,,j]%*%sigma%*%t(A[,,j]))[1:(n-1),1:(n-1)])%*%(x[-j,i]-(A[,,j]%*%mu)[1:(n-1)])
      b<-(A[,,j]%*%sigma%*%t(A[,,j]))[n,n]-t((A[,,j]%*%sigma%*%t(A[,,j]))[n,1:(n-1)])%*%solve((A[,,j]%*%sigma%*%t(A[,,j]))[1:(n-1),1:(n-1)])%*%(A[,,j]%*%sigma%*%t(A[,,j]))[1:(n-1),n]
      if(b<0){
       b<-0
      }
      x[j,i]<-rnorm(1,a,b)
    }
  }
  return(x)
}

BQRNN<-function(Y,X,k,tau=0.5,prior,stepsize,m,initialvalue,burnin,l){
  #prior (list) 包含beta_0   sigma_0^2  gamma_j0   sigma_1^2   a   b  
  #         beta~N(beta_0,sigma_0^2 I_k+1)
  #         gamma_j~N(gamma_j0,sigma_1^2 I_p+1)
  #         sigma~IG(a/2,b/2)相互独立
  #stepsize指MH算法中正态提议分布的方差
  #m指mcmc 链长
  #initialvalue (list)  beta,gamma的初始值
  #burn 预烧期
  #l指在链中每L个值取1个
  Y<-as.matrix(Y)
  X<-as.matrix(X)
  n<-nrow(X)
  p<-ncol(X)
  ##先验
  beta_0<-prior[[1]]  #k+1
  sigma_0<-prior[[2]] #1
  gamma_0<-prior[[3]] #p+1*k
  sigma_1<-prior[[4]] #1
  a<- prior[[5]]      #1
  b<- prior[[6]]      #1
  
  beta<-matrix(0,k+1,m)
  gamma<-array(0,dim=c(p+1,k,m))
  psi<-function(x){
    1/(1+exp(-x))
  }
  sigma<-numeric(m)
  v<-matrix(0,n,m)
  ## 初始化
  beta[,1]<-initialvalue[[1]]
  gamma[,,1]<-initialvalue[[2]]
  sigma[1]<-b/(a-2)
  v[,1]<- b/(a-2)/tau/(1-tau)
  #gamma 接受率
  acc.prob<-matrix(0,p+1,k)
  logb<-array(1,dim=c(p+1,k,m))
  
  
  for(i in 2:m){
    beta[,i]<-beta[,i-1]
    gamma[,,i]<-gamma[,,i-1]
    sigma[i]<-sigma[i-1]
    v[,i]<-v[,i-1]
    
    #beta
    x<-cbind(1,X)
    L<-cbind(1,psi(x %*% gamma[,,i]))
    V<-diag(1/v[,i])
    mu<-solve((t(L)%*%V%*%L)/2/sigma[i]+diag(k+1)/sigma_0^2)%*%((t(L)%*%V%*%(Y-(1-2*tau)*v[,i]))/2/sigma[i]+beta_0/sigma_0^2)
    c<-solve((t(L)%*%V%*%L)/2/sigma[i]+diag(k+1)/sigma_0^2)
    beta[,i]<-rmnGibbsR(2,beta[,i],mu,c)[,2]
    for(j in 1:(k+1)){
      if(is.na(beta[j,i])){
        
        return(list(beta[,i-1],mu,c,i))
        
        stop("beta")
      }
    }
    
    #gamma  MH
    
    for(j in 1:k){
      for(t in 1:(p+1)){
        prop.gamma<-gamma[,,i]
        prop.gamma[t,j]<-rnorm(1,prop.gamma[t,j],stepsize)
        prop.L<-cbind(1,psi(x %*% prop.gamma))
        L<-cbind(1,psi(x %*% gamma[,,i]))
        loga<-(t(Y-L%*%beta[,i]-(1-2*tau)*v[,i])%*%V%*%(Y-L%*%beta[,i]-(1-2*tau)*v[,i]))/4/sigma[i]+t(gamma[,j,i]-gamma_0[,j])%*%(gamma[,j,i]-gamma_0[,j])/2/sigma_1^2-(t(Y-prop.L%*%beta[,i]-(1-2*tau)*v[,i])%*%V%*%(Y-prop.L%*%beta[,i]-(1-2*tau)*v[,i]))/4/sigma[i]-t(prop.gamma[,j]-gamma_0[,j])%*%(prop.gamma[,j]-gamma_0[,j])/2/sigma_1^2
        logb[t,j,i]<-loga
        if(is.na(loga)){
          print(c(i,j,t))
          stop("loga")
        }
        u<-runif(1)
        u<-log(u)
        if(u<loga){
          gamma[t,j,i]<-prop.gamma[t,j]
          acc.prob[t,j]<-acc.prob[t,j]+1
        }
        if(is.na(gamma[t,j,i])){
          print(c(i,j,t))
          stop("gamma")
        }
      }
    }
    #sigma
    L<-cbind(1,psi(x %*% gamma[,,i]))
    shape<-(3*n+a)/2
    scale<-(t(Y-L%*%beta[,i]-(1-2*tau)*v[,i])%*%V%*%(Y-L%*%beta[,i]-(1-2*tau)*v[,i]))/4+tau*(1-tau)*sum(v[,i])+b/2
    if(is.na(scale)){
      print(i)
      stop("scale")
    }
    sigma[i]<-actuar::rinvgamma(1,shape,scale)
    if(is.na(sigma[i])){
      print(i)
      stop("sigma")
    }
    #v
    for(j in 1:n){
     nu<-1/2
     rho1<-((Y[j]-L[j,]%*%beta[,i])^2)/2/sigma[i]
     if(is.na(rho1)){
       print(c(sigma[i],i))
       stop("rho1")
     }
     rho2<-1/2/sigma[i]
     if(is.na(rho2)){
       print(c(sigma[i],i))
       stop("rho2")
     }
     v[j,i]<-GIGrvg::rgig(1,nu,rho1,rho2)   
          
        
        
    }
  }
  #output
 
  beta<-beta[,(burnin+1):m]
  d<-seq(1,ncol(beta),by=l)
  beta<-beta[,d]
  beta.hat<-apply(beta,1,mean)
  library(coda)
  Beta<-list(0)
  for(i in 1:(k+1)){
    Beta[[i]]<-beta[i,]
  }
  Beta<-lapply(Beta,mcmc)
  GRstatistic<-gelman.diag(as.mcmc.list(Beta))$psrf[,1]
  
  gamma<-gamma[,,(burnin+1):m]
  gamma<-gamma[,,d]
  gamma.hat<-apply(gamma,c(1,2),mean)
  acc.prob<-acc.prob/m
  
  return(list(beta.hat,GRstatistic,gamma.hat,acc.prob))
}


##simulations
##1
set.seed(1)
x<-matrix(runif(600,0,5),200,3)
beta1<-c(2,4,6)
beta2<-c(0.1,0.3,0.5)
e<-rnorm(200,0,1)
y<-x%*%beta1+x%*%beta2*e
k<-4
p<-3
prior<-list(rep(0,k+1),10,matrix(0,p+1,k),10,3,0.1)
stepsize<-0.01^2
m<-100000
initialvalue<-list(rep(0,k+1),matrix(0,p+1,k))
burnin<-m/2
l<-10
a<-BQRNN(y,x,k,0.05,prior,stepsize,m,initialvalue,burnin,l)
a

