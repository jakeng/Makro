##  Part I: Writing complex functions in R
### Task 1: Functions to estimate regression models using the SSVS prior

#Loading the data
library(BMS)
#data(datafls)
Y<-matrix(datafls$y)
X<-as.matrix(datafls[2:42])

bayes <- function(X,Y,nsave=1000,nburn=1000) {
  #Start with actual estimation and Gibbs preliminaries
  ntot <- nsave+nburn
  #Prior prelims for SSVS
  tau0 <- 0.01
  tau1 <- 10
  
  s0 <- 0.01
  S0 <- 0.01
  #Construct Y and X
  N <- nrow(Y)
  K <- ncol(X)
  
  #Get OLS quantities
  A.OLS <- solve(crossprod(X))%*%crossprod(X,Y)
  SSE <- crossprod(Y-X%*%A.OLS)
  SIG.OLS <- SSE/(N-K)
  #In the next step, create storage matrices for Gibbs loop and initialize prior
  gamma <- matrix(1,K,1) #indicators, start with full model
  sigma2.draw <- as.numeric(SIG.OLS)
  V.prior <- diag(as.numeric(gamma*tau1+(1-gamma)*tau0))
  
  ALPHA.store <- matrix(NA,nsave,K)
  SIGMA.store <- matrix(NA,nsave,1)
  Gamma.store <- matrix(NA,nsave,K)
  
  for (irep in 1:ntot){
    #Draw ALPHA given rest from multivariate normal
    V.post <- solve(crossprod(X)*1/sigma2.draw+diag(1/diag(V.prior)))
    A.post <- V.post%*%(crossprod(X,Y)*1/sigma2.draw)
    A.draw <- A.post+t(chol(V.post))%*%rnorm(K)
    
    #Draw indicators conditional on ALPHA
    for (jj in 1:K){
      p0 <- dnorm(A.draw[[jj]],0,sqrt(tau0))
      p1 <- dnorm(A.draw[[jj]],0,sqrt(tau1))
      p11 <- p1/(p0+p1)
      
      if (p11>runif(1)) gamma[[jj]] <- 1 else gamma[[jj]] <- 0
    }
    #Construct prior VC matrix conditional on gamma
    V.prior <- diag(as.numeric(gamma*tau1+(1-gamma)*tau0))
    
    #Simulate sigma2 from inverse Gamma
    S.post <- crossprod(Y-X%*%A.draw)/2+S0
    s.post <- S0+N/2
    sigma2.draw <- 1/rgamma(1,s.post,S.post)  
    
    if (irep>nburn){
      ALPHA.store[irep-nburn,] <- A.draw
      SIGMA.store[irep-nburn,] <- sigma2.draw
      Gamma.store[irep-nburn,] <- gamma
    }
    print(irep)
  }
  #Calculate posterior inclusion probabilities
  PIP.mean <- apply(Gamma.store,2,mean)
  A.mean <- apply(ALPHA.store,2,mean)
  SIG.mean <- apply(SIGMA.store,2,mean)
  
  return(cbind(PIP.mean,A.mean,SIG.mean))
}

bayes(X,Y)
