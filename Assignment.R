##  Part I: Writing complex functions in R
### Task 1: Functions to estimate regression models using the SSVS prior

#Loading required libraries
library(BMS)
library(ggplot2)
library(reshape2)

#Loading data
Y<-matrix(datafls$y)
X<-as.matrix(datafls[2:42])

bayes <- function(X,Y,nsave=1000,nburn=1000, tau0=0.01, tau1=10) {
  #Start with actual estimation and Gibbs preliminaries
  ntot <- nsave+nburn
  #Prior prelims for SSVS
  #tau0 <- 0.01
  #tau1 <- 10
  
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

## Testing different values for posterior inclusion probabilities (PIP.mean)
pip1 <- bayes(X,Y)[,1] #tau0=0.01, tau1=10
pip2 <- bayes(X,Y,tau0=0.001,tau1=50)[,1]
pip3 <- bayes(X,Y,tau0=0.00001,tau1=1000)[,1]
pip4 <- bayes(X,Y,tau0=1*10^-15)[,1]
pip_df <- as.data.frame(cbind(scale=1:41,pip1,pip2,pip3,pip4))
pip_df <- melt(pip_df, id.vars='scale',variable.name = 'series')
pip_plot <- ggplot(pip_df, aes(scale,value)) + geom_line(aes(colour = series))+ylab("PIP mean")+xlab("")+ggtitle("Visualizing PIP mean for different values of tau0 and tau1")
pip_plot

pip_df <- as.data.frame(cbind(scale=1:41,bayes(X,Y)))
X_norm <- scale(X)
Y_norm <- scale(Y)
pip_df_stand <- as.data.frame(cbind(scale=1:41,bayes(X_norm,Y_norm)))
pip_plot_diff <- ggplot(pip_df_stand, aes(scale, PIP.mean)) +
  geom_line() +
  geom_line(data=pip_df, aes(y=PIP.mean), colour='red')
pip_plot_diff

matr_norm <- matrix(cbind(mean(X_norm),mean(Y_norm),sd(X_norm),sd(Y_norm)),nrow=2)
colnames(matr_norm) <- c("Std","Mean")
rownames(matr_norm) <- c("X","Y")
matr_norm

matr <- matrix(cbind(mean(X),mean(Y),sd(X),sd(Y)),nrow=2)
colnames(matr) <- c("Std","Mean")
rownames(matr) <- c("X","Y")
matr
