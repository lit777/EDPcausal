####################
#### Scenario 1 ####
####################


#------ required libraries
library(Rcpp)
library(lpSolve)
library(mcclust)
library(cluster)
library(extraDistr)
library(truncnorm)
library(dplyr)
sourceCpp("clustering.cpp", rebuild=T)
source("function.R")


y <- a <- Y1 <- Y0 <- Pi0 <- Pi1 <- p0 <- p1 <- Y1M0 <- p10 <- Pi10 <- NULL


for(process in 1:40){
  tryCatch({

ll <- 0
Expit <- function(x){
  exp(x)/(1+exp(x))
}

x1 <- rbinom(600, 1, 0.7)
x2 <- rep(c(-2,-1,0,1,2,3), each=100)
x3 <- rnorm(600, 0, 1)  
x4 <- rnorm(600, 0, 1)  
x5 <- rnorm(600, 0, 1)  

# data generation
for(i in 1:600){
  a[i] <- rbinom(1,1,Expit(1-1*x1[i]+2*x2[i]-1*abs(x3[i]-1)+0.4*x4[i]*x5[i]))
  p0[i] <- Expit(-2.5+0.2*x1[i]-0.5*x3[i]-3.5*I(x2[i]>0)*0+0.2*x4[i]+0.2*x5[i])
  p1[i] <- Expit(-2.5+0.2*x1[i]-0.5*x3[i]-3.5*I(x2[i]>0)*1+0.2*x4[i]+0.2*x5[i])
  Pi0[i] <- sample(c(0,1), 1, prob=c(1-p0[i], p0[i]))
  Pi1[i] <- sample(c(0,1), 1, prob=c(1-p1[i], p1[i]))
  rn1 <- rnorm(1,  mean=0.4*x1[i]+2*x2[i]+1*x2[i]^2*1-2.5*1+0.3*abs(x3[i]+1)+0.4*x4[i]+exp(0.1*x5[i]), sd=0.1)
  rn0 <- rnorm(1,  mean=0.4*x1[i]+2*x2[i]+1*x2[i]^2*0-2.5*0+0.3*abs(x3[i]+1)+0.4*x4[i]+exp(0.1*x5[i]), sd=0.1)
  Y1[i] <- ifelse(Pi1[i]==1, 0, ifelse(rn1 < 0, 0, rn1))
  Y0[i] <- ifelse(Pi0[i]==1, 0, ifelse(rn0 < 0, 0, rn0))
}
# observed Y
Y <- a*Y1+(1-a)*Y0

# generate clusters
Sy <- sample(c(1,2,3,4,5), 600, replace=TRUE)
Sx <- sample(c(1,2),600, replace=TRUE)

uniqueYX <- distinct(tbl_df(cbind(Sy,Sx)))
uniqueYX <- uniqueYX[order(uniqueYX$Sy,uniqueYX$Sx),]

n <- length(Sy)
y <- a*Y1+(1-a)*Y0
x <- cbind(1, x1, x2, x3, x4, x5)
xa <- cbind(1, a, x1, x2, x3, x4, x5)

yPAR_beta <- matrix(rnorm(7*5, 0, 1), ncol=5, nrow=7)
yPAR_r <- matrix(rnorm(7*5, 0, 1), ncol=5, nrow=7)
yPAR_sig <- rep(1, 5)

xPAR_p0 <- matrix(rep(0.5, 2*10), ncol=2, nrow=10)
xPAR_p1 <- matrix(rep(0.5, 2*10), ncol=2, nrow=10)
xPAR_mu <- matrix(rnorm(40,0,1), ncol=10, nrow=4)
xPAR_sig <- matrix(1, ncol=10, nrow=4)

# (hyper-)parameters
Gamma_0 <- 1
a0 <- 2; b0 <- 1
alpha_theta <- 2
alpha_omega <- 2
c0 <- 1
nu0 <- 2; tau0 <- 1
alpha00 <- c(0.1,0.1)
alpha01 <- c(0.1,0.1)

mu.r=0; Sigma.r=10
mu.beta=0; Sigma.beta=10

M <- 1
mu0 <- 0

# list objects to save MCMC outputs
E <- list()
C <- list()
C[[1]] <- list(Sx=Sx, Sy=Sy)
P <- list()
P[[1]] <- list(alpha_omega=alpha_omega,  alpha_theta=alpha_theta,xPAR_p0=xPAR_p0, xPAR_p1=xPAR_p1, xPAR_mu=xPAR_mu, xPAR_sig=xPAR_sig, 
               yPAR_beta=yPAR_beta, yPAR_r=yPAR_r, yPAR_sig=yPAR_sig)

#------ Main MCMC run
for(l in 2:20000){
  C[[l]] <- clustering_c( C[[l-1]]$Sy, C[[l-1]]$Sx, 
                          P[[l-1]]$xPAR_p0, P[[l-1]]$xPAR_p1, P[[l-1]]$xPAR_mu, P[[l-1]]$xPAR_sig, 
                          P[[l-1]]$yPAR_beta, P[[l-1]]$yPAR_r, P[[l-1]]$yPAR_sig, 
                          P[[l-1]]$alpha_omega, P[[l-1]]$alpha_theta, xa, y)
  P[[l]] <- update_parameters(C[[l]]$Sx, C[[l]]$Sy, 
                              C[[l]]$xPAR_p0, C[[l]]$xPAR_p1, 
                              C[[l]]$xPAR_mu, C[[l]]$xPAR_sig, 
                              C[[l]]$yPAR_beta, C[[l]]$yPAR_r, C[[l]]$yPAR_sig)
  print(l)
}

#------- Post-processing
for(l in seq(10001, 20000, by=10)){
  ll <- ll + 1
  EY11 <- EY00 <- NULL
  numY <- length(unique(C[[l]]$Sy))
  cluster <- cbind(C[[l]]$Sy, C[[l]]$Sx)
  nj_i <- table(factor(C[[l]]$Sy, levels=1:numY));
  nlj_i <- table(C[[l]]$Sy, C[[l]]$Sx);
  numXj <- apply(nlj_i, 1, function(x) length(which(x!=0)))
  
  Nlj <- matrix(nrow=numY, ncol=2)
  Nlj[1,1:2] <- c(1, numXj[1])
  if(numY!=1){
    for(c in 2:numY){
      Nlj[c, 1:2] <- c(cumsum(numXj)[c-1]+1, cumsum(numXj)[c])
    }
  }
  
  nnlj_i <- c(t(nlj_i))
  nnlj_i <- nnlj_i[which(nnlj_i!=0)]
  
  # EY11
  ZEROS0 <- ZEROS1 <- rep(0, n)
  X1 <- cbind(1,1,x1,x2,x3,x4,x5)
  
  prob <- matrix(nrow=n, ncol=length(nj_i))
  for(clus.x in 1:length(nj_i)){
    prob[,clus.x] <- nj_i[clus.x]/(P[[l]]$alpha_theta+n)*
      (rowSums(sapply(seq(Nlj[clus.x, 1],Nlj[clus.x, 2]) , function(x) nnlj_i[x]/(P[[l]]$alpha_omega+nj_i[clus.x])*dcat(1+1,P[[l]]$xPAR_p0[x,])*dcat(X1[,3]+1,P[[l]]$xPAR_p1[x,])*dnorm(X1[,4], P[[l]]$xPAR_mu[1,x], sqrt(P[[l]]$xPAR_sig[1,x]))*dnorm(X1[,5], P[[l]]$xPAR_mu[2,x], sqrt(P[[l]]$xPAR_sig[2,x]))*dnorm(X1[,6], P[[l]]$xPAR_mu[3,x], sqrt(P[[l]]$xPAR_sig[3,x]))*dnorm(X1[,7], P[[l]]$xPAR_mu[4,x], sqrt(P[[l]]$xPAR_sig[4,x])) )))
  }
  prob <- prob/rowSums(prob)
  
  
  for(ii in 1:n){
    EY11.temp <- NULL
    ZEROS1.temp <- NULL
    for(iii in 1:length(nj_i)){
      aa=0
      bb=Inf
      sig=sqrt(P[[l]]$yPAR_sig[iii, 1])
      MU = X1[ii,]%*%P[[l]]$yPAR_beta[,iii]
      EY11.temp[iii] <-mean(rtruncnorm(1000,aa,bb,mean=MU,sd=sig) )*(1-expit(P[[l]]$yPAR_r[,iii], X1[ii,]))
      ZEROS1.temp[iii] <- expit(P[[l]]$yPAR_r[,iii], X1[ii,])
    }
    EY11[ii] <-sum(prob[ii,]*EY11.temp)
    ZEROS1[ii] <- sum(prob[ii,]*ZEROS1.temp)
  }
  
  X0 <- cbind(1,0,x1,x2,x3,x4,x5)
  
  prob <- matrix(nrow=n, ncol=length(nj_i))
  for(clus.x in 1:length(nj_i)){
    prob[,clus.x] <- nj_i[clus.x]/(P[[l]]$alpha_theta+n)*
      (rowSums(sapply(seq(Nlj[clus.x, 1],Nlj[clus.x, 2]) , function(x) nnlj_i[x]/(P[[l]]$alpha_omega+nj_i[clus.x])*dcat(1+0,P[[l]]$xPAR_p0[x,])*dcat(X0[,3]+1,P[[l]]$xPAR_p1[x,])*dnorm(X0[,4], P[[l]]$xPAR_mu[1,x], sqrt(P[[l]]$xPAR_sig[1,x]))*dnorm(X0[,5], P[[l]]$xPAR_mu[2,x], sqrt(P[[l]]$xPAR_sig[2,x]))*dnorm(X0[,6], P[[l]]$xPAR_mu[3,x], sqrt(P[[l]]$xPAR_sig[3,x]))*dnorm(X0[,7], P[[l]]$xPAR_mu[4,x], sqrt(P[[l]]$xPAR_sig[4,x])) )))
  }
  prob <- prob/rowSums(prob)
  
  
  for(ii in 1:n){
    EY00.temp <- NULL
    ZEROS0.temp <- NULL
    for(iii in 1:length(nj_i)){
      aa=0
      bb=Inf
      sig=sqrt(P[[l]]$yPAR_sig[iii, 1])
      MU = X0[ii,]%*%P[[l]]$yPAR_beta[,iii]
      EY00.temp[iii] <-mean(rtruncnorm(1000,aa,bb,mean=MU,sd=sig) )*(1-expit(P[[l]]$yPAR_r[,iii], X0[ii,]))
      ZEROS0.temp[iii] <- expit(P[[l]]$yPAR_r[,iii], X0[ii,])
    }
    EY00[ii] <-sum(prob[ii,]*EY00.temp)
    ZEROS0[ii] <- sum(prob[ii,]*ZEROS0.temp)
  }
  
  
  E[[ll]] <- list(EPI1=ZEROS1, EPI0=ZEROS0, EY11=EY11,EY00=EY00, Ey = EY11*a+(1-a)*EY00)  # individual causal effects
  
}


#------ Estimating the effects
E.final <- matrix(nrow=600, ncol=length(E))
for(l in 1:(length(E))){
  E.final[,l] <- E[[l]]$EY11-E[[l]]$EY00 
}

E.y <- matrix(nrow=600, ncol=length(E))
for(l in 1:length(E)){
  E.y[,l] <- E[[l]]$Ey 
}

E.pi <- matrix(nrow=600, ncol=length(E))
for(l in 1:(length(E))){
  E.pi[,l] <- E[[l]]$EPI1-E[[l]]$EPI0
}
print(process)
save(E.final,E,E.y,E.pi,Y1,Y0,Pi1,Pi0,p1,p0, file=paste("out_edp_v1",process,".RData",sep=""))
  }, error=function(e){})
}


