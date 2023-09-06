####################
#### Scenario 2 ####
####################

E.F <- matrix(ncol=200, nrow=600)
a<-y<-Y1<-Y0<-Pi0<-Pi1<-p0<-p1<-Y1M0<-p10<-Pi10<-NULL

load("covariate.RData")

for(test_case in 1:200){
#------ required libraries
library(bcf)
library(lpSolve)
library(mcclust)
library(cluster)
library(extraDistr)
library(truncnorm)
library(dplyr)


#------ Preparation of the simulation data / initial set-ups
  ll <- 0
  
  expit <- function(x){
    exp(x)/(1+exp(x))
  }
  
  
  for(i in 1:600){
  # version 2
    a[i] <- rbinom(1,1,expit(1-1*x1[i]+5*x2[i]-1*abs(x3[i]-1)+0.4*x4[i]*x5[i]))
    p0[i] <- expit(-3.5+0.2*x1[i]-0.5*x3[i]-3.5*I(x2[i]>0)*0+0.2*x4[i]+0.2*x5[i])
    p1[i] <- expit(-3.5+0.2*x1[i]-0.5*x3[i]-3.5*I(x2[i]>0)*1+0.2*x4[i]+0.2*x5[i])
    Pi0[i] <- sample(c(0,1), 1, prob=c(1-p0[i], p0[i]))
    Pi1[i] <- sample(c(0,1), 1, prob=c(1-p1[i], p1[i]))
    rn1 <- rnorm(1,  mean=0.4*x1[i]+2*x2[i]+1*x2[i]^2*1-2.5*1+0.3*abs(x3[i]+1)+0.4*x4[i]+exp(0.1*x5[i]), sd=0.1)
    rn0 <- rnorm(1,  mean=0.4*x1[i]+2*x2[i]+1*x2[i]^2*0-2.5*0+0.3*abs(x3[i]+1)+0.4*x4[i]+exp(0.1*x5[i]), sd=0.1)
    Y1[i] <- ifelse(Pi1[i]==1, 0, ifelse(rn1 < 0, 0, rn1))
    Y0[i] <- ifelse(Pi0[i]==1, 0, ifelse(rn0 < 0, 0, rn0))
  }
  
  pi.fit <- glm(a~x1+x2+x3+x4+x5, family=binomial)
  pihat <- predict(pi.fit, type="response")
  
  y <- a*Y1+(1-a)*Y0
  x <- cbind(1, x1, x2, x3, x4, x5)
  fit <- bcf(y, a, cbind(x1,x2,x3,x4,x5), cbind(x1,x2,x3,x4,x5), pihat, nburn=3000, nsim=3000)
  
  plot(x2, Y1-Y0, ylim=c(-10,15))
  points(x2, colMeans(fit$tau), col="red", cex=0.2)

  #------ Collect Posterior Samples      
  E.final <- apply(fit$tau, 2, mean)
  E.F[,test_case] <- E.final
  
  save(E.final,Y1,Y0, fit, file=paste("out_bcf",test_case,".RData",sep=""))
}

