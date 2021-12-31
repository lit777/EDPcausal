####################
#### Scenario 1 ####
####################

#------ required libraries
library(lpSolve)
library(mcclust)
library(cluster)
library(extraDistr)
library(truncnorm)
library(dplyr)

#------ Preparation of the simulation data / initial set-ups
  BIN.final1 <- list()
  BIN.final2 <- list()
  BIN.final3 <- list()
  BIN.final4 <- list()
  BIN.final5 <- list()
  BIN.final6 <- list()
  ll <- 0
  
  x1<-x2<-x3<-x4<-x5<-a<-y<-Y1<-Y0<-Pi0<-Pi1<-p0<-p1<-Y1M0<-p10<-Pi10<-NULL
  
  cluster.x <- c(rep(1,150), rep(2,150), rep(3,150), rep(4,150))
  cluster.y <- c(rep(1,150), rep(2,150), rep(3,150), rep(4,150))
  
  Expit <- function(x){
    exp(x)/(1+exp(x))
  }
  
  for(i in 1:600){
    if(cluster.x[i]==1){
      x1[i] <- sample(c(0,1), 1, replace=FALSE,prob=c(0.2,0.8))
      cov <- matrix(0.1, nrow=4, ncol=4)
      diag(cov) <- rep(0.2,4)
      x.temp <- mnormt::rmnorm(1, c(0.5,0.5,0.2,5), cov)
      x2[i] <- x.temp[1]
      x3[i] <- x.temp[2]
      x4[i] <- x.temp[3]
      x5[i] <- x.temp[4]
    }
    
    if(cluster.x[i]==2){
      x1[i] <- sample(c(0,1), 1, replace=FALSE,prob=c(0.8,0.2))
      cov <- matrix(0.1, nrow=4, ncol=4)
      diag(cov) <- rep(0.2,4)
      x.temp <- mnormt::rmnorm(1, c(1.5,2.5,0.7,0.5), cov)
      x2[i] <- x.temp[1]
      x3[i] <- x.temp[2]
      x4[i] <- x.temp[3]
      x5[i] <- x.temp[4]
    }
    
    if(cluster.x[i]==3){
      x1[i] <- sample(c(0,1), 1, replace=FALSE,prob=c(0.8,0.2))
      cov <- matrix(0.1, nrow=4, ncol=4)
      diag(cov) <- rep(0.2,4)
      x.temp <- mnormt::rmnorm(1, c(1.5,2.5,0.5,-1), cov)
      x2[i] <- x.temp[1]
      x3[i] <- x.temp[2]
      x4[i] <- x.temp[3]
      x5[i] <- x.temp[4]
    }
    
    
    if(cluster.x[i]==4){
      x1[i] <- sample(c(0,1), 1, replace=FALSE,prob=c(0.2,0.8))
      cov <- matrix(0.1, nrow=4, ncol=4)
      diag(cov) <- rep(0.2,4)
      x.temp <- mnormt::rmnorm(1, c(2.5,4.5,0.4,1), cov)
      x2[i] <- x.temp[1]
      x3[i] <- x.temp[2]
      x4[i] <- x.temp[3]
      x5[i] <- x.temp[4]
    }
  }
  
  for(i in 1:600){
    if(cluster.y[i]==1){
      a[i] <- rbinom(1,1,Expit(-1.5-1*x1[i]+0.5*x3[i]-1*x2[i]-0.5*x4[i]+0.5*x5[i]))
      p0[i] <- Expit(-3-1*0-0.2*x1[i]-1*x3[i]-0.5*x2[i]+0.5*x4[i]+0.5*x5[i])
      p1[i] <- Expit(-3-1*1-0.2*x1[i]-1*x3[i]-0.5*x2[i]+0.5*x4[i]+0.5*x5[i])
      Pi0[i] <- sample(c(0,1), 1, prob=c(1-p0[i], p0[i]))
      Pi1[i] <- sample(c(0,1), 1, prob=c(1-p1[i], p1[i]))
      Y1[i] <- ifelse(Pi1[i]==1, 0, rtruncnorm(1, a=0, mean=1+1*1-2*x1[i]+0.8*x2[i]*x3[i]+0.5*x4[i]^2+1.5*x1[i]*x4[i], sd=0.2))
      Y0[i] <- ifelse(Pi0[i]==1, 0, rtruncnorm(1, a=0, mean=1+1*0-2*x1[i]+0.8*x2[i]*x3[i]+0.5*x4[i]^2+1.5*x1[i]*x4[i], sd=0.2))
    }
    if(cluster.y[i]==2){
      a[i] <- rbinom(1,1,Expit(-0.1-0.3*x1[i]+0.3*x3[i]-0.5*x2[i]-0.6*x4[i]+0.6*x5[i]))
      p1[i] <- Expit(-3.5+0.5*1+0.1*x1[i]+0.5*x3[i]-0.5*x2[i]+0.3*x4[i]-0.1*x5[i])
      Pi1[i] <- sample(c(0,1), 1, prob=c(1-p1[i], p1[i]))
      p0[i] <- Expit(-3.5+0.5*0+0.1*x1[i]+0.5*x3[i]-0.5*x2[i]+0.3*x4[i]-0.1*x5[i])
      Pi0[i] <- sample(c(0,1), 1, prob=c(1-p0[i], p0[i]))
      Y1[i] <- ifelse(Pi1[i]==1, 0, rtruncnorm(1, a=0, mean=5-3*1+0.8*x1[i]-0.4*x2[i]+0.5*x3[i]+0.4*(x4[i])^2-2*x5[i], sd=0.2))
      Y0[i] <- ifelse(Pi0[i]==1, 0, rtruncnorm(1, a=0, mean=5-3*0+0.8*x1[i]-0.4*x2[i]+0.5*x3[i]+0.4*(x4[i])^2-2*x5[i], sd=0.2))
    }
    if(cluster.y[i]==3){
      a[i] <- rbinom(1,1,Expit(0.5-0.3*x1[i]+0.3*x3[i]-0.5*x2[i]-0.6*x4[i]+0.6*x5[i]))
      p1[i] <- Expit(-1.5-0.8*1+0.9*x1[i]-0.3*x3[i]+0.3*x2[i]-0.9*x4[i]-0.5*x5[i])
      Pi1[i] <- sample(c(0,1), 1, prob=c(1-p1[i], p1[i]))
      p0[i] <- Expit(-1.5-0.8*0+0.9*x1[i]-0.3*x3[i]+0.3*x2[i]-0.9*x4[i]-0.5*x5[i])
      Pi0[i] <- sample(c(0,1), 1, prob=c(1-p0[i], p0[i]))
      Y1[i] <- ifelse(Pi1[i]==1, 0, rtruncnorm(1, a=0, mean=2+1.8*1-0.5*x1[i]+0.5*x2[i]+0.5*x5[i]+0.8*(x4[i])-0.3*x4[i]*x3[i], sd=0.2))
      Y0[i] <- ifelse(Pi0[i]==1, 0, rtruncnorm(1, a=0, mean=2+1.8*0-0.5*x1[i]+0.5*x2[i]+0.5*x5[i]+0.8*(x4[i])-0.3*x4[i]*x3[i], sd=0.2))
    }
    if(cluster.y[i]==4){
      a[i] <- rbinom(1,1,Expit(-1.5-0.8*x1[i]+0.3*x3[i]+0.2*x2[i]-1*x4[i]+0.5*x5[i]))
      p1[i] <- Expit(-2-1*1+1*x1[i]-0.1*x3[i]+0.3*x2[i]-1*x4[i]-0.5*x5[i])
      Pi1[i] <- sample(c(0,1), 1, prob=c(1-p1[i], p1[i]))
      p0[i] <- Expit(-2-1*0+1*x1[i]-0.1*x3[i]+0.3*x2[i]-1*x4[i]-0.5*x5[i])
      Pi0[i] <- sample(c(0,1), 1, prob=c(1-p0[i], p0[i]))
      Y1[i] <- ifelse(Pi1[i]==1, 0, rtruncnorm(1,a=0,  mean=4.5+4.5*1+0.6*x1[i]-1.5*x2[i]+0.7*x3[i]*x4[i]+0.4*x5[i]*x3[i],sd=0.2))
      Y0[i] <- ifelse(Pi0[i]==1, 0, rtruncnorm(1, a=0, mean=4.5+4.5*0+0.6*x1[i]-1.5*x2[i]+0.7*x3[i]*x4[i]+0.4*x5[i]*x3[i],sd=0.2))
    }
  }
  
  
  Sy <- sample(c(1,2,3,4,5), 600, replace=TRUE)
  Sx <- sample(c(1,2),600, replace=TRUE)

  uniqueYX <- distinct(tbl_df(cbind(Sy,Sx)))
  uniqueYX <- uniqueYX[order(uniqueYX$Sy,uniqueYX$Sx),]

  n <- length(Sy)
  y <- a*Y1+(1-a)*Y0
  x <- cbind(1, x1, x2, x3, x4, x5)
  xa <- cbind(1, a, x1, x2, x3, x4, x5)
  
  yPAR <- list()
  yPAR_alt <- list()
  
  xPAR <- list()
  for(j in 1:5){
    yPAR[[j]] <- list(beta=matrix(rnorm(7,0,2),ncol=1), r=matrix(rnorm(7,0,2),ncol=1),  sigma=0.5)
  }
  for(j in 1:10){
    xPAR[[j]] <- list(p0=matrix(rep(0.5,2), nrow=1), p1=matrix(rep(0.5,2),nrow=1),  p2.mu=0, p2.sig=1,   p3.mu=0, p3.sig=1, p4.mu=0, p4.sig=1, p5.mu=0, p5.sig=1)
  }
  
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
  
  source("function.R")
  
  E <- list()
  E[[1]] <- list(EY=0)
  C <- list()
  C[[1]] <- list(Sx=Sx, Sy=Sy, xPAR=xPAR, yPAR=yPAR)
  P <- list()
  P[[1]] <- list(alpha_omega=alpha_omega,  alpha_theta=alpha_theta, xPAR=xPAR, yPAR=yPAR)
  
  E <- list()
  E[[1]] <- list(EY=0)
  
  E1 <- NULL
  E2 <- NULL
  E3 <- NULL
  E4 <- NULL

#------ Main MCMC run
  for(l in 2:4000){
    C[[l]] <- clustering(C[[l-1]]$Sx, C[[l-1]]$Sy, P[[l-1]]$xPAR, P[[l-1]]$yPAR, P[[l-1]]$alpha_omega, P[[l-1]]$alpha_theta)
    P[[l]] <- update_parameters(C[[l]]$Sx, C[[l]]$Sy, C[[l]]$xPAR, C[[l]]$yPAR)

#------ Estimating the effect
    EY11 <- EY00 <- NULL
    numY <- length(unique(C[[l]]$Sy))
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
  
    y.beta <- t(as.matrix(mnormt::rmnorm(1, mu.beta, diag(Sigma.beta,7))))
    y.r <- t(as.matrix(mnormt::rmnorm(1, mu.r, diag(Sigma.r,7))))
    y.sigma <- rgamma(1, 1, 5)
    
    a.p <- as.matrix(rdirichlet(200, alpha00))
    x1.p <- as.matrix(rdirichlet(200, alpha01))
    x2.var <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
    x2.mu <- rnorm(1, mean(x2), sqrt(x2.var/c0)) # c0: 0.5
    x3.var <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
    x3.mu <- rnorm(200, mean(x3), sqrt(x3.var/c0)) # c0: 0.5
    
    x4.var <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
    x4.mu <- rnorm(1, mean(x4), sqrt(x4.var/c0)) # c0: 0.5
    
    x5.var <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
    x5.mu <- rnorm(1, mean(x5), sqrt(x5.var/c0)) # c0: 0.5
    
    X0 <- cbind(1,1,x1,x2,x3,x4,x5)
    XX0 <- cbind(1,1,x1,x2,x3,x4,x5)
    
    prob <- matrix(nrow=n, ncol=length(nj_i)+1)
    prob[,1] <- P[[l]]$alpha_theta/(P[[l]]$alpha_theta+n)*mean(dcat(1+1, a.p))*
      sapply(1:n, function(x) mean(dcat(XX0[x,3]+1, x1.p)*
                                     dnorm(XX0[x,4], x2.mu, sqrt(x2.var))*
                                     dnorm(XX0[x,5], x3.mu, sqrt(x3.var))*
                                     dnorm(XX0[x,6], x4.mu, sqrt(x4.var))*
                                     dnorm(XX0[x,7], x5.mu, sqrt(x5.var))
      ))
    
    for(clus.x in 1:length(nj_i)){
      prob[,clus.x+1] <- nj_i[clus.x]/(P[[l]]$alpha_theta+n)*
        (P[[l]]$alpha_omega/(P[[l]]$alpha_omega+nj_i[clus.x])*mean(dcat(1+1, a.p))*sapply(1:n, function(x) mean(dcat(XX0[x,3]+1, x1.p)*
                                                                                                                  dnorm(XX0[x,4], x2.mu, sqrt(x2.var))*dnorm(XX0[x,5], x3.mu, sqrt(x3.var))*dnorm(XX0[x,6], x4.mu, sqrt(x4.var))*dnorm(XX0[x,7], x5.mu, sqrt(x5.var))))+
           rowSums(sapply(seq(Nlj[clus.x, 1],Nlj[clus.x, 2]) , function(x) nnlj_i[x]/(P[[l]]$alpha_omega+nj_i[clus.x])*dcat(1+1,P[[l]]$xPAR[[x]]$p0)*dnorm(XX0[,4], P[[l]]$xPAR[[x]]$p2.mu, sqrt(P[[l]]$xPAR[[x]]$p2.sig))*dnorm(XX0[,5], P[[l]]$xPAR[[x]]$p3.mu, sqrt(P[[l]]$xPAR[[x]]$p3.sig))*dnorm(XX0[,6], P[[l]]$xPAR[[x]]$p4.mu, sqrt(P[[l]]$xPAR[[x]]$p4.sig))*dnorm(XX0[,7], P[[l]]$xPAR[[x]]$p5.mu, sqrt(P[[l]]$xPAR[[x]]$p5.sig)) )))
    }
    
    EY11 <- NULL
    first.ind <- apply(prob, 1, function(x) rcat(1, x))
    zeros <- rbern(n, expit(t(y.r), XX0))
    for(temp.ind in 1:n){
      if(first.ind[temp.ind]==1){
        EY11[temp.ind] <- ifelse(zeros[temp.ind]==1, 0, rnorm(1, mean=XX0[temp.ind,]%*%t(y.beta), sd=y.sigma))
      }else{
        clus.x <- first.ind[temp.ind]-1
        zero <- rbern(1, expit(P[[l]]$yPAR[[clus.x]]$r, XX0[temp.ind,]))
        EY11[temp.ind] <- ifelse(zero==1, 0, rnorm(1, XX0[temp.ind,]%*%P[[l]]$yPAR[[clus.x]]$beta, P[[l]]$yPAR[[clus.x]]$sigma))
      }
    }
  
    X0 <- cbind(1,0,x1,x2,x3,x4,x5)
    XX0 <- cbind(1,0,x1,x2,x3,x4,x5)
    
    prob <- matrix(nrow=n, ncol=length(nj_i)+1)
    prob[,1] <- P[[l]]$alpha_theta/(P[[l]]$alpha_theta+n)*mean(dcat(1+0, a.p))*
      sapply(1:n, function(x) mean(dcat(XX0[x,3]+1, x1.p)*
                                     dnorm(XX0[x,4], x2.mu, sqrt(x2.var))*
                                     dnorm(XX0[x,5], x3.mu, sqrt(x3.var))*
                                     dnorm(XX0[x,6], x4.mu, sqrt(x4.var))*
                                     dnorm(XX0[x,7], x5.mu, sqrt(x5.var))
      ))
    
    for(clus.x in 1:length(nj_i)){
      prob[,clus.x+1] <- nj_i[clus.x]/(P[[l]]$alpha_theta+n)*
        (P[[l]]$alpha_omega/(P[[l]]$alpha_omega+nj_i[clus.x])*mean(dcat(1+0, a.p))*sapply(1:n, function(x) mean(dcat(XX0[x,3]+1, x1.p)*
                                                                                                                  dnorm(XX0[x,4], x2.mu, sqrt(x2.var))*dnorm(XX0[x,5], x3.mu, sqrt(x3.var))*dnorm(XX0[x,6], x4.mu, sqrt(x4.var))*dnorm(XX0[x,7], x5.mu, sqrt(x5.var))))+
           rowSums(sapply(seq(Nlj[clus.x, 1],Nlj[clus.x, 2]) , function(x) nnlj_i[x]/(P[[l]]$alpha_omega+nj_i[clus.x])*dcat(1+0,P[[l]]$xPAR[[x]]$p0)*dnorm(XX0[,4], P[[l]]$xPAR[[x]]$p2.mu, sqrt(P[[l]]$xPAR[[x]]$p2.sig))*dnorm(XX0[,5], P[[l]]$xPAR[[x]]$p3.mu, sqrt(P[[l]]$xPAR[[x]]$p3.sig))*dnorm(XX0[,6], P[[l]]$xPAR[[x]]$p4.mu, sqrt(P[[l]]$xPAR[[x]]$p4.sig))*dnorm(XX0[,7], P[[l]]$xPAR[[x]]$p5.mu, sqrt(P[[l]]$xPAR[[x]]$p5.sig)) )))
    }
    
    
    EY00 <- NULL
    first.ind <- apply(prob, 1, function(x) rcat(1, x))
    zeros <- rbern(n, expit(t(y.r), XX0))
    for(temp.ind in 1:n){
      if(first.ind[temp.ind]==1){
        EY00[temp.ind] <- ifelse(zeros[temp.ind]==1, 0, rnorm(1, mean=XX0[temp.ind,]%*%t(y.beta), sd=y.sigma))
      }else{
        clus.x <- first.ind[temp.ind]-1
        zero <- rbern(1, expit(P[[l]]$yPAR[[clus.x]]$r, XX0[temp.ind,]))
        EY00[temp.ind] <- ifelse(zero==1, 0, rnorm(1, XX0[temp.ind,]%*%P[[l]]$yPAR[[clus.x]]$beta, P[[l]]$yPAR[[clus.x]]$sigma))
      }
    }

    E[[l]] <- list(EY11=EY11,EY00=EY00, Ey = EY11*a+(1-a)*EY00)  # individual causal effects
    
    
#------ Find the posterior similarity matrix    
    pair <- cbind(C[[l-1]]$Sy,C[[l-1]]$Sx)
    colnames(pair) <- c("Sy", "Sx")
    unique.pair <- distinct(tbl_df(pair))
    bin <- matrix(NA,nrow=length(unlist(unique(unique.pair[,1]))), ncol=600)
    
    for(i in 1:length(unlist(unique(unique.pair[,1])))){
      ind1 <- which(pair[,1]==as.numeric(unique.pair[i,1]))
      for(j in 1:length(which(unique.pair[,1]==i))){
        ind2 <- which(pair[,1]==i & pair[,2]==j)
        bin[i,ind2] <- j
      }
    }
    BIN <- array(0,dim=c(600,600,dim(bin)[1]))
    for(t in 1:dim(bin)[1]){
      for(i in 1:600){
        if(is.na(bin[t,i])){
          BIN[i,,t] <- 0
          BIN[i,i,t] <- 1
        }else{
          BIN[i,which(bin[t,]==bin[t,i]),t] <- 1
          BIN[i,which(bin[t,]!=bin[t,i]),t] <- 0
        }
      }
      
    }
    temp.seq1 <- seq(2001,2500, by=1)
    temp.seq2 <- seq(2501,3000, by=1)
    temp.seq3 <- seq(3001,3500, by=1)
    temp.seq4 <- seq(3501,4000, by=1)

    if(l %in% temp.seq1){
      BIN.final1[[(l-2000)]] <- apply(BIN, 1:2, function(x) mean(x, na.rm=TRUE))
      if(l==max(temp.seq1)){
        BIN.final1 <- apply(simplify2array(BIN.final1), 1:2, function(x) mean(x, na.rm=TRUE))
      }
    }else{
      if(l %in% temp.seq2){
        BIN.final2[[(l-2500)]] <- apply(BIN, 1:2, function(x) mean(x, na.rm=TRUE))
        if(l==max(temp.seq2)){
          BIN.final2 <- apply(simplify2array(BIN.final2), 1:2, function(x) mean(x, na.rm=TRUE))
        }
      }else{
        if(l %in% temp.seq3){
          BIN.final3[[(l-3000)]] <- apply(BIN, 1:2, function(x) mean(x, na.rm=TRUE))
          if(l==max(temp.seq3)){
            BIN.final3 <- apply(simplify2array(BIN.final3), 1:2, function(x) mean(x, na.rm=TRUE))
          }
        }else{
          if(l %in% temp.seq4){
          BIN.final4[[(l-3500)]] <- apply(BIN, 1:2, function(x) mean(x, na.rm=TRUE))
          if(l==max(temp.seq4)){
            BIN.final4 <- apply(simplify2array(BIN.final4), 1:2, function(x) mean(x, na.rm=TRUE))
          }
        }
      }      
      }
    }
  }
  
  BIN.final <- list()
  BIN.final[[1]] <- BIN.final1
  BIN.final[[2]] <- BIN.final2
  BIN.final[[3]] <- BIN.final3
  BIN.final[[4]] <- BIN.final4
  BIN.final <- apply(simplify2array(BIN.final), 1:2, function(x) mean(x, na.rm=TRUE))
  
  
  BIN.final <- ifelse(is.na(BIN.final), 0, BIN.final)
  sim<-BIN.final
  
  best.sim<-NULL
  for(i in 1:9){
    best.sim[i] <- pam(1-sim, i, diss=TRUE)$silinfo$avg.width
  }
  num.clus <- which.max(best.sim)
  best <- pam(1-sim, num.clus, diss=TRUE)
  
  id <- best$clustering
  
  temp <- matrix(nrow=2000, ncol=length(unique(id)))
  
  for(t in 1:2000){
    for(tt in 1:length(unique(id))){
      temp[t, tt] <- mean(E[[t+2000]]$EY11[id==tt]-E[[t+2000]]$EY00[id==tt])
    }
  }
  t1 <- as.numeric(names(sort(table(id[1:150]),decreasing=TRUE)[1]))
  t2 <- as.numeric(names(sort(table(id[151:300]),decreasing=TRUE)[1]))
  t3 <- as.numeric(names(sort(table(id[301:450]),decreasing=TRUE)[1]))
  t4 <- as.numeric(names(sort(table(id[451:600]),decreasing=TRUE)[1]))
  
  temp <- temp[,c(t1,t2,t3,t4)]
#  colMeans(temp)
  


