#------ Random number function and for the (scaled) inverse-chi-squared distribution.
rinvchisq <- function (n, df, scale = 1/df){
  df <- rep(df, len = n)
  scale <- rep(scale, len = n)
  if (any(df <= 0)) 
    stop("The df parameter must be positive.")
  if (any(scale <= 0)) 
    stop("The scale parameter must be positive.")
  x <- (df * scale)/rchisq(n, df = df)
  return(x)
}

#------ Expit function
expit <- function(beta,x){
  p <- exp(x%*%beta)/(1+exp(x%*%beta))
  return(p)
}

#------ Function to (sub)clustering Y, A, and X
clustering <- function(Sx, Sy, xPAR, yPAR,  alpha_omega,  alpha_theta){
  
  # Loop through every person and change cluster memberships
  for(i in 1:n){ # the number of subjects : n
    
    # original copies for each iteration
    SSy <- Sy; SSx <- Sx; 
    
    uniqueYX <- distinct(tbl_df(cbind(Sy,Sx)))
    uniqueYX <- uniqueYX[order(uniqueYX$Sy,uniqueYX$Sx),]
    
    lenY <- length(which(Sy==Sy[i]))-1
    lenX <- length(which(Sx==Sx[i] & Sy==Sy[i]))-1
    indY <- which(order(unique(Sy))==Sy[i])
    indYX <- which(uniqueYX[,1]==Sy[i])
    
    indX <- which(uniqueYX[,1]==Sy[i] & uniqueYX[,2]==Sx[i])
    
    indS1 <- which(uniqueYX[,1]==Sy[i])
    indS11 <- which(Sy==Sy[i] & Sx==Sx[i])
    
    # temp. storages of parameters
    yPAR.temp <- yPAR;
    xPAR.temp <- xPAR;
    
    if(length(which(Sy == Sy[i] & Sx == Sx[i]))==1 ){
      
      if(length(which(Sy == Sy[i]))==1){
        yPAR.temp <- yPAR.temp[-Sy[i]] #yPAR has a list format
        xPAR.temp <- xPAR.temp[-indX]
      } else {
        if(length(indS11)==1){
          # yPAR (coefficients and sigmas) has a list format
          xPAR.temp <- xPAR.temp[-indX] 
          
          # relabel Y cluster
          indYY <- which(Sy==Sy[i])
          indYX <- which(indYY == i)
          Sx[indYY[-indYX]] <- as.numeric(factor(Sx[indYY[-indYX]], labels=1:length(unique(Sx[indYY[-indYX]])), levels=sort(unique(Sx[indYY[-indYX]]))))
          indS <- which(uniqueYX[,1]==Sy[i])
          uniqueYX[indS[-Sx[i]],2] <- as.numeric(factor(uniqueYX[indS[-Sx[i]],2]$Sx, labels=1:length(unique(uniqueYX[indS[-Sx[i]],2]$Sx)), levels=sort(unique(uniqueYX[indS[-Sx[i]],2]$Sx))))
        }
        
      }
      
      # relabel Y cluster (if needed)
      if(length(which(Sy == Sy[i]))==1 & length(which(Sy > Sy[i])) > 0){
        indYY <- which(Sy>Sy[i])
        Sy[indYY] <- Sy[indYY]-1 
        indS <- which(uniqueYX[,1] > Sy[i])
        uniqueYX[indS,1] <- uniqueYX[indS,1]-1
      }
      # Get rid of the row for i-th individual (uniqueYX & uniqueXA) # may not need this step
      if(length(indS11)==1){      
        uniqueYX <- uniqueYX[-indX,]
      }
      
    } 
    
    # temp variables for Sy[i], Sx[i]
    St <- c(Sy[i], Sx[i])
    
    # delete row of Sy, Sx
    Sx <- Sx[-i]; Sy <- Sy[-i]
    
    # recalculate number of unique clusters
    numY <- length(yPAR.temp); 
    numX <- length(xPAR.temp); # <- same as numXj?
    
    # counts for # in X,Y cluster excluding the ith person
    nj_i <- table(factor(Sy, levels=1:numY));
    nlj_i <- table(Sy, Sx);
    
    # of X clusters within each Y cluster
    numXj <- apply(nlj_i, 1, function(x) length(which(x!=0)))
  
    # The current value of S[i] is in one of the remaining clusters
    if(lenY*lenX != 0){
      y.coef.r <- matrix(mnormt::rmnorm(M, rep(mu.r,7), diag(Sigma.r,7)), ncol=M, byrow=TRUE) # mu.r=rep(0,7), Sigma.r=100*diag(0,7)
      y.coef.beta <- matrix(mnormt::rmnorm(M, rep(mu.beta,7), diag(Sigma.beta,7)), ncol=M, byrow=TRUE) # mu.beta=rep(0,7), Sigma.beta=100*diag(0,7)
      y.coef.sigma <- rgamma(M, 1, 5)
      
      for(m in 1:M){
        x.cat0 <- MCMCpack::rdirichlet(1, alpha00)
        x.cat1 <- MCMCpack::rdirichlet(1, alpha01) # p11: # of categorical X1, alpha01: initial hyperparameters (vector) 
        x.var2 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu2 <- rnorm(1, 0, sqrt(x.var2/c0)) # c0: 0.5
        x.var3 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu3 <- rnorm(1, 0, sqrt(x.var3/c0)) # c0: 0.5
        x.var4 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu4 <- rnorm(1, 0, sqrt(x.var4/c0)) # c0: 0.5
        x.var5 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu5 <- rnorm(1, 0, sqrt(x.var5/c0)) # c0: 0.5
        x.var6 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu6 <- rnorm(1, 0, sqrt(x.var6/c0)) # c0: 0.5
        x.var7 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu7 <- rnorm(1, 0, sqrt(x.var7/c0)) # c0: 0.5
        xPAR.temp[[(numX+m)]] <- list(p0=x.cat0, p1=x.cat1, p2.mu=x.mu2, p2.sig=x.var2, p3.mu=x.mu3,  p3.sig=x.var3, p4.mu=x.mu4, p4.sig=x.var4, p5.mu=x.mu5, p5.sig=x.var5, p6.mu=x.mu6, p6.sig=x.var6, p7.mu=x.mu7, p7.sig=x.var7)
        yPAR.temp[[(numY+m)]] <- list(beta=matrix(y.coef.beta[,m], ncol=1), r=matrix(y.coef.r[,m], ncol=1),sigma=y.coef.sigma[m])
      }
      
      for(m in 1:M){
        for(k in 1:numY){
          
          x.cat0 <- MCMCpack::rdirichlet(1, alpha00)
          x.cat1 <- MCMCpack::rdirichlet(1, alpha01) # p11: # of categorical X1, alpha01: initial hyperparameters (vector) 
          x.var2 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu2 <- rnorm(1, 0, sqrt(x.var2/c0)) # c0: 0.5
          x.var3 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu3 <- rnorm(1, 0, sqrt(x.var3/c0)) # c0: 0.5
          x.var4 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu4 <- rnorm(1, 0, sqrt(x.var4/c0)) # c0: 0.5
          x.var5 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu5 <- rnorm(1, 0, sqrt(x.var5/c0)) # c0: 0.5
          x.var6 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu6 <- rnorm(1, 0, sqrt(x.var6/c0)) # c0: 0.5
          x.var7 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu7 <- rnorm(1, 0, sqrt(x.var7/c0)) # c0: 0.5
          xPAR.temp <- rlist::list.insert(xPAR.temp, cumsum(numXj+m)[k], list(p0=x.cat0, p1=x.cat1, p2.mu=x.mu2, p2.sig=x.var2, p3.mu=x.mu3,  p3.sig=x.var3, p4.mu=x.mu4, p4.sig=x.var4, p5.mu=x.mu5, p5.sig=x.var5, p6.mu=x.mu6, p6.sig=x.var6, p7.mu=x.mu7, p7.sig=x.var7) )
        }
      }
    }
    
    # The current value of S[i] is in an existing x-cluster, but not in an existing y-subcluster, and in an existing a-subcluster
    if(lenY != 0 & lenX == 0){
      
      y.coef.r <- matrix(mnormt::rmnorm(M, rep(mu.r,7), diag(Sigma.r,7)), ncol=M, byrow=TRUE) # mu.r=rep(0,7), Sigma.r=100*diag(0,7)
      y.coef.beta <- matrix(mnormt::rmnorm(M, rep(mu.beta,7), diag(Sigma.beta,7)), ncol=M, byrow=TRUE) # mu.beta=rep(0,7), Sigma.beta=100*diag(0,7)
      y.coef.sigma <- rgamma(M, 1, 5)
      
      for(m in 1:M){
        x.cat0 <- MCMCpack::rdirichlet(1, alpha00)
        x.cat1 <- MCMCpack::rdirichlet(1, alpha01) # p11: # of categorical X1, alpha01: initial hyperparameters (vector) 
        x.var2 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu2 <- rnorm(1, 0, sqrt(x.var2/c0)) # c0: 0.5
        x.var3 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu3 <- rnorm(1, 0, sqrt(x.var3/c0)) # c0: 0.5
        x.var4 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu4 <- rnorm(1, 0, sqrt(x.var4/c0)) # c0: 0.5
        x.var5 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu5 <- rnorm(1, 0, sqrt(x.var5/c0)) # c0: 0.5
        x.var6 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu6 <- rnorm(1, 0, sqrt(x.var6/c0)) # c0: 0.5
        x.var7 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu7 <- rnorm(1, 0, sqrt(x.var7/c0)) # c0: 0.5
        xPAR.temp[[(numX+m)]] <- list(p0=x.cat0, p1=x.cat1, p2.mu=x.mu2, p2.sig=x.var2, p3.mu=x.mu3,  p3.sig=x.var3, p4.mu=x.mu4, p4.sig=x.var4, p5.mu=x.mu5, p5.sig=x.var5, p6.mu=x.mu6, p6.sig=x.var6, p7.mu=x.mu7, p7.sig=x.var7)
        yPAR.temp[[(numY+m)]] <- list(beta=matrix(y.coef.beta[,m], ncol=1),  r=matrix(y.coef.r[,m], ncol=1),sigma=y.coef.sigma[m])
      }
      for(m in 1:M){ 
        for(k in 1:numY){
          x.cat0 <- MCMCpack::rdirichlet(1, alpha00)
          x.cat1 <- MCMCpack::rdirichlet(1, alpha01) # p11: # of categorical X1, alpha01: initial hyperparameters (vector) 
          x.var2 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu2 <- rnorm(1, 0, sqrt(x.var2/c0)) # c0: 0.5
          x.var3 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu3 <- rnorm(1, 0, sqrt(x.var3/c0)) # c0: 0.5
          x.var4 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu4 <- rnorm(1, 0, sqrt(x.var4/c0)) # c0: 0.5
          x.var5 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu5 <- rnorm(1, 0, sqrt(x.var5/c0)) # c0: 0.5
          x.var6 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu6 <- rnorm(1, 0, sqrt(x.var6/c0)) # c0: 0.5
          x.var7 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu7 <- rnorm(1, 0, sqrt(x.var7/c0)) # c0: 0.5
          xPAR.temp <- rlist::list.insert(xPAR.temp, cumsum(numXj+m)[k], list(p0=x.cat0, p1=x.cat1, p2.mu=x.mu2, p2.sig=x.var2, p3.mu=x.mu3,  p3.sig=x.var3, p4.mu=x.mu4, p4.sig=x.var4, p5.mu=x.mu5, p5.sig=x.var5, p6.mu=x.mu6, p6.sig=x.var6, p7.mu=x.mu7, p7.sig=x.var7) )
        }
      }
      xPAR.temp <- xPAR.temp[-(cumsum(numXj+M)[St[1]]-(M-1))]
      xPAR.temp <- rlist::list.insert(xPAR.temp, (cumsum(numXj+M)[St[1]]-(M-1)), xPAR[[indX]] )
    }
    
    # The current value of S[i] is not in an existing x-cluster
    if(lenY == 0){
      y.coef.r <- matrix(mnormt::rmnorm(M, rep(mu.r,7), diag(Sigma.r,7)), ncol=M, byrow=TRUE) # mu.r=rep(0,7), Sigma.r=100*diag(0,7)
      y.coef.beta <- matrix(mnormt::rmnorm(M, rep(mu.beta,7), diag(Sigma.beta,7)), ncol=M, byrow=TRUE) # mu.beta=rep(0,7), Sigma.beta=100*diag(0,7)
      y.coef.sigma <- rgamma(M, 1, 5)
      
      for(m in 1:M){
        x.cat0 <- MCMCpack::rdirichlet(1, alpha00)
        x.cat1 <- MCMCpack::rdirichlet(1, alpha01) # p11: # of categorical X1, alpha01: initial hyperparameters (vector) 
        x.var2 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu2 <- rnorm(1, 0, sqrt(x.var2/c0)) # c0: 0.5
        x.var3 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu3 <- rnorm(1, 0, sqrt(x.var3/c0)) # c0: 0.5
        x.var4 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu4 <- rnorm(1, 0, sqrt(x.var4/c0)) # c0: 0.5
        x.var5 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu5 <- rnorm(1, 0, sqrt(x.var5/c0)) # c0: 0.5
        x.var6 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu6 <- rnorm(1, 0, sqrt(x.var6/c0)) # c0: 0.5
        x.var7 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
        x.mu7 <- rnorm(1, 0, sqrt(x.var7/c0)) # c0: 0.5
        xPAR.temp[[(numX+m)]] <- list(p0=x.cat0, p1=x.cat1, p2.mu=x.mu2, p2.sig=x.var2, p3.mu=x.mu3,  p3.sig=x.var3, p4.mu=x.mu4, p4.sig=x.var4, p5.mu=x.mu5, p5.sig=x.var5, p6.mu=x.mu6, p6.sig=x.var6, p7.mu=x.mu7, p7.sig=x.var7)
        yPAR.temp[[(numY+m)]] <- list(beta=matrix(y.coef.beta[,m], ncol=1),  r=matrix(y.coef.r[,m], ncol=1),sigma=y.coef.sigma[m])
      }  
      for(m in 1:M){  
        for(k in 1:numY){
          x.cat0 <- MCMCpack::rdirichlet(1, alpha00)
          x.cat1 <- MCMCpack::rdirichlet(1, alpha01) # p11: # of categorical X1, alpha01: initial hyperparameters (vector) 
          x.var2 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu2 <- rnorm(1, 0, sqrt(x.var2/c0)) # c0: 0.5
          x.var3 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu3 <- rnorm(1, 0, sqrt(x.var3/c0)) # c0: 0.5
          x.var4 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu4 <- rnorm(1, 0, sqrt(x.var4/c0)) # c0: 0.5
          x.var5 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu5 <- rnorm(1, 0, sqrt(x.var5/c0)) # c0: 0.5
          x.var6 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu6 <- rnorm(1, 0, sqrt(x.var6/c0)) # c0: 0.5
          x.var7 <- rinvchisq(1, nu0, tau0) # p2: # of continuous X, nu0: dfs (2), tau0: scale (2)
          x.mu7 <- rnorm(1, 0, sqrt(x.var7/c0)) # c0: 0.5
          xPAR.temp <- rlist::list.insert(xPAR.temp, cumsum(numXj+m)[k], list(p0=x.cat0, p1=x.cat1, p2.mu=x.mu2, p2.sig=x.var2, p3.mu=x.mu3,  p3.sig=x.var3, p4.mu=x.mu4, p4.sig=x.var4, p5.mu=x.mu5, p5.sig=x.var5, p6.mu=x.mu6, p6.sig=x.var6, p7.mu=x.mu7, p7.sig=x.var7) )
        }
      }
      yPAR.temp <- yPAR.temp[-(numY+1)]
      yPAR.temp <- rlist::list.insert(yPAR.temp, (numY+1), yPAR[[indY]] )  
      xPAR.temp <- xPAR.temp[-(cumsum(numXj+M)[numY]+1)]
      xPAR.temp <- rlist::list.insert(xPAR.temp, (cumsum(numXj+M)[numY]+1), xPAR[[indX]] )      
    }
    p <- list()
    xx <- cbind(1, a, x1, x2, x3, x4, x5, x6, x7)
    xxx <- cbind(1, a, x1, x2, x3, x4, x5)
    
    for(j in 1:(numY)){
      p[[j]] <- list(ya=matrix(nrow=numXj[j]+M,ncol=1))
      if(j==1){indBaseX <- 0}else{indBaseX <- cumsum(numXj+M)[j-1]}
      
      for(l in 1:(numXj[j]+M)){
        if(l <= numXj[j]){
          p[[j]]$ya[l,1] <- nj_i[j]*nlj_i[j,l]/(nj_i[j]+alpha_omega)*
            dcat(a[i]+1, xPAR.temp[[indBaseX+l]]$p0)*
            dcat(x1[i]+1, xPAR.temp[[indBaseX+l]]$p1)*
            dnorm(x2[i],xPAR.temp[[indBaseX+l]]$p2.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p2.sig[1]))*
            dnorm(x3[i],xPAR.temp[[indBaseX+l]]$p3.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p3.sig[1]))*
            dnorm(x4[i],xPAR.temp[[indBaseX+l]]$p4.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p4.sig[1]))*
            dnorm(x5[i],xPAR.temp[[indBaseX+l]]$p5.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p5.sig[1]))*
            dnorm(x6[i],xPAR.temp[[indBaseX+l]]$p6.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p6.sig[1]))*
            dnorm(x7[i],xPAR.temp[[indBaseX+l]]$p7.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p7.sig[1]))*
            (I(y[i]==0)*expit(yPAR.temp[[j]]$r,xxx[i,])+(1-expit(yPAR.temp[[j]]$r,xxx[i,]))*
               dnorm(y[i], mean=xxx[i,]%*%(yPAR.temp[[j]]$beta), sd=yPAR.temp[[j]]$sigma) )
        } else{
          p[[j]]$ya[l,1] <- nj_i[j]*(alpha_omega/M)/(nj_i[j]+alpha_omega)*
            dcat(a[i]+1, xPAR.temp[[indBaseX+l]]$p0)*
            dcat(x1[i]+1, xPAR.temp[[indBaseX+l]]$p1)*
            dnorm(x2[i],xPAR.temp[[indBaseX+l]]$p2.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p2.sig[1]))*
            dnorm(x3[i],xPAR.temp[[indBaseX+l]]$p3.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p3.sig[1]))*
            dnorm(x4[i],xPAR.temp[[indBaseX+l]]$p4.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p4.sig[1]))*
            dnorm(x5[i],xPAR.temp[[indBaseX+l]]$p5.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p5.sig[1]))*
            dnorm(x6[i],xPAR.temp[[indBaseX+l]]$p6.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p6.sig[1]))*
            dnorm(x7[i],xPAR.temp[[indBaseX+l]]$p7.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p7.sig[1]))*
            (I(y[i]==0)*expit(yPAR.temp[[j]]$r,xxx[i,])+(1-expit(yPAR.temp[[j]]$r,xxx[i,]))*
               dnorm(y[i], mean=xxx[i,]%*%(yPAR.temp[[j]]$beta), sd=yPAR.temp[[j]]$sigma) )
        }
      }
    }
    
    for(j in (numY+1):(numY+M)){
      jj <- j - numY
      p[[j]] <- list(ya=1)
      p[[j]]$ya <- (alpha_theta/M)*dcat(a[i]+1, xPAR.temp[[j]]$p0)*
        dcat(a[i]+1, xPAR.temp[[indBaseX+l]]$p0)*
        dcat(x1[i]+1, xPAR.temp[[indBaseX+l]]$p1)*
        dnorm(x2[i],xPAR.temp[[indBaseX+l]]$p2.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p2.sig[1]))*
        dnorm(x3[i],xPAR.temp[[indBaseX+l]]$p3.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p3.sig[1]))*
        dnorm(x4[i],xPAR.temp[[indBaseX+l]]$p4.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p4.sig[1]))*
        dnorm(x5[i],xPAR.temp[[indBaseX+l]]$p5.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p5.sig[1]))*
        dnorm(x6[i],xPAR.temp[[indBaseX+l]]$p6.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p6.sig[1]))*
        dnorm(x7[i],xPAR.temp[[indBaseX+l]]$p7.mu[1], sqrt(xPAR.temp[[indBaseX+l]]$p7.sig[1]))*
        (I(y[i]==0)*expit(yPAR.temp[[j]]$r,xxx[i,])+(1-expit(yPAR.temp[[j]]$r,xxx[i,]))*
           dnorm(y[i], mean=xxx[i,]%*%(yPAR.temp[[j]]$beta), sd=yPAR.temp[[j]]$sigma) )
    }

    p.len <- length(p)
    p.lenlen <- sapply(1:p.len, function(xx) length(p[[xx]]$ya))
    
    p.Sy <- unlist(sapply(1:p.len, function(xx) rep(xx,p.lenlen[xx])))
    p.Sx <- unlist(sapply(1:p.len, function(xx) rep(1:p.lenlen[xx], 1)))
    
    sSyx <- sample(1:sum(p.lenlen), 1, replace=FALSE, prob=unlist(p))
    s.Sy <- p.Sy[sSyx]
    s.Sx <- p.Sx[sSyx]
    
    if(i==1){Sy <- c(s.Sy, Sy)}else{Sy <- c(Sy[1:(i-1)], s.Sy, Sy[-(1:(i-1))])}
    if(s.Sy <= numY){yPAR <- yPAR.temp[-(numY+M)]}else{yPAR <- yPAR.temp}
    if(i==1){Sx <- c(s.Sx, Sx)}else{Sx <- c(Sx[1:(i-1)],  s.Sx, Sx[-(1:(i-1))])}
    if(s.Sy <= numY){
      if(s.Sx <= numXj[s.Sy]){
        xPAR <- xPAR.temp[-c(cumsum(numXj+M),cumsum(numXj+M)[numY]+M)]
      }else{
        xPAR <- xPAR.temp[-c(cumsum(numXj+M)[-s.Sy],cumsum(numXj+M)[numY]+M)]
      }
    }else{
      xPAR <- xPAR.temp[-c(cumsum(numXj+M))];
    }
    
    # recalculate number of unique clusters
    numY <- length(yPAR); 
    numX <- length(xPAR); # <- same as numYj?
    
    # counts for # in X,Y,A cluster excluding the ith person
    nj_i <- table(factor(Sy, levels=1:numY));
    nlj_i <- table(Sy, Sx);
    
    # of Y or A clusters within each X cluster
    numXj <- apply(nlj_i, 1, function(x) length(which(x!=0)))
    
    uniqueYX <- distinct(tbl_df(cbind(Sy,Sx)))
    uniqueYX <- uniqueYX[order(uniqueYX$Sy,uniqueYX$Sx),]
  }
  return(list(Sx=Sx, Sy=Sy, xPAR=xPAR, yPAR=yPAR))
}



update_parameters <- function(Sx, Sy, xPAR, yPAR){
  
  # recalculate number of unique clusters
  numY <- length(yPAR); 
  numX <- length(xPAR); # <- same as numYj?
  
  uniqueYX <- distinct(tbl_df(cbind(Sy,Sx)))
  uniqueYX <- uniqueYX[order(uniqueYX$Sy,uniqueYX$Sx),]
  
  uniqueY <- unique(Sy)
  
  # counts for # in X,Y,A cluster excluding the ith person
  nj_i <- table(factor(Sy, levels=1:numY));
  nlj_i <- table(Sy, Sx);

  mu0 <- 0
  for(j in 1:length(uniqueY)){
    
    X_lj <- matrix(subset(cbind(1,a,x1, x2, x3, x4, x5), Sy==j), ncol=7)
    prop <- rgamma(1, 10000*yPAR[[j]]$sigma, 10000)
    num <- sum(log(expit(yPAR[[j]]$r, X_lj)*identity(subset(y, Sy==j)==0)+(1-expit(yPAR[[j]]$r, X_lj))*dnorm(subset(y, Sy==j),  mean=X_lj%*%yPAR[[j]]$beta, sd=prop) ))+dgamma(prop, 1,5,log=TRUE) + dgamma(yPAR[[j]]$sigma,10000*prop,10000,log=TRUE)
    den <- sum(log(expit(yPAR[[j]]$r, X_lj)*identity(subset(y, Sy==j)==0)+(1-expit(yPAR[[j]]$r, X_lj))*dnorm(subset(y, Sy==j), mean=X_lj%*%yPAR[[j]]$beta, sd=yPAR[[j]]$sigma) ))+dgamma(yPAR[[j]]$sigma, 1,5,log=TRUE) + dgamma(prop,10000*yPAR[[j]]$sigma,10000,log=TRUE)
    if (log(runif(1))< (num-den)){
      yPAR[[j]]$sigma <- prop
    }
    prop <-  yPAR[[j]]$beta
    for(cc in 1:length(yPAR[[j]]$beta)){   
      prop[cc] <- yPAR[[j]]$beta[cc] + rnorm(1, 0, 0.05)
      num <- sum(log(expit(yPAR[[j]]$r, X_lj)*identity(subset(y, Sy==j)==0)+(1-expit(yPAR[[j]]$r, X_lj))*dnorm(subset(y, Sy==j), mean=X_lj%*%prop, sd=yPAR[[j]]$sigma) ))+dnorm(prop[cc], mu.beta, Sigma.beta, log=TRUE)
      den <- sum(log(expit(yPAR[[j]]$r, X_lj)*identity(subset(y, Sy==j)==0)+(1-expit(yPAR[[j]]$r, X_lj))*dnorm(subset(y, Sy==j), mean=X_lj%*%yPAR[[j]]$beta, sd=yPAR[[j]]$sigma) ))+dnorm(yPAR[[j]]$beta[cc], mu.beta, Sigma.beta, log=TRUE)
      if (log(runif(1))< (num-den)){
        yPAR[[j]]$beta[cc] <- prop[cc]
      }else{ prop[cc] <- yPAR[[j]]$beta[cc]}
    }    
    prop <-  yPAR[[j]]$r
    for(cc in 1:length(yPAR[[j]]$r)){ 
      prop[cc] <- yPAR[[j]]$r[cc] + rnorm(1, 0, 0.05)
      num <- sum(log(expit(prop, X_lj)*identity(subset(y, Sy==j)==0)+(1-expit(prop, X_lj))*dnorm(subset(y, Sy==j),mean=X_lj%*%yPAR[[j]]$beta, sd=yPAR[[j]]$sigma)   ))+dnorm(prop[cc], mu.r, Sigma.r, log=TRUE)
      den <- sum(log(expit(yPAR[[j]]$r, X_lj)*identity(subset(y, Sy==j)==0)+(1-expit(yPAR[[j]]$r, X_lj))*dnorm(subset(y, Sy==j), mean=X_lj%*%yPAR[[j]]$beta, sd=yPAR[[j]]$sigma)   ))+dnorm(yPAR[[j]]$r[cc], mu.r, Sigma.r, log=TRUE)
      if (log(runif(1))< (num-den)){
        yPAR[[j]]$r[cc] <- prop[cc]
      }else{ prop[cc] <- yPAR[[j]]$r[cc]}
    }  

    for(l in 1:length(which(uniqueYX$Sy==j))){
      ind <- which(uniqueYX$Sy==j & uniqueYX$Sx==l)
      xPAR[[ind]]$p0 <- MCMCpack::rdirichlet(1, alpha00+c(table(factor(subset(a+1, Sy==j & Sx==l), levels=1:length(alpha00)))))
      xPAR[[ind]]$p1 <- MCMCpack::rdirichlet(1, alpha01+c(table(factor(subset(x1+1, Sy==j & Sx==l), levels=1:length(alpha01)))))
      sd2 <- ifelse(is.na(sd(subset(x2, Sy==j & Sx==l))^2), 0, sd(subset(x2, Sy==j & Sx==l))^2)
      xPAR[[ind]]$p2.sig <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x2, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR[[ind]]$p2.sig
      xPAR[[ind]]$p2.mu <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x2, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
      
      sd2 <- ifelse(is.na(sd(subset(x3, Sy==j & Sx==l))^2), 0, sd(subset(x3, Sy==j & Sx==l))^2)
      xPAR[[ind]]$p3.sig <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x3, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR[[ind]]$p3.sig
      xPAR[[ind]]$p3.mu <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x3, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
      
      sd2 <- ifelse(is.na(sd(subset(x4, Sy==j & Sx==l))^2), 0, sd(subset(x4, Sy==j & Sx==l))^2)
      xPAR[[ind]]$p4.sig <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x4, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR[[ind]]$p4.sig
      xPAR[[ind]]$p4.mu <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x4, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
      
      sd2 <- ifelse(is.na(sd(subset(x5, Sy==j & Sx==l))^2), 0, sd(subset(x5, Sy==j & Sx==l))^2)
      xPAR[[ind]]$p5.sig <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x5, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR[[ind]]$p5.sig
      xPAR[[ind]]$p5.mu <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x5, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
      
      sd2 <- ifelse(is.na(sd(subset(x6, Sy==j & Sx==l))^2), 0, sd(subset(x6, Sy==j & Sx==l))^2)
      xPAR[[ind]]$p6.sig <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x6, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR[[ind]]$p6.sig
      xPAR[[ind]]$p6.mu <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x6, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
      
      sd2 <- ifelse(is.na(sd(subset(x7, Sy==j & Sx==l))^2), 0, sd(subset(x7, Sy==j & Sx==l))^2)
      xPAR[[ind]]$p7.sig <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x7, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR[[ind]]$p7.sig
      xPAR[[ind]]$p7.mu <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x7, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
    }
  }
  
  # Gamma prior
  eta <- rbeta(1, alpha_theta+1, n)
  pi <- (1+length(uniqueY)-1)/(n*(1-log(eta))) / (1 + (1+length(uniqueY)-1)/(n*(1-log(eta))))
  pi.temp <- sample(c(1,0), 1,prob=c(pi,1-pi))
  alpha_theta <- pi.temp*rgamma(1,1+length(uniqueY), 1-log(eta))+(1-pi.temp)*rgamma(1,1+length(uniqueY)-1, 1-log(eta))
  
  # counts for # in X,Y,A cluster excluding the ith person
  nj_i <- table(factor(Sy, levels=1:length(uniqueY)));
  nlj_i <- table(Sy, Sx);
  
  # of Y or A clusters within each X cluster
  numXj <- apply(nlj_i, 1, function(x) length(which(x!=0)))

  prop <- rgamma(1,alpha_omega*1000, 1000)
  rat <- log(dgamma(prop, 1, 1)*prop^(sum(numXj-1))*prod((prop+nj_i)*beta(prop+1,nj_i)))+log(dgamma(alpha_omega,prop*1000,1000))-log(dgamma(alpha_omega, 1, 1)*alpha_omega^(sum(numXj-1))*prod((alpha_omega+nj_i)*beta(alpha_omega+1,nj_i)))-log(dgamma(prop,alpha_omega*1000,1000))
  
  if(log(runif(1)) > rat) {
    prop <- alpha_omega
  }else{alpha_omega <- prop}
  
  return(list(alpha_omega=alpha_omega,  alpha_theta=alpha_theta, xPAR=xPAR, yPAR=yPAR))
}
