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

#------ Update parameters based on full conditionals
update_parameters <- function(Sx, Sy, xPAR_p0, xPAR_p1, xPAR_mu, xPAR_sig, yPAR_beta, yPAR_r, yPAR_sig){
  
  # recalculate number of unique clusters
  numY <- dim(yPAR_beta)[2]; 
  numX <- dim(xPAR_mu)[2];
  
  sysx <- cbind(Sy,Sx)
  colnames(sysx) <- c("Sy", "Sx")
  uniqueYX <- distinct(tbl_df(sysx))
  uniqueYX <- uniqueYX[order(uniqueYX$Sy,uniqueYX$Sx),]
  
  uniqueY <- unique(Sy)
  
  # counts for # in X,Y,A cluster excluding the ith person
  nj_i <- table(factor(Sy, levels=1:numY));
  nlj_i <- table(Sy, Sx);
  
  mu0 <- 0
  for(j in 1:length(uniqueY)){
    
    X_lj <- matrix(subset(xa, Sy==j), ncol=7)
    prop <- rgamma(1, 10000*yPAR_sig[j,1], 10000)
    num <- sum(log(expit(yPAR_r[,j], X_lj)*identity(subset(y, Sy==j)==0)+identity(subset(y, Sy==j)!=0)*(1-expit(yPAR_r[,j], X_lj))*dtruncnorm(subset(y, Sy==j), a=0, b=Inf, mean=X_lj%*%yPAR_beta[,j], sd=sqrt(prop)) ))+invgamma::dinvgamma(prop, 10,10000,log=TRUE) + dgamma(yPAR_sig[j,1],10000*prop,10000,log=TRUE)
    den <- sum(log(expit(yPAR_r[,j], X_lj)*identity(subset(y, Sy==j)==0)+identity(subset(y, Sy==j)!=0)*(1-expit(yPAR_r[,j], X_lj))*dtruncnorm(subset(y, Sy==j), a=0, b=Inf, mean=X_lj%*%yPAR_beta[,j], sd=sqrt(yPAR_sig[j,1])) ))+invgamma::dinvgamma(yPAR_sig[j,1], 10,10000,log=TRUE) + dgamma(prop,10000*yPAR_sig[j,1],10000,log=TRUE)
    if (log(runif(1))< (num-den) & !is.nan(num-den)){
      yPAR_sig[j,1] <- prop
    }
    
    prop <-  yPAR_beta[,j]
    for(cc in 1:length(yPAR_beta[,j])){   
      prop[cc] <- yPAR_beta[cc,j] + rnorm(1, 0, 0.05)
      num <- sum(log(expit(yPAR_r[,j], X_lj)*identity(subset(y, Sy==j)==0)+identity(subset(y, Sy==j)!=0)*(1-expit(yPAR_r[,j], X_lj))*dtruncnorm(subset(y, Sy==j),  a=0, b=Inf, mean=X_lj%*%prop, sd=sqrt(yPAR_sig[j,1])) ))+dnorm(prop[cc], mu.beta, sqrt(Sigma.beta), log=TRUE)
      den <- sum(log(expit(yPAR_r[,j], X_lj)*identity(subset(y, Sy==j)==0)+identity(subset(y, Sy==j)!=0)*(1-expit(yPAR_r[,j], X_lj))*dtruncnorm(subset(y, Sy==j),  a=0, b=Inf, mean=X_lj%*%yPAR_beta[,j], sd=sqrt(yPAR_sig[j,1])) ))+dnorm(yPAR_beta[cc,j], mu.beta, sqrt(Sigma.beta), log=TRUE)
      if (log(runif(1))< (num-den) & !is.nan(num-den)){
        yPAR_beta[cc,j] <- prop[cc]
      
      }else{ 
        prop[cc] <- yPAR_beta[cc,j]
        }
    }
    
    prop <-  yPAR_r[,j]
    for(cc in 1:length(yPAR_r[,j])){ 
      prop[cc] <- yPAR_r[cc,j] + rnorm(1, 0, 0.05)
      num <- sum(log(expit(prop, X_lj)*identity(subset(y, Sy==j)==0)+identity(subset(y, Sy==j)!=0)*(1-expit(prop, X_lj))*dtruncnorm(subset(y, Sy==j), a=0, b=Inf, mean=X_lj%*%yPAR_beta[,j], sd=sqrt(yPAR_sig[j,1]))   ))+dnorm(prop[cc], mu.r, sqrt(Sigma.r), log=TRUE)
      den <- sum(log(expit(yPAR_r[,j], X_lj)*identity(subset(y, Sy==j)==0)+identity(subset(y, Sy==j)!=0)*(1-expit(yPAR_r[,j], X_lj))*dtruncnorm(subset(y, Sy==j),  a=0, b=Inf, mean=X_lj%*%yPAR_beta[,j], sd=sqrt(yPAR_sig[j,1]))   ))+dnorm(yPAR_r[cc,j], mu.r, sqrt(Sigma.r), log=TRUE)
      if (log(runif(1))< (num-den) & !is.nan(num-den)){
        yPAR_r[cc,j] <- prop[cc]
      }else{ 
        prop[cc] <- yPAR_r[cc,j]
        }
    }  

    for(l in 1:length(which(uniqueYX$Sy==j))){
      ind <- which(uniqueYX$Sy==j & uniqueYX$Sx==l)
      xPAR_p0[ind,] <- MCMCpack::rdirichlet(1, alpha00+c(table(factor(subset(a+1, Sy==j & Sx==l), levels=1:length(alpha00)))))
      xPAR_p1[ind,] <- MCMCpack::rdirichlet(1, alpha01+c(table(factor(subset(x1+1, Sy==j & Sx==l), levels=1:length(alpha01)))))
      sd2 <- ifelse(is.na(sd(subset(x2, Sy==j & Sx==l))^2), 0, sd(subset(x2, Sy==j & Sx==l))^2)
      xPAR_sig[1,ind] <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x2, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR_sig[1,ind]
      xPAR_mu[1,ind] <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x2, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
      sd2 <- ifelse(is.na(sd(subset(x3, Sy==j & Sx==l))^2), 0, sd(subset(x3, Sy==j & Sx==l))^2)
      xPAR_sig[2,ind] <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x3, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR_sig[2,ind]
      xPAR_mu[2,ind] <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x3, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
      sd2 <- ifelse(is.na(sd(subset(x4, Sy==j & Sx==l))^2), 0, sd(subset(x4, Sy==j & Sx==l))^2)
      xPAR_sig[3,ind] <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x4, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR_sig[3,ind]
      xPAR_mu[3,ind] <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x4, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
  
      sd2 <- ifelse(is.na(sd(subset(x5, Sy==j & Sx==l))^2), 0, sd(subset(x5, Sy==j & Sx==l))^2)
      xPAR_sig[4,ind] <- rinvchisq(1, nu0 + length(which(Sy==j & Sx==l)), (nu0*tau0 + (length(which(Sy==j & Sx==l))-1) * sd2  +  c0*length(which(Sy==j & Sx==l))/(c0+length(which(Sy==j & Sx==l))) * (mean(subset(x5, Sy==j & Sx==l)) - mu0)^2  )/(nu0 + length(which(Sy==j & Sx==l))) )  # mu0 : 0 
      term1 <- c0 / tau0
      term2 <- length(which(Sy==j & Sx==l)) / xPAR_sig[4,ind]
      xPAR_mu[4,ind] <- rnorm(1, (term1 * mu0 + term2 * mean(subset(x5, Sy==j & Sx==l))) / (term1 + term2), sqrt(1 / (term1 + term2)))
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
  }else{
    alpha_omega <- prop
    }
  
  return(list(alpha_omega=alpha_omega,  alpha_theta=alpha_theta, xPAR_p0=xPAR_p0, xPAR_p1=xPAR_p1, xPAR_mu=xPAR_mu, xPAR_sig=xPAR_sig, yPAR_beta=yPAR_beta, yPAR_r=yPAR_r, yPAR_sig=yPAR_sig))
}
