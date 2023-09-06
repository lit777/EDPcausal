####################
#### Scenario 1 ####
####################

library(MASS)
library(msm)
library(mvtnorm)
library(LaplacesDemon)
library(Rcpp)
library(RcppArmadillo)
library(igraph)
library(extraDistr)
library(truncnorm) 
sourceCpp('SourceCode/GcompSupport.cpp')
source("SourceCode/DPMix_Full_SourceCode_cpp_new.R")
source("SourceCode/class_update.R")
source("SourceCode/G_comp_Vectorized_FunctionSourceCode.R")

a<-y<-Y1<-Y0<-Pi0<-Pi1<-p0<-p1<-Y1M0<-p10<-Pi10<-NULL

load("covariate.RData")

for(process in 1:200){
        
n<-600 

Expit <- function(x){
  exp(x)/(1+exp(x))
}

for(i in 1:600){
  a[i] <- rbinom(1,1,Expit(1-1*x1[i]+2*x2[i]-1*abs(x3[i]-1)+0.4*x4[i]*x5[i]))
  p0[i] <- Expit(-3.5+0.2*x1[i]-0.5*x3[i]-3.5*I(x2[i]>0)*0+0.2*x4[i]+0.2*x5[i])
  p1[i] <- Expit(-3.5+0.2*x1[i]-0.5*x3[i]-3.5*I(x2[i]>0)*1+0.2*x4[i]+0.2*x5[i])
  Pi0[i] <- sample(c(0,1), 1, prob=c(1-p0[i], p0[i]))
  Pi1[i] <- sample(c(0,1), 1, prob=c(1-p1[i], p1[i]))
  rn1 <- rnorm(1,  mean=0.4*x1[i]+2*x2[i]+1*x2[i]^2*1-2.5*1+0.3*abs(x3[i]+1)+0.4*x4[i]+exp(0.1*x5[i]), sd=0.1)
  rn0 <- rnorm(1,  mean=0.4*x1[i]+2*x2[i]+1*x2[i]^2*0-2.5*0+0.3*abs(x3[i]+1)+0.4*x4[i]+exp(0.1*x5[i]), sd=0.1)
  Y1[i] <- ifelse(Pi1[i]==1, 0, ifelse(rn1 < 0, 0, rn1))
  Y0[i] <- ifelse(Pi0[i]==1, 0, ifelse(rn0 < 0, 0, rn0))
}


y <- a*Y1+(1-a)*Y0
x <- cbind(x1,x2,x3,x4, x5)

d <- data.frame(L=x, Y = y, A=a)

####------------------------------------------------------------------------####
####                    Run DP Zero-Inflated Model                          ####
####------------------------------------------------------------------------####

burnin <- 10000
gibbs_iter <- 20000

reg <- lm(data=d[d$Y>0,], Y ~ L.x1+L.x2+ L.x3+L.x4+L.x5++ A)
beta_prior_mean <- reg$coefficients
beta_prior_var <- 10000*diag(vcov(reg))

pr_mean <- c(0,0,0,0,0,0,0)
pr_var <- c(2,2,2,2,2,2,2)

DPglm_res<-DPglmMix(d = d, y = 'Y',
                    x = c('L.x1','L.x2','L.x3','L.x4','L.x5','A'),
                    x_type=c('binary','numeric','numeric','numeric','numeric','binary'),
                    burnin=burnin, gibbs_iter = gibbs_iter, trim=1,
                    K = 5, a = 1,
                    g1=10, b1=10000,
                    beta_prior_mean = beta_prior_mean,
                    beta_prior_var = beta_prior_var,
                    gamma_prior_mean = pr_mean, gamma_prior_var = pr_var,
                    eta_prior_mean = pr_mean[-1], eta_prior_var = pr_var[-1] )

## compute posterior model cluster assignment and poster mode adjacency matrix
adjmat <- calc_adjmat(DPglm_res$c_shell)

temp1 <- matrix(nrow=10000, ncol=600)
temp0 <- matrix(nrow=10000, ncol=600)
for(t in 1:10000){
  for(tt in 1:2){
    eval(parse(text=(paste0('temp1[t, ]<-DPglm_res$pp$"',t+10000,'"$y_isp1'))))
    eval(parse(text=(paste0('temp0[t, ]<-DPglm_res$pp$"',t+10000,'"$y_isp0'))))
  }
}

print(process)
save(temp1,temp0, Y1,Y0, file=paste("out_jason",process,".RData",sep=""))
}

