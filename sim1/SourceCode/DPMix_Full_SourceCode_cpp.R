rcond_post_beta <- function(nparams, y, xm, psi, Imat, 
                            beta_prior_mean, beta_prior_var){
  
  mu_beta <- beta_prior_mean
  v_beta <- diag(beta_prior_var)
  v_beta_inv <- diag(1/beta_prior_var)
  
  xtx <- t(xm)%*%xm
  
  if(length(y)>0 ){
    
    post_cov <- chol2inv(chol( v_beta_inv + (1/psi)*xtx))
    post_mean <- post_cov %*% (v_beta_inv %*% mu_beta + (1/psi)*t(xm)%*%y )
    
  }else{
    post_mean <- mu_beta
    post_cov <- v_beta
  }
  
  draw <- rmvnorm(n = 1, mean = post_mean, sigma = post_cov)
  return(draw)  
}



rcond_post_psi <- function(beta, y, xm, g1, b1){
  n_k <- length(y)
  if(n_k==0){ 
    shape_k <- g1
    rate_k <- b1
  }else{
    shape_k <- g1 + n_k/2
    rate_k <- .5*(  sum( (y - xm %*% beta)^2 ) )  + b1
  }
  draw <- invgamma::rinvgamma(n = 1, shape = shape_k, rate = rate_k)
  return(draw)
}

rcond_post_mu_x <- function(x_ck, phi_x, lambda, tau ){
  nvec_k <- length(x_ck)
  mu_mean <- (1/(1/tau + nvec_k/phi_x))*(lambda/tau + sum(x_ck)/phi_x)
  mu_sd <- sqrt( (1/tau + nvec_k/phi_x)^(-1) )
  draw <- rnorm(n = 1, mean = mu_mean, sd =  mu_sd)
  return(draw)
}

rcond_post_phi <- function(mu_x, x, g2, b2){
  n_k <- length(x)
  if(n_k==0){ 
    shape_k <- g2
    rate_k <- b2
  }else{
    shape_k <- g2 + n_k/2
    rate_k <- .5*(  sum( (x - mu_x )^2 ) )  + b2
  }
  phi_post <- invgamma::rinvgamma(n = 1, shape = shape_k, rate = rate_k)
  return(phi_post)
}

rcond_post_mu_trt <- function(x_ck){
  draw<-rbeta(n = 1,shape1 = 1 + sum(x_ck), 
              shape2 = 1 + length(x_ck) - sum(x_ck))
  return(draw)
}


cond_post_alpha <- function(alpha, n, K, nvec){
  lik <- lgamma(alpha) - lgamma(n + alpha) + sum(lgamma(nvec + alpha/K)) - K*lgamma(alpha/K)
  #pr <- invgamma::dinvgamma(x = alpha, 1, .001,log = T)
  pr <- invgamma::dinvgamma(x = alpha, 1, 1,log = T)
  #pr <- invgamma::dinvgamma(x = alpha, 10, 1000,log = T)
  post <- lik + pr
  return(post)
}


cond_post_pz <- function(gamma, y, xp, 
                         gamma_prior_mean,gamma_prior_var){

  eta <- xp %*% gamma
  eta <- ifelse(eta < -10, -10, ifelse(eta>10, 10, eta))
  p <- invlogit(eta)
  
  lk <- sum(dbinom(x = y, size = 1, prob = p, log = T))
  pr <- dmvnorm(x = gamma, gamma_prior_mean, diag(gamma_prior_var), log = T)
  
  return(lk + pr)
}

metrop_hastings<-function(x_0, iter=1, log_post_density,
                          proposal_dist = function(x, prop_sigma){ 
                            mvrnorm(1, mu = x, Sigma = prop_sigma )
                          }, 
                          lower=-Inf, upper=Inf, prop_sigma,
                          ... ){
  for(i in 1:iter){
    # draw from proposal distribution
    x_star <- proposal_dist(x_0, prop_sigma) 
    
    # calculate ratio of conditional posterior densities
    r_num <- do.call(log_post_density, c(list(x_star), list(...)) )
    r_denom <- do.call(log_post_density, c(list(x_0), list(...)) )
    r <- exp(r_num - r_denom)
    rmin<-min(r,1)
    if(is.na(rmin)) browser()
    # accept / reject proposal
    if(rbinom(1,1,rmin)==1){ 
      x_0<-x_star
    }
  }
  
  res<-list(x_0 = x_0, accept_prob = rmin )
  return(res)
}

create_shells <- function(xall_names, all_type_map, iter_store){
  
  param_shell <- vector(mode = 'list', length = iter_store)
  n_covars <- length(xall_names)
  
  n_bin <- sum(all_type_map=='binary')
  n_num <- sum(all_type_map=='numeric')
  
  for(l in 1:iter_store ){
    param_shell[[l]] <- vector(mode = 'list', length = 4)
    names(param_shell[[l]]) <- c('param_y','param_z', 'param_trt','param_c')
    
    ## shells for outcome model, zero prob model, and covariate dists.
    param_shell[[l]][['param_y']] <- vector(mode = 'list', length = 2)
    names(param_shell[[l]][['param_y']]) <- c('beta', 'psi')
    
    param_shell[[l]][['param_c']] <- vector(mode = 'list', length = 2)
    names(param_shell[[l]][['param_c']]) <- c('binary','numeric')
    
    param_shell[[l]][['param_c']][['binary']] <- vector(mode = 'list', length = n_bin )
    names(param_shell[[l]][['param_c']][['binary']]) <- xall_names[all_type_map=='binary']
    
    param_shell[[l]][['param_c']][['numeric']] <- vector(mode='list', length = n_num)
    names(param_shell[[l]][['param_c']][['numeric']]) <- xall_names[all_type_map=='numeric']
  } 
  
  return(param_shell)
}

DPglmMix<-function(d, y, x, trt_name='A', burnin, gibbs_iter, trim, K, a=1,
                   x_type,
                   g1=1, b1=.1,
                   beta_prior_mean, beta_prior_var,
                   gamma_prior_mean, gamma_prior_var,
                   eta_prior_mean, eta_prior_var){
  
  ###------------------------------------------------------------------------###
  #### 0 - Parse User Inputs                                                ####
  ###------------------------------------------------------------------------###
  store_l <- seq(burnin, gibbs_iter, trim)
  iter_store <- length(store_l)
  y <- d[,y]
  z <- as.numeric(y==0)
  trt <- d[,trt_name]
  x_names <- x
  x <- as.matrix(cbind(1, d[,x_names]))
  
  nparams <- ncol(x)
  n<-nrow(x) 
  
  xall_names <- x_names
  nparams_all <- length(xall_names)
  xall <- d[,xall_names]
  
  all_type_map <- c(x_type)
  names(all_type_map) <- c(x_names)
  xall_type <- all_type_map[as.character(xall_names)]
  
  xall_names_bin <- xall_names[xall_type == 'binary']
  xall_names_num <- xall_names[xall_type == 'numeric']
  
  # parse numeric versus binary variables for outcome model
  n_bin_p <- sum(x_type=='binary')
  n_num_p <- sum(x_type=='numeric')
  
  names_binminA <- xall_names_bin[xall_names_bin!='A']
  names_minA <- setdiff(xall_names, 'A')
  trt_mat_subset <- c('1',names_minA)
  
  cat_idx <- match(names_minA, colnames(x))[ all_type_map[names_minA]=='binary'  ] - 1
  num_idx <- match(names_minA, colnames(x))[ all_type_map[names_minA]=='numeric'  ] - 1
  
  ###------------------------------------------------------------------------###
  #### 1 - Create Shells for storing Gibbs Results                          ####
  ###------------------------------------------------------------------------###
  curr_class_list <- 1:K
  class_names <- paste0('c',curr_class_list)
  
  c_shell<-matrix(NA, nrow=n, ncol=1) # shell for indicators
  
  class_shell <- matrix(NA, nrow=n, ncol=iter_store)
  colnames(class_shell) <- store_l
  
  psi_shell <- matrix(NA, nrow= 1, ncol = K) # shell for Y's cond. var
  colnames(psi_shell) <- class_names
  
  beta_shell <- matrix(NA, nrow = nparams, ncol = K)
  colnames(beta_shell) <- class_names
  rownames(beta_shell) <- colnames(x)
  
  gamma_shell <- matrix(NA, nrow = nparams, ncol = K)
  colnames(gamma_shell) <- class_names
  rownames(gamma_shell) <- colnames(x)
  
  eta_shell <- matrix(NA, nrow = nparams-1, ncol = K)
  colnames(eta_shell) <- class_names
  rownames(eta_shell) <- colnames(x)[colnames(x)!=trt_name]
  
  c_b_shell <- vector(mode = 'list', length = n_bin_p-1) # -1 for trt
  names(c_b_shell) <- xall_names_bin[xall_names_bin!='A']
  
  c_n_shell <- vector(mode = 'list', length = n_num_p)
  names(c_n_shell) <- xall_names_num
  
  for( p in xall_names_bin[xall_names_bin!='A'] ){
    c_b_shell[[p]] <- matrix(NA, nrow = 1, ncol = K)
    
    c_b_shell[[p]][1,] <- rep(.5, K)
    
    colnames(c_b_shell[[p]]) <- class_names
  }
  
  for( p in xall_names_num ){
    c_n_shell[[p]] <- list(matrix(NA, nrow=1, ncol=K),
                           matrix(NA, nrow=1, ncol=K))
    
    c_n_shell[[p]][[2]][1, ] <- rep(0, K) 
    c_n_shell[[p]][[2]][1, ] <- rep(20^2, K)
    
    colnames(c_n_shell[[p]][[1]]) <- class_names
    colnames(c_n_shell[[p]][[2]]) <- class_names
  }
  
  c_b_new <- numeric(length = n_bin_p-1)
  names(c_b_new) <- names_binminA
  
  c_n_new <- matrix(NA, nrow = 2, ncol = n_num_p)
  colnames(c_n_new) <- xall_names_num
  
  param_shell <- create_shells(xall_names, xall_type, iter_store)
  names(param_shell) <- store_l
  
  pp_shell <- vector(mode = 'list', length = iter_store)
  names(pp_shell) <- store_l
  
  y_isp <- numeric(length = n)
  z_isp <- numeric(length = n)
  pz_isp <- numeric(length = n)
  
  prop_score <- matrix(NA, nrow = n, ncol = iter_store)
  colnames(prop_score) <- store_l
  
  
  ###------------------------------------------------------------------------###
  #### 2 - Set Initial Values                                               ####
  ###------------------------------------------------------------------------###
  
  # initially assign everyone to one of K clusters randomly with uniform prob
  c_shell[,1] <- sample(x = class_names, size = n, prob = rep(1/K,K),
                        replace = T)
  
  beta_start <- matrix(0, nrow=nparams, ncol=K)
  colnames(beta_start) <- class_names
  
  gamma_start <- matrix(0, nrow=nparams, ncol=K)
  colnames(gamma_start) <- class_names
  
  eta_start <- matrix(0, nrow=nparams-1, ncol=K)
  colnames(eta_start) <- class_names
  
  for(k in class_names){
    beta_shell[,k] <- beta_prior_mean
    gamma_shell[,k] <- gamma_start[,k]
    eta_shell[,k] <- eta_start[,k]
  }
  
  for(p in xall_names_bin[xall_names_bin!='A'] ){
    c_b_shell[[p]][1,] <- rep(1, K)
  }
  
  for(p in xall_names_num ){
    c_n_shell[[p]][[2]][ 1 , ] <- rep(0, K) 
    c_n_shell[[p]][[2]][ 1 , ] <- rep(20^2, K)
  }
  
  prop_sigma_z <- diag(rep(.005, nparams))
  prop_sigma_trt <- diag(rep(.05, nparams-1 ))
  
  psi_shell[1,1:K] <- rep(500000, K)
  
  # want to change this later to update alpha
  alpha <- numeric(length = iter_store)
  names(alpha) <- store_l
  
  prior_means <- apply(x[,xall_names_num, drop=F], 2, mean)
  names(prior_means) <- xall_names_num
  
  prior_var <- apply(x[,xall_names_num, drop=F], 2, var)
  names(prior_var) <- xall_names_num
  
  g2 <- rep(2, n_num_p)
  names(g2) <- xall_names_num
  
  b2 <- apply(x[,xall_names_num, drop=F], 2, var)
  names(b2) <- xall_names_num
  
  
  ###------------------------------------------------------------------------###
  #### 3 - Gibbs Sampler                                                    ####
  ###------------------------------------------------------------------------###
  
  for(i in 2:gibbs_iter){
    
    #print(i)
    
    # compute size of each existing cluster
    class_ind<-c_shell[,1]
    nvec<-table(factor(class_ind,levels = class_names) )

    ###----------------------------------------------------------------------###
    #### 3.1 - update Concentration Parameter                               ####
    ###----------------------------------------------------------------------###
    
    a_star <- rnorm(n = 1, mean =  a, sd = 1)
    if(a_star>0){
      r_num <- cond_post_alpha(a_star, n, K, nvec)
      r_denom <- cond_post_alpha(a, n, K, nvec)
      r <- exp(r_num - r_denom)
      accept <- rbinom(1,1, min(r,1) )==1
      if(accept){ a <- a_star }
    }
    
    ###----------------------------------------------------------------------###
    #### 3.2 - update parameters conditional on classification              ####
    ###----------------------------------------------------------------------###
    K <- length(unique(class_names))
    for(k in class_names){ # cycle through existing clusters - update parameters
      
      ck_ind <- class_ind == k
      
      y_ck <- y[ck_ind]
      
      x_ck <- x[ck_ind,, drop=FALSE]
      trt_ck <- trt[ck_ind]
      
      psi_shell[1, k]<-rcond_post_psi(beta = beta_shell[,k, drop=F], 
                                      y = y_ck[y_ck>0], xm = x_ck[y_ck>0, , drop=F], 
                                      g1 = g1, b1 = b1)
      
      beta_shell[,k] <- rcond_post_beta(y=y_ck[y_ck>0], 
                                        xm=x_ck[y_ck>0, , drop=F], 
                                        psi = psi_shell[1,k],
                                        beta_prior_mean = beta_prior_mean, 
                                        beta_prior_var = beta_prior_var)
      
      for( p in names_binminA ){
        c_b_shell[[p]][1, k] <- rcond_post_mu_trt(x_ck = x_ck[, p])
      }
      
      for( p in xall_names_num ){
        c_n_shell[[p]][[1]][1,k] <- rcond_post_mu_x(x_ck = x_ck[, p], 
                                                    phi_x = c_n_shell[[p]][[2]][1,k], 
                                                    lambda = prior_means[p], 
                                                    tau = prior_var[p])
        
        c_n_shell[[p]][[2]][1,k] <- rcond_post_phi(mu_x = c_n_shell[[p]][[1]][1, k], 
                                                   x = x_ck[ , p], 
                                                   g2 = g2[p], b2 = b2[p])
      
      }
      
      gamma_shell[ ,k] <- metrop_hastings(x_0 = gamma_shell[ ,k], iter = 1,
                                          log_post_density = cond_post_pz,
                                          y= y_ck==0, xp = x_ck, 
                                          prop_sigma = prop_sigma_z,
                                          gamma_prior_mean = gamma_prior_mean, 
                                          gamma_prior_var = gamma_prior_var)[[1]]
      
      eta_shell[ ,k] <- metrop_hastings(x_0 = eta_shell[ ,k], iter = 1,
                                        log_post_density = cond_post_pz,
                                        y = trt_ck, xp = x_ck[,trt_mat_subset], 
                                        prop_sigma = prop_sigma_trt, 
                                        gamma_prior_mean = eta_prior_mean, 
                                        gamma_prior_var = eta_prior_var)[[1]]
    }
    
    ###----------------------------------------------------------------------###
    #### 3.3 - update classification conditional on clusters                ####
    ###----------------------------------------------------------------------###
    
    ## draw parameters for potentially new cluster from priors
    ### draw beta, gamma, phi, covariate parameters
    
    beta_new <- rmvnorm(n = 1, mean = beta_prior_mean, 
                        sigma = diag(beta_prior_var) )
    psi_new <- invgamma::rinvgamma(n = 1, shape = g1, rate = b1)
    gamma_new <- rmvnorm(n = 1, mean = gamma_prior_mean, 
                         sigma = diag(gamma_prior_var))
    eta_new <- rmvnorm(n = 1, mean = eta_prior_mean, 
                         sigma = diag(eta_prior_var))
    
    for(p in names_binminA){
      c_b_new[p] <- rbeta(n = 1, shape1 = 1, shape2 = 1)
    }
    
    for(p in xall_names_num){
      c_n_new[1, p] <- rnorm(n = 1, prior_means[p], sqrt(prior_var[p]) )
      c_n_new[2, p] <- invgamma::rinvgamma(n = 1, shape = g2[p], rate = b2[p])
    }
    
    name_new <- paste0('c', max(as.numeric(substr(class_names, 2,10))) + 1)
    
    ## compute post prob of membership to existing and new cluster, 
    ## then classify
    
    #print(table(class_ind))
    weights <- class_update(n = n, K = length(class_names), alpha = a, 
                            name_new =  name_new,
                            uniq_clabs = colnames(beta_shell), clabs = class_ind, 
                            
                            y = y, x = x, z = z, 
                            trt = trt, x_trt=x[, trt_mat_subset ],
                            x_cat_shell = c_b_shell , x_num_shell = c_n_shell,
                            
                            ## col number of each covar...index with base 0
                            cat_idx=names_binminA, num_idx=xall_names_num,
                            
                            beta_shell = beta_shell, 
                            psi_shell = psi_shell, 
                            gamma_shell = gamma_shell, 
                            eta_shell = eta_shell,
                            
                            beta_new=beta_new, psi_new=psi_new, 
                            cat_new = c_b_new, num_new=c_n_new, 
                            gamma_new=gamma_new, eta_new=eta_new)
    
    c_shell[,1] <- apply(weights, 1, FUN = sample, 
                         x=c(colnames(beta_shell), name_new), size=1, replace=T)
    
    #browser()
    ###----------------------------------------------------------------------###
    #### FOR TESTING: Output Posterior Predictive Draws                     ####
    ###----------------------------------------------------------------------###
    tabb <- table(c_shell)
    pc <- sample(names(tabb), 1000, T, prob = tabb )

    x_pp <- matrix(NA, nrow = 1000, ncol = ncol(x) )
    colnames(x_pp) <- colnames(x)
    x_pp[,'1'] <- 1

    A_pp <- numeric(length = 1000)
    y_pp <- numeric(length = 1000)
    zp_pp <- numeric(length = 1000)

    new_idx <- pc==name_new
    n_new <- sum(new_idx)
    exist_idx <- pc!=name_new
    n_exist <- sum(exist_idx)

    exist_labs <- unique(pc[exist_idx])


    ## draw from posterior predictive of covariates
    for(p in names_binminA ){
      x_pp[new_idx, p] <- rbinom(n = n_new, size = 1, prob = c_b_new[p] )

      x_pp[exist_idx, p] <- rbinom(n = n_exist,
                                   size = 1,
                                   prob = c_b_shell[[p]][,pc[exist_idx] ] )
    }

    for(p in xall_names_num ){
      x_pp[new_idx, p] <- rnorm(n = n_new,
                                mean = c_n_new[1, p],
                                sd =  sqrt(c_n_new[2, p]) )
      x_pp[exist_idx, p] <- rnorm(n = n_exist,
                                  mean = c_n_shell[[p]][[1]][,pc[exist_idx]],
                                  sd =  sqrt(c_n_shell[[p]][[2]][,pc[exist_idx]]) )
    }

    ## draw posterior predictive treatment
    A_pp[new_idx] <- rbinom(n = n_new, 1, prob =  invlogit( x_pp[new_idx, trt_mat_subset] %*% t(eta_new)  ) )
    for( k in exist_labs  ){
      A_pp[pc==k] <- rbinom(n = sum(pc==k), 1,
                            prob =  invlogit( x_pp[pc==k, trt_mat_subset] %*% eta_shell[,k, drop=F]  ) )
    }
    x_pp[,'A'] <- A_pp

    # draw posterior predictive cost
    y_pp[new_idx] <- rnorm(n = n_new, mean =  x_pp %*% t(beta_new) ,sd = sqrt(psi_new) )
    for( k in exist_labs  ){
      y_pp[pc==k] <- rnorm(n = sum(pc==k),
                           mean = x_pp[pc==k,] %*% beta_shell[,k,drop=F],
                           sd = sqrt(psi_shell[1,k]) )
    }

    # draw posterior predictive cost
    zp_pp[new_idx] <- rbinom(n = n_new, 1, prob =  invlogit( x_pp[new_idx,] %*% t(gamma_new)  ) )
    for( k in exist_labs  ){
      zp_pp[pc==k] <- rbinom(n = sum(pc==k), 1,
                             prob =  invlogit( x_pp[pc==k, ] %*% gamma_shell[,k, drop=F]  ) )
    }
    #if(i==1500) browser()
    y_pp[zp_pp==1] <- 0
         
    ###----------------------------------------------------------------------###
    #### 4.0 - Store Results                                                ####
    ###----------------------------------------------------------------------###
    
    ###----------------------------------------------------------------------###
    #### 5.0 - update shell dimensions as clusters die/ are born            ####
    ###----------------------------------------------------------------------###
    
    new_class_names <- unique(c_shell[,1])
    K_new <- length(new_class_names)
    
    new_class_name <- setdiff(new_class_names, class_names)
    if(length(new_class_name)>0){ # new cluster was born
      
      beta_shell <- cbind(beta_shell, t(beta_new) )
      colnames(beta_shell)[K+1] <- new_class_name
      
      psi_shell <- cbind(psi_shell, psi_new)
      colnames(psi_shell)[K+1] <- new_class_name
      
      gamma_shell <- cbind(gamma_shell, t(gamma_new) )
      colnames(gamma_shell)[K+1] <- new_class_name
      
      eta_shell <- cbind(eta_shell, t(eta_new) )
      colnames(eta_shell)[K+1] <- new_class_name
      
      for(p in names_binminA ){
        c_b_shell[[p]] <- cbind( c_b_shell[[p]], c_b_new[p] )
        colnames(c_b_shell[[p]])[K+1] <- new_class_name
      }
      
      for(p in xall_names_num ){
        c_n_shell[[p]][[1]] <- cbind(c_n_shell[[p]][[1]], c_n_new[1, p])
        colnames(c_n_shell[[p]][[1]])[K+1] <- new_class_name
        
        c_n_shell[[p]][[2]] <- cbind(c_n_shell[[p]][[2]], c_n_new[2, p])
        colnames(c_n_shell[[p]][[2]])[K+1] <- new_class_name
      }
      
    }
    
    dead_classes <- setdiff(class_names, new_class_names)
    if(length(dead_classes) > 0){ # clusters have died
      beta_shell <- beta_shell[, ! colnames(beta_shell) %in% dead_classes, drop=F]
      gamma_shell <- gamma_shell[, ! colnames(gamma_shell) %in% dead_classes, drop=F]
      eta_shell <- eta_shell[, ! colnames(eta_shell) %in% dead_classes, drop=F]
      
      psi_shell <- psi_shell[, ! colnames(psi_shell) %in% dead_classes, drop=F]
      
      for(p in names_binminA){
        c_b_shell[[p]] <- c_b_shell[[p]][, ! colnames(c_b_shell[[p]]) %in% dead_classes, drop=F]
      }
      
      for(p in xall_names_num){
        c_n_shell[[p]][[1]] <- c_n_shell[[p]][[1]][, ! colnames(c_n_shell[[p]][[1]]) %in% dead_classes, drop=F]
        c_n_shell[[p]][[2]] <- c_n_shell[[p]][[2]][, ! colnames(c_n_shell[[p]][[2]]) %in% dead_classes, drop=F]
      }
      
    }
    
    K <- K_new
    class_names <- new_class_names
    
    for(k in unique(c_shell[,1])){
      y_isp[c_shell[,1]==k] <- rnorm(n = sum(c_shell[,1]==k), 
                                     mean = x[c_shell[,1]==k, , drop=F] %*% beta_shell[ , k , drop=F], 
                                     sd = sqrt(psi_shell[1,k]) )
      
      z_isp[c_shell[,1]==k] <- rbinom(n = sum(c_shell[,1]==k), 1,
                                      prob =  invlogit( x[c_shell[,1]==k, ] %*% gamma_shell[, k, drop=F]  ) )
      
      pz_isp[c_shell[,1]==k] <- invlogit( x[c_shell[,1]==k, ] %*% gamma_shell[, k, drop=F]  )
    }
    y_isp[z_isp==1] <- 0 
    
    
    if( i >= burnin ){
      colcount <- as.character(i)
      
      param_shell[[ colcount ]]$param_y$beta <- beta_shell
      param_shell[[ colcount ]]$param_y$psi <- psi_shell
      param_shell[[ colcount ]]$param_z <- gamma_shell
      param_shell[[ colcount ]]$param_trt <- eta_shell
      
      for(p in names_binminA){
        param_shell[[ colcount ]]$param_c$binary[[p]] <- c_b_shell[[p]][1, , drop=F]
      }
      
      for(p in xall_names_num){
        param_shell[[ colcount ]]$param_c$numeric[[p]]$mu <- c_n_shell[[p]][[1]][1, , drop=F]
        param_shell[[ colcount ]]$param_c$numeric[[p]]$phi <- c_n_shell[[p]][[2]][1, , drop=F]
      }
      
      class_shell[, colcount] <- c_shell[,1]
      
      alpha[colcount] <- a
      
      pp_shell[[colcount]] <- list(y_pp = y_pp, x_pp = x_pp, 
                                   y_isp=y_isp, z_isp=z_isp, pz_isp = pz_isp)
      
      for(t in 1:n){ 
        if(c_shell[t,1] %in% colnames(eta_shell)){
          ptrt <- x[t, trt_mat_subset] %*% eta_shell[, c_shell[t,1] , drop=F]
        }else{
          ptrt <- x[t, trt_mat_subset] %*% t(eta_new)
        }
        prop_score[t, colcount] <- invlogit(ptrt)
      }
      
      
    }
    
    
  }
  
  results <- list(c_shell=class_shell,
                  param_shell = param_shell,
                  model_info = list(gibbs_iter = gibbs_iter, 
                                    burnin = burnin, trim=trim,
                                    K=K,
                                    x=x_names,
                                    x_type=x_type,
                                    g1=g1, b1=b1, g2=g2, b2=b2,
                                    beta_prior_mean=beta_prior_mean, 
                                    beta_prior_var = beta_prior_var,
                                    prior_means = prior_means,
                                    prior_vars = prior_var,
                                    gamma_prior_mean=gamma_prior_mean,
                                    gamma_prior_var=gamma_prior_var,
                                    eta_prior_mean=eta_prior_mean, 
                                    eta_prior_var=eta_prior_var),
                  pp = pp_shell,
                  prop_score=prop_score,
                  alpha = alpha)
  
  return(results)
}

calc_adjmat <- function(c_shell){
  iter<-ncol(c_shell)
  n <- nrow(c_shell)
  
  adjmat<-matrix(0, nrow=n, ncol=n)
  
  ## compute frequency matrix 
  for(i in 1:iter){
    adjmat_i <- matrix(0, nrow=n, ncol=n)
    class <- c_shell[,i]
    
    for(j in 1:n){
      adjmat_i[j,] <- class==class[j]
    }
    adjmat <- adjmat + adjmat_i
  }
  
  ## choose best cluster assignment
  err_vec <- numeric(length = ncol(c_shell))
  for(i in 1:iter){
    adjmat_i <- matrix(0, nrow=n, ncol=n)
    class <- c_shell[,i]
    for(j in 1:n){
      adjmat_i[j,] <- class==class[j]
    }
    err_vec[i] <- sum((adjmat_i - adjmat)^2)
  }
  
  class_mem <- c_shell[,c(1:ncol(c_shell))[ err_vec==min(err_vec) ] ]
  
  adjmat <- adjmat/iter
  return(list(adjmat=adjmat, class_mem=class_mem))
}