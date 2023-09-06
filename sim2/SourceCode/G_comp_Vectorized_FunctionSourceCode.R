sim_new_cluster <- function(n_new, xall_names, xall_type, 
                            g1, b1, g2, b2, prior_means, prior_vars, 
                            beta_prior_mean, beta_prior_var,
                            gamma_prior_mean, gamma_prior_var,
                            interv_var, interv_val, interv_ref){
  
  beta_new <- rmvnorm(n = n_new, mean = beta_prior_mean, 
                      sigma = diag(beta_prior_var))
  psi_new <- invgamma::rinvgamma(n = n_new, shape = g1, rate = b1)
  gamma_new <- rmvnorm(n = n_new, mean = gamma_prior_mean, 
                       sigma = diag(gamma_prior_var))
  
  xmat <- matrix(NA, nrow = n_new, ncol = length(xall_names) + 1)
  xmat[,1] <- 1
  colnames(xmat) <- c('constant', xall_names)
  
  for(p in setdiff(xall_names,interv_var) ){
    if(xall_type[p]=='binary'){
      c_b_new <- rbeta(n = n_new, shape1 = 1, shape2 = 1)
      xmat[ , p] <- rbinom(n = n_new, size = 1, prob = c_b_new)
    }else{
      c_n_new_mu <- rnorm(n = n_new, prior_means[p], sqrt(prior_vars[p]) )
      c_n_new_phi <- invgamma::rinvgamma(n = n_new, shape = g2, rate = b2)
      xmat[ , p] <- rnorm(n = n_new, mean = c_n_new_mu, sd = sqrt(c_n_new_phi) )
    }
  }
  
  xmat_interv <- xmat_ref <- xmat
  
  xmat_interv[, interv_var] <- interv_val
  xmat_ref[, interv_var] <- interv_ref
  
  y_mu_interv <- dot_prod(xmat_interv, t(beta_new), n_new, ncol(xmat))
  y_mu_ref <- dot_prod(xmat_ref, t(beta_new), n_new, ncol(xmat))
  
  pz_interv <- invlogit(dot_prod(xmat_interv, t(gamma_new), n_new, ncol(xmat)))
  pz_ref <- invlogit(dot_prod(xmat_ref, t(gamma_new), n_new, ncol(xmat)))
  
  y_pos_draw_interv <-rnorm(n = n_new, mean = y_mu_interv, sd = sqrt( psi_new ))
  y_pos_draw_ref <- rnorm(n = n_new, mean = y_mu_ref, sd = sqrt( psi_new ))
  
  y_draw_interv <- y_pos_draw_interv*(1-pz_interv)
  y_draw_ref <- y_pos_draw_ref*(1-pz_ref)
  
  delta <- y_draw_interv - y_draw_ref
  
  res <- cbind(y_interv = y_draw_interv, 
               y_ref = y_draw_ref, 
               delta = delta)
  
  return(res)
}


sim_exist_cluster <- function(c_vec,param_shell_i, 
                              xall_names, xall_type,
                              interv_var, interv_val, interv_ref){
  n_sim <- length(c_vec)
  
  xmat <- matrix(NA, nrow = n_sim, ncol = length(xall_names) + 1)
  xmat[,1] <- 1
  colnames(xmat) <- c('constant', xall_names)
  
  for(p in setdiff(xall_names,interv_var) ){
    if(xall_type[p]=='binary'){ 
      xmat[,p] <- rbinom(n = n_sim, size = 1, prob = param_shell_i$param_c$binary[[p]][1,c_vec])
    }else{
      xmat[,p] <- rnorm(n = n_sim, 
                        mean = param_shell_i$param_c$numeric[[p]]$mu[1, c_vec], 
                        sd = sqrt(param_shell_i$param_c$numeric[[p]]$phi[1,c_vec]) )
    }
  }
  
  xmat_interv <- xmat_ref <- xmat
  
  xmat_interv[, interv_var] <- interv_val
  xmat_ref[, interv_var] <- interv_ref
  
  y_mu_interv <- dot_prod(xmat_interv, 
                          param_shell_i$param_y$beta[,c_vec, drop=F], 
                          n_sim, ncol(xmat))
  
  y_mu_ref <- dot_prod(xmat_ref, 
                       param_shell_i$param_y$beta[,c_vec, drop=F], 
                       n_sim, ncol(xmat))

  
  pz_interv <- invlogit(dot_prod(xmat_interv, 
                                 param_shell_i$param_z[ , c_vec, drop=F], 
                                 n_sim, ncol(xmat)))
  
  pz_ref <- invlogit(dot_prod(xmat_ref, 
                              param_shell_i$param_z[ , c_vec, drop=F], 
                              n_sim, ncol(xmat)))
  
  y_pos_draw_interv <-rnorm(n = n_sim, 
                            mean = y_mu_interv, 
                            sd = sqrt( param_shell_i$param_y$psi[1,c_vec] ))
  
  y_pos_draw_ref <- rnorm(n = n_sim, 
                          mean = y_mu_ref, 
                          sd = sqrt( param_shell_i$param_y$psi[1,c_vec] ))
  
  y_draw_interv <- y_pos_draw_interv*(1-pz_interv)
  y_draw_ref <- y_pos_draw_ref*(1-pz_ref)
  delta <- y_draw_interv - y_draw_ref
  
  res <- cbind(y_interv = y_draw_interv, 
               y_ref = y_draw_ref, 
               delta = delta)
  return(res)  
}


integrate_confounders <-function(iter, pvec, param_shell_i,
                                 beta_prior_mean, beta_prior_var,
                                 gamma_prior_mean, gamma_prior_var,
                                 g1, b1, g2, b2,
                                 prior_means, prior_vars, 
                                 xall_names, xall_type,
                                 interv_var, interv_val, interv_ref){
  
  c_vec <- sample(x = names(pvec), size = iter, replace = T, prob = pvec)
  n_new <- sum(c_vec=='new')
  
  if(n_new>0){
    ### integrate over new clusters
    new_res <- sim_new_cluster(n_new, xall_names, xall_type, 
                               g1, b1, g2, b2, prior_means, prior_vars, 
                               beta_prior_mean, beta_prior_var,
                               gamma_prior_mean, gamma_prior_var,
                               interv_var, interv_val, interv_ref)
    
    ### integrate over existing clusters
    exist_res <- sim_exist_cluster(c_vec = c_vec[c_vec!='new'],
                                   param_shell_i, 
                                   xall_names, xall_type,
                                   interv_var, interv_val, interv_ref)
    
    res <- colMeans(rbind(exist_res, new_res))

  }else{
  
    ### integrate over existing clusters
    exist_res <- sim_exist_cluster(c_vec = c_vec[c_vec!='new'],
                                   param_shell_i, 
                                   xall_names, xall_type,
                                   interv_var, interv_val, interv_ref)
    
    res <- colMeans(exist_res)
  }
  
  return(res)
}


bnp_standardization <- function(DPmix_res, iter,
                                interv_var, interv_val, interv_ref){
  
  ###------------------------------------------------------------------------###
  #### 0 - Parse User Inputs                                                ####
  ###------------------------------------------------------------------------###
  
  c_shell <- DPmix_res$c_shell ## n X iter_store matrix
  
  iter_store <- length(DPmix_res$param_shell)
  
  nparams <- length(DPmix_res$model_info$x_type)
  
  ## continuous covariate prior means and variances
  prior_means <- DPmix_res$model_info$prior_means
  prior_vars <- DPmix_res$model_info$prior_vars
  
  n <- nrow(c_shell) # number of subjects
  
  x_type <- DPmix_res$model_info$x_type
  
  x_names <- DPmix_res$model_info$x
  
  n_bin_p <- sum(DPmix_res$model_info$x_type=='binary')
  n_num_p <- sum(DPmix_res$model_info$x_type=='numeric')
  
  names(x_type) <- x_names
  
  xall_names <- x_names
  nparams_all <- length(xall_names)
  
  all_type_map <- c(x_type)
  names(all_type_map) <- c(x_names)
  xall_type <- all_type_map[as.character(xall_names)]
  
  n_bin_all <- sum(xall_type == 'binary')
  n_num_all <- sum(xall_type == 'numeric')
  
  xall_names_bin <- xall_names[xall_type == 'binary']
  xall_names_num <- xall_names[xall_type == 'numeric']
  
  a <- DPmix_res$alpha
  
  beta_prior_mean <- DPmix_res$model_info$beta_prior_mean
  beta_prior_var <- DPmix_res$model_info$beta_prior_var
  
  gamma_prior_mean <- DPmix_res$model_info$gamma_prior_mean
  gamma_prior_var <- DPmix_res$model_info$gamma_prior_var
  
  g1 <- DPmix_res$model_info$g1
  b1 <- DPmix_res$model_info$b1
  
  g2 <- DPmix_res$model_info$g2
  b2 <- DPmix_res$model_info$b2
  
  ###------------------------------------------------------------------------###
  #### 1 - Create Shells for storage                                        ####
  ###------------------------------------------------------------------------###
  
  res <- matrix(NA, nrow=iter_store-1, ncol=3)  
  
  ###------------------------------------------------------------------------###
  #### 2 - Cycle Through Posterior Draws                                    ####
  ###------------------------------------------------------------------------###
  

  
  for( i in 2:iter_store){
    
    pvec <- c(table(factor(c_shell[,i])), a[i])
    names(pvec)[length(pvec)] <- 'new'

    ###----------------------------------------------------------------------###
    #### 3 - Integrate over Confounders                                     ####
    ###----------------------------------------------------------------------###
    
    ## integrate using Monte Carlo procedure with "iter" iterations.
    #browser()
    res[i-1,] <-integrate_confounders(iter, pvec, 
                                param_shell_i = DPmix_res$param_shell[[i]],
                                beta_prior_mean, beta_prior_var,
                                gamma_prior_mean, gamma_prior_var,
                                g1, b1, g2, b2,
                                prior_means, prior_vars,
                                xall_names, xall_type,
                                interv_var, interv_val, interv_ref)
  }
  
  
  return(res)
}
