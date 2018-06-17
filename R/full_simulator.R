reg_horseshoe = function(zb, aux1_global, aux2_global, aux1_local, aux2_local, 
                         caux, scale_global, slab_scale) {
  K = length(zb)
  
  lambda = aux1_local * sqrt ( aux2_local );
  tau = aux1_global * sqrt ( aux2_global ) * scale_global * 1 ;
  c = slab_scale * sqrt ( caux );
  lambda_tilde = sqrt ( c ^2 * lambda ^2 / (c ^2 + tau ^2* lambda^2) );
  zb * lambda_tilde * tau ;  
}

horseshoe_hyperprior_params <- function(tau) {
  #Using data.frame because it avoids the need to name params. list(tau = tau) would work as well
  data.frame(type = "horseshoe",tau) 
}

default_intercept_params <- function(mean, sigma) {
  data.frame(mean, sigma)
}

normal_likelihood_params <- function(xi_sigma) {
  data.frame(type = "normal", xi_sigma)
}
  
exp_link_params <- function() {
  data.frame(type = "exp")
}

linear_link_params <- function() {
  data.frame(type = "linear")
}

full_simulator = function(n_genes, n_tubes, hyperprior_params, intercept_params, link_params, likelihood_params) {
  observed = list()
  true = list()
  
  # 0. Draw covarietes - assuming disease progression here
  observed$X = runif(n_tubes, 0, 1)
  observed$G = n_genes
  observed$T = n_tubes
  
  # 1. Draw from hyperprior
  if(hyperprior_params$type == "horseshoe") {
    true$horseshoe_sigma = abs(rnorm(1, 0, hyperprior_params$tau))
    true$beta1 = rnorm(n_genes, 0, true$horseshoe_sigma)
    
    observed$horseshoe_tau = hyperprior_params$tau
  } else {
    stop("Unrecognized hyperprior type: ",hyperprior_params$type)
  }
  
  # 2. Draw the intercepts
  true$intercept = rnorm(n_genes, intercept_params$mean, intercept_params$sigma)
  
  observed$intercept_mean = intercept_params$mean
  observed$intercept_sigma = intercept_params$sigma
  
  # 3. Link function
  linear_predictor = true$intercept + observed$X %*% t(true$beta1)
  true$linear_predictor = linear_predictor
  if(link_params$type == "linear") {
    y_hat = linear_predictor
    observed$link_type = 0
  } else if(link_params$type == "exp") {
    y_hat = exp(linear_predictor)
    observed$link_type = 1
  } else {
    stop("Unrecognized link function:", link_params$type)
  }
  
  if(likelihood_params$type == "normal") {
    true$xi = abs(rnorm(1, 0, likelihood_params$xi_sigma))
    
    observed$xi_sigma = likelihood_params$xi_sigma
    observed$y_real = rnorm(n_tubes * n_genes, y_hat, true$xi) %>% matrix(ncol = n_genes, nrow = n_tubes)
    observed$likelihood_type = 0
  } else {
    stop("Unrecognized likelihood type:", likelihood_params$type)
  }
  
  #scale_global =  par_ratio / sqrt(n_tubes)

  # beta1[1,] = reg_horseshoe(beta1_z[1,], aux1_global ,
  #                           aux2_global,
  #                           aux1_local ,
  #                           aux2_local ,
  #                           caux,
  #                           scale_global,
  #                           slab_scale
  # )
  list(true = true, observed = observed)
}