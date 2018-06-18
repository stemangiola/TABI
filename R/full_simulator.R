##########################
#' Calculating log(1+exp(x)) accurately
#' 
#' @description Calculates \code{log(1+exp(x))} in a numerically stable fashion.
#'  
#' @param x a numeric vector. 
#' @return A numeric vector where the i-th entry is equal to \code{log(1+exp(x[i]))}, but computed more stably.
#' @details We follow the recipe of Machler (2012), that is formula (10) page 7. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @references Machler, M. (2012). Accurately computing log(1-exp(-|a|)). 
#'             URL: \url{https://cran.r-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf}.
#' @examples
#' set.seed(141)
#' library(qgam); 
#' x <- rnorm(100, 0, 100)
#' log1pexp(x) - log1p(exp(x))
log1pexp <- function(x)
{
  indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), T)
  
  kk <- which(indx==1)
  if( length(kk) ){  x[kk] <- exp(x[kk])  }
  
  kk <- which(indx==2)
  if( length(kk) ){  x[kk] <- log1p( exp(x[kk]) ) }
  
  kk <- which(indx==3)
  if( length(kk) ){  x[kk] <- x[kk] + exp(-x[kk]) }
  
  return(x)
}

## from http://tr.im/hH5A
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

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

no_hyperprior_params <- function(sigma) {
  data.frame(type = "no", sigma)
}

default_intercept_params <- function(mean, sigma) {
  data.frame(mean, sigma)
}

normal_likelihood_params <- function(xi_sigma) {
  data.frame(type = "normal", xi_sigma)
}

lognormal_likelihood_params <- function(xi_sigma) {
  data.frame(type = "lognormal", xi_sigma)
}

neg_binomial_likelihood_params <- function(xi_sigma) {
  data.frame(type = "neg_binomial", xi_sigma)
}

dirichlet_multinom_likelihood_params <- function(xi_sigma) {
  data.frame(type = "dirichlet_multinom", xi_sigma)
}


linear_link_params <- function() {
  data.frame(type = "linear")
}

exp_link_params <- function() {
  data.frame(type = "exp")
}

#' With link sigmoid intercept is interpreted as the location of the inflection point
sigmoid_link_params <- function(mean_expression_mean, mean_expression_sigma) {
  data.frame(type = "sigmoid", mean_expression_mean, mean_expression_sigma)
}

logsigmoid_link_params <- function(inversion_mean, inversion_sigma) {
  data.frame(type = "logsigmoid", inversion_mean, inversion_sigma)
}

##########################
#' Simulate all possible hyperprior-link-likelihood combinations

full_simulator = function(n_genes, n_tubes, hyperprior_params, intercept_params, link_params, likelihood_params) {
  observed = list()
  true = list()
  
  true$hyperprior_params = hyperprior_params
  true$intercept_params = intercept_params
  true$link_params = link_params
  true$likelihood_params = likelihood_params
  
  # 0. Draw covarietes - assuming disease progression here
  observed$X = runif(n_tubes, 0, 1)
  observed$G = n_genes
  observed$T = n_tubes
  
  # 1. Draw from hyperprior
  recognized_hyperprior = FALSE
  if(hyperprior_params$type == "no") {
    recognized_hyperprior = TRUE
    true$beta1 = rnorm(n_genes, 0, hyperprior_params$sigma)
    
    observed$beta_sigma = array(hyperprior_params$sigma, 1)
    observed$hyperprior_type = 0
    
  } else {
    observed$beta_sigma = numeric(0)
  }
  
  if(hyperprior_params$type == "horseshoe") {
    recognized_hyperprior = TRUE
    true$horseshoe_sigma = abs(rnorm(1, 0, hyperprior_params$tau))
    true$beta1 = rnorm(n_genes, 0, true$horseshoe_sigma)
    
    observed$horseshoe_tau = array(hyperprior_params$tau, 1)
    observed$hyperprior_type = 1
    
  } else {
    observed$horseshoe_tau = numeric(0)
  }

  if(!recognized_hyperprior) {
    stop("Unrecognized hyperprior type: ",hyperprior_params$type)
  }
  
    
  # 2. Draw the intercepts
  true$intercept = rnorm(n_genes, intercept_params$mean, intercept_params$sigma)
  if(link_params$type == "sigmoid") { 
    #Interpet intercept as the location of the inflection point
    true$intercept_raw = true$intercept
    true$intercept = - true$intercept_raw * true$beta1
  }
  
  observed$intercept_mean = intercept_params$mean
  observed$intercept_sigma = intercept_params$sigma
  
  # 3. Linear predictor (X_beta in model)
  X_with_intercept = array(-Inf, c(n_tubes, 2))
  X_with_intercept[,1] = 1
  X_with_intercept[,2] = observed$X
  
  beta = array(-Inf, c(2, n_genes))
  beta[1,] = true$intercept
  beta[2,] = true$beta1
  
  true$linear_predictor = X_with_intercept %*% beta
  
  
  # 4. Link function
  recognized_link = FALSE
  if(link_params$type == "linear") {
    recognized_link = TRUE
    y_hat = true$linear_predictor
    observed$link_type = 0
  }    
  
  if(link_params$type == "exp") {
    recognized_link = TRUE
    y_hat = exp(true$linear_predictor)
    observed$link_type = 1
  } 
  
    if (link_params$type == "logsigmoid") {
    recognized_link = TRUE
    
    true$inversion = rnorm(n_genes, link_params$inversion_mean, link_params$inversion_sigma)
    y_hat <- matrix(-Inf, nrow = n_tubes, ncol = n_genes)
    for(g in 1:n_genes) {
      #TODO: there probably should be log(true$inversion[g]) AND inversion should be forced to be positive
      #TODO: the same in model code
      y_hat[, g] = true$inversion[g] + log1pexp(-true$intercept[g]) -  log1pexp(-true$linear_predictor[,g])
    }

    observed$link_type = 2    
    observed$inversion_mean = array(link_params$inversion_mean, 1)
    observed$inversion_sigma = array(link_params$inversion_sigma, 1)
  } else {
    observed$inversion_mean = numeric(0)
    observed$inversion_sigma = numeric(0)
  }

  if (link_params$type == "sigmoid") {
    recognized_link = TRUE
    #TODO check
    #true$mean_expression = rlnorm(n_genes, link_params$mean_expression_mean, link_params$mean_expression_sigma)
    true$mean_expression = exp(rnorm(n_genes) * link_params$mean_expression_sigma + link_params$mean_expression_mean)

    y_hat <- matrix(-Inf, nrow = n_tubes, ncol = n_genes)
    for(g in 1:n_genes) {
      sigmoid_out = 1 / (1 + exp(-true$linear_predictor[,g])  )
      plateau = true$mean_expression[g] / mean(sigmoid_out)
      y_hat[, g] = plateau * sigmoid_out
    }
    
    observed$link_type = 3    
    observed$mean_expression_mean = array(link_params$mean_expression_mean, 1)
    observed$mean_expression_sigma = array(link_params$mean_expression_sigma, 1)
  } else {
    observed$mean_expression_mean = numeric(0)
    observed$mean_expression_sigma = numeric(0)
  }  
  
    
  true$y = y_hat
  
  if(!recognized_link) {
    stop("Unrecognized link function:", link_params$type)
  }
  
  recognized_likelihood = FALSE
  if(likelihood_params$type  %in% c("normal","lognormal")) {
    true$xi = abs(rnorm(1, 0, likelihood_params$xi_sigma))

    if(likelihood_params$type == "normal") {
      recognized_likelihood = TRUE
      rfunc = rnorm
      observed$likelihood_type = 0
    } else if (likelihood_params$type == "lognormal") {
      recognized_likelihood = TRUE
      rfunc = rlnorm
      observed$likelihood_type = 1
    }
        
    observed$xi_sigma = likelihood_params$xi_sigma
    observed$y_real = rfunc(n_tubes * n_genes, y_hat, true$xi) %>% matrix(ncol = n_genes, nrow = n_tubes)
  } else {
    observed$y_real = array(0,c(0,0))
  }

  
  if(likelihood_params$type %in% c("neg_binomial","dirichlet_multinom")) {
    true$xi = abs(rnorm(1, 0, likelihood_params$xi_sigma))
    
    observed$xi_sigma = likelihood_params$xi_sigma
  } else {
    observed$y = array(as.integer(0), c(0,0))
  }
  
  if(likelihood_params$type == "neg_binomial") {
    recognized_likelihood = TRUE
    observed$y = rnbinom(n_tubes * n_genes, mu = y_hat, size = true$xi)  %>% matrix(ncol = n_genes, nrow = n_tubes)  
    observed$likelihood_type = 2
  } 
  
  if(likelihood_params$type == "dirichlet_multinom") {
    recognized_likelihood = TRUE
    
    observed$y = array(.Machine$integer.max, c(n_tubes, n_genes))
    
    #100 reads/gene on average
    observed$exposure = (50 * n_genes) + rgeom(n_tubes, 1 / (n_genes * 50) )
    for(t in 1:n_tubes) {
      #TODO: Add dependnecy on MCMCPack
      
      dirichlet_draw = MCMCpack::rdirichlet(1, softmax(y_hat[t,]))
      observed$y[t,] = rmultinom(1, observed$exposure[t], dirichlet_draw)
    }
    observed$likelihood_type = 3
  } else {
    observed$exposure = integer(0)
  }
  
  
  if(!recognized_likelihood) {
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


plot_simulated <- function(data, n_genes_to_show = min(data$observed$G, 9)) {
  genes_to_show = sample(1:data$observed$G, n_genes_to_show)
  
  if(data$true$likelihood_params$type %in% c("normal","lognormal")) {
    data_observed <- data$observed$y_real
  } else {
    data_observed <- data$observed$y
  }
  
  matrix_to_plot_data <- function(data, matrix_to_show, genes_to_show) {
    rownames(matrix_to_show) <- paste("Tube",1:data$observed$T)
    colnames(matrix_to_show) <- paste("Gene",1:data$observed$G)
    plot_data <- matrix_to_show[, genes_to_show] %>% as.data.frame() %>%
      rownames_to_column("tube") %>%
      mutate(X = data$observed$X) %>%
      gather("gene","expression", -X,-tube)
    
    plot_data
  }
  
  plot_data_observed <- matrix_to_plot_data(data, data_observed, genes_to_show)
  
  plot_data_true <- matrix_to_plot_data(data, data$true$y, genes_to_show)
  
  plot_data_observed %>% 
    ggplot(aes(x = X, y = expression)) + geom_point() + 
    geom_line(data = plot_data_true, color = "blue") +
    facet_wrap(~gene, scales = "free_y")
}