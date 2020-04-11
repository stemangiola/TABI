#' Get the sigmoid model (intended for easy transition between development
#' and release)
#'
#' @return The model
#'
get_sigmoid_model = function() {
  #stanmodels$DE_sigmoid ##This should be the variant for release
  rstan::stan_model(here::here("inst","stan","DE_sigmoid_try.stan"))
    }

#' Perform generalised linear model on RNA seq data
#'
#' @param X A design matrix
#' @param y A matrix of expression
#' @param model can override the default model for tests/development
#' @param house_keeping_genes A character string
#' @return A list
#'
#' @export
#'
sigmoid_link = function(
	X,
	y,
	multiplier,
	prior,
	iter,
	warmup,
	model = get_sigmoid_model(),
	control=list(
		adapt_delta=0.8,
		stepsize = 0.1,
		max_treedepth =10
	),
	prior_only = 0
){

	#######################################
	## Model 1
	#######################################

	# Create dimensions for model
	G = ncol(y)
	T = nrow(y)
	R_1 = ncol(X)-1

	# Exposure terms
	exposure = y %>% rowSums()

	# Horseshoe
	nu_local = 1
	nu_global = prior$nu_global
	par_ratio = prior$prop_DE
	slab_df = prior$slab_df
	slab_scale = prior$scale_DE
	
	# Run model
	fit =
		sampling(
		  model, #model
			iter =   iter,
			warmup = warmup,
			chains = 4,
			cores = 4,
			control = control
		)

	#######################################
	## Return
	#######################################

	list(

		# All fits
		fit = fit,

		# Produce output table posterior
		posterior_df = fit %>%
			gather_samples(beta[covariate_idx, gene_idx]) %>%
			left_join(
				tibble(
					gene_idx = 1:ncol(y),
					gene = colnames(y)
				)
				, by= "gene_idx"
			),

		# Generated quantities
		generated_quantities = fit %>%
			gather_samples(y_gen[sample_idx, gene_idx] ) %>%
			mean_qi()

	)

}


#' Simulate data
#'
#' @param delta_magnitude the magnitude of the slope
#' @param n_genes Number of genes to be simulated
#' @param changing_genes Number of changing genes
#' @return A list
#'
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom dplyr bind_cols
#' @importFrom dplyr mutate_if
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#'
#' @export
#'
simulate_from_sigmoid = function(delta_magnitude = 5, n_genes = 100, changing_genes = round(n_genes*0.3), hkg = round(n_genes*0.3), n_samples = 13, precision_NB =20){

	# Custom sigmoid
	inv_logit_gen = function(x, k)     k / ( exp( - x  ) + 1 )

	# Design matrix
	X = model.matrix( ~ sort(runif(n_samples, -1, 1)))

	# beta
	beta =
		data.frame(
			"(Intercept)" = rnorm(n_genes, 0, 1),
			"slope" = c(
				rnorm(hkg, 0, 0.001),
				rep(0, n_genes - changing_genes - hkg),
				rep( (delta_magnitude), changing_genes)
			)
		)

	# Max expression
	k = abs(rnorm(n_genes, 7, 2))

	# Calculate baselines
	y = t( X %*% t(beta) )
	y_sigmoid  = apply(y, 2, function(mc) inv_logit_gen(mc, k) )

	# Produce tissue samples
	y_real = foreach(n = 1:n_samples, .combine=bind_cols) %do% {
		tibble(    rnbinom(n_genes, mu = exp(y_sigmoid[,n]), size = precision_NB)    )
	} %>%
		setNames(sprintf("s%s", 1:n_samples)) %>%
		mutate_if(is.numeric, as.integer)

	# produce probability values
	y_prob = foreach(n = 1:n_samples, .combine=bind_cols) %do% { y_real[,n]/sum(y_real[,n]) } %>% as_tibble()

	# sample from multinomial
	y_observed = foreach(n = 1:n_samples, .combine=bind_cols) %do% {
		tibble(as.numeric(rmultinom(n = 1, prob = y_prob %>% pull(n), size = n_genes * 100)) )
	} %>%
		setNames(colnames(y_real)) %>%
		mutate_if(is.numeric, as.integer)

	list(
		G =                           nrow(y_observed),
		T =                           ncol(y_observed),
		F =                           n_genes - changing_genes,
		R =                           ncol(beta),
		y =                           t(y_observed),
		cond =                        cov,
		X =                           X,
		y_prob =                      y_prob,
		k =                           k,
		beta =                        beta,
		y_sigmoid =                   y_sigmoid,
		y_real =                      y_real
	)

}

