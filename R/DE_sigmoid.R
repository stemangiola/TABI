#' Perform generalised linear model on RNA seq data
#'
#' @param X A design matrix
#' @param y A matrix of expression
#' @param house_keeping_genes A character string
#' @return A list
#'
#'
sigmoid_link = function(
	X,
	y,
	prior,
	iter,
	warmup
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
	nu_global = 45
	par_ratio = prior$prop_DE
	slab_df = 4
	slab_scale = prior$scale_DE

	# Set inits
	init.fn <- function(chain) list(xi_z=runif(1, 1, 2))

	# Run model
	fit =
		sampling(
			stanmodels$DE_sigmoid,
			#rstan::stan(file = "~/PhD/TABI/src/stan_files/DE_sigmoid.stan",
			iter =   iter,
			warmup = warmup,
			chains = 4,
			cores = 4,
			init = init.fn
			#, control=list(adapt_delta=0.95, stepsize = 0.05, max_treedepth =15)
		)

	#######################################
	## Return
	#######################################

	list(

		# All fits
		fit = fit,

		# Produce output table posterior
		posterior_df = fit %>% gather_samples(beta[covariate_idx, gene_idx]),

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
