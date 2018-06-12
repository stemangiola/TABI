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
	nu_local = array(rep(1,R_1), dim=R_1)
	nu_global = array(rep(45, R_1), dim=R_1)
	par_ratio = array(prior$prop_DE, dim=R_1)
	slab_df = array(rep(4, R_1), dim=R_1)
	slab_scale = array(prior$scale_DE, dim=R_1)

	# Run model
	fit =
		#sampling(
			#stanmodels$DE_sigmoid,
			rstan::stan(file = "~/PhD/TABI/src/stan_files/DE_sigmoid.stan",
			iter =   iter,
			warmup = warmup,
			chains = 4,
			cores = 4
			#, control=list(adapt_delta=0.95, stepsize = 0.05, max_treedepth =15)
		)

	# Produce output table
	CI_df = fit %>%
		gather_samples(beta[covariate, gene]) %>%
		filter(covariate==2) %>%
		mean_qi(.prob =  1 - (0.05 / G) )

	#######################################
	## Return
	#######################################

	list(

		# All fits
		fit = fit,

		# Produce output table CI
		CI_df = CI_df,

		# Generated quantities
		generated_quantities = fit %>%
			gather_samples(y_gen[sample_idx, gene_idx] ) %>%
			mean_qi(),

		# Plot of coefficient
		p =
			CI_df %>%
			mutate(sig = !(conf.low <= 0 & conf.high >=0 )) %>%
			filter(gene > F) %>%
			ggplot(aes(x=gene, y=estimate, color=sig)) +
			geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width = 0)

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
simulate_from_sigmoid = function(delta_magnitude = 5, n_genes = 100, changing_genes = round(n_genes*0.3), hkg = round(n_genes*0.3)){

	# Custom sigmoid
	inv_logit_gen = function(x, k)     k / ( exp( - x  ) + 1 )

	# number of changing genes
	#changing_genes = round(n_genes * prop_canging)

	# numbe rof samples
	n_samples = 13

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
		tibble(    rnbinom(n_genes, mu = exp(y_sigmoid[,n]), size = 100)    )
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
