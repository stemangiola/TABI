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
	house_keeping_genes
){

	#######################################
	## Model 1
	#######################################

	# Create dimensions for model
	G = ncol(y)
	T = nrow(y)
	F = length(house_keeping_genes %in% colnames(y))
	R = ncol(X)

	# Exposure terms
	exposure = y %>% rowSums()

	# Run model
	fit1 =
		sampling(
			stanmodels$DE_sigmoid,
			#rstan::stan(file = "~/PhD/TABI/src/stan_files/DE_sigmoid.stan",
			iter =   500,
			warmup = 250,
			chains = 4,
			cores = 4
			#, control=list(adapt_delta=0.95, stepsize = 0.05, max_treedepth =15)
		)

	# Produce output table
	CI_df = fit1 %>%
		gather_samples(beta[covariate, gene]) %>%
		filter(covariate==2) %>%
		mean_qi(.prob =  1 - (0.05 / G) )

	#######################################
	## Return
	#######################################

	list(

		# All fits
		fit1 = fit1,

		# Table of confidence intervals
		CI_df = CI_df,

		# Plot of coefficient
		p =
			CI_df %>%
			mutate(sig = !(conf.low <= 0 & conf.high >=0 )) %>%
			filter(gene > F) %>%
			ggplot(aes(x=gene, y=estimate, color=sig)) +
			geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width = 0)

	)

}

