#' The 'TABI' package.
#'
#' @description A DESCRIPTION OF THE PACKAGE
#'
#' @docType package
#' @name TABI-package
#' @aliases TABI
#' @useDynLib TABI, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr everything
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr distinct
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#' @importFrom dplyr enquo
#'
#' @importFrom purrr is_numeric
#'
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#'
#' @importFrom tidybayes gather_samples
#' @importFrom tidybayes mean_qi
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#'
#'
#'
#'
NULL



#' Perform generalised linear model on RNA seq data
#'
#' @param formula A formula
#' @param .data A tibble
#' @param link A character string
#' @param .sample A column name
#' @param .abundance A column name
#' @param  .transcript A column name
#' 
#' @importFrom tidybulk tidybulk
#' @importFrom dplyr mutate_if
#' @importFrom dplyr enquo
#' @importFrom tidyr pivot_wider
#' @importFrom tidybulk scale_abundance
#' @importFrom tidybulk aggregate_duplicates
#' 
#' @return A tibble
#' 
#' @export
#'
TABI_glm = function(
  .data,
  .formula,
	.sample, # Sample ID column
	.transcript, # Transcript name/ID column
	.abundance, # Raw expression read column
	link = "sigmoid",
	prop_DE =0.05,
	scale_DE = 5,
	nu_global = 5,
	slab_df = 40	,
	iter = 500,
	warmup = round(iter/2),
	model = get_sigmoid_model(),
	control=list(
		adapt_delta=0.8,
		stepsize = 0.1,
		max_treedepth =10
	),
	prior_only = 0,
  shards = 1
){
 
  Sys.setenv("STAN_NUM_THREADS" = shards)
   
	# Set prior
	prior = list(
		prop_DE =prop_DE,
		scale_DE = scale_DE,
		nu_global = nu_global,
		slab_df = slab_df
	)

	# Make col names
	.sample = dplyr::enquo(.sample)
	.transcript = dplyr::enquo(.transcript)
	.abundance = dplyr::enquo(.abundance)
	
	# Covariate column
	cov_columns =
	  parse_formula(.formula)$covariates %>%
	  map_chr(~ .x %>% gsub("censored\\(|\\)| ", "", .) %>% str_split("\\,") %>% `[[` (1) %>% `[` (1)) %>%
	  ifelse_pipe((.) %>% is.null, ~ c())
	
	# Censoring column
	.cens_label_column = parse_formula(.formula)$covariates %>% grep("censored(", ., fixed = T, value = T)  %>% gsub("censored\\(|\\)| ", "", .) %>% str_split("\\,") %>% ifelse_pipe(length(.)>0, ~.x %>% `[[` (1) %>% `[` (-1), ~NULL)
	.cens_value_column = parse_formula(.formula)$covariates %>% grep("censored(", ., fixed = T, value = T)  %>% gsub("censored\\(|\\)| ", "", .) %>% str_split("\\,") %>% ifelse_pipe(length(.)>0, ~.x %>% `[[` (1) %>% `[` (1), ~NULL)
	
	if(length(.cens_label_column) == 1) {
	  cens = .data %>% select(sample, .cens_label_column) %>% distinct %>% arrange(sample) %>% pull(2)
	  
	  sd_survival_months = 29.3
	  
	  .data = .data %>% mutate(!!.cens_value_column := !!as.symbol(.cens_value_column) / sd_survival_months)
	  my_formula = as.formula( paste("~",paste(cov_columns, collapse = "+")))
	}
	else{
	  cens = NULL
	  my_formula = .formula
	} 
	
	# Create design matrix
	X =
	  model.matrix(object = my_formula, 
	               data = .data %>%  select(one_of(cov_columns), !!.sample) %>% distinct()) %>%
	  as_tibble(rownames="sample_idx") 
	#%>%
	#  scale_design(.formula)

	# tt Object to run through tidybulk
	multiplier =
	  .data %>%
	  tidybulk::tidybulk( !!.sample, !!.transcript, !!.abundance) %>%
	   tidybulk::aggregate_duplicates() %>%
	  tidybulk::scale_abundance(action="get") %>%
	  arrange(!!.sample) %>%
	  select(multiplier)
	
	# Create tibble of expression data
	# Each column is a different gene
	y = 
	  .data %>%
	  select( !!.sample, !!.transcript, !!.abundance) %>%
	  spread(!!.transcript, !!.abundance) %>%
	  select(-!!.sample)

	# Return
	c(

		# Return the inputs to the model
		input = list(
			X = X,
			y = y,
			.data = .data
		),

		# Return the outcome of the model
		switch(
			link,
			"sigmoid" =
				sigmoid_link(
					X, y, multiplier,
					prior,
					iter,
					warmup,
					model = model,
					control = control,
					prior_only = prior_only,
					shards = shards,
					cens
				)
		)
	)
}

#' Plots the observed aganst the modelled read counts for a gene
#'
#' @param TABi_obj A TABI output object
#' @param gene A character string
#' @param covariate A character string
#' @param CI A real number
#'
#' @return A tibble
#'
#' @export
#'
plot_generated_gene = function(TABi_obj, gene, covariate = colnames(TABi_obj$input.X)[2], CI = 1 - (0.05 / ncol(TABi_obj$input.y))){

	TABi_obj$generated_quantities %>%

		# Add gene names
		left_join(
			tibble(
				gene_idx = 1:ncol(TABi_obj$input.y),
				gene =     colnames(TABi_obj$input.y)
			),
			by = "gene_idx"
		) %>%

		# Fiter gene
		filter(gene == !!gene) %>%

		# Attach X
		bind_cols(
			TABi_obj$input.X %>%
				select(!!covariate) %>%
				setNames("x")
		) %>%

		# Attach observed y
		bind_cols(
			TABi_obj$input.y %>% select(!!gene) %>%
				setNames("observed")
		) %>%

		# Plot
		ggplot(aes(x=x, y=observed)) +
		geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width = 0) +
		geom_point(color="red") +
		theme_bw()



}

#' Plots the posterior credible interval across all genes
#'
#' @param TABi_obj A TABI output object
#' @param covariate A character string
#' @param CI A real number
#'
#' @return A tibble
#'
#' @export
#'
plot_posterior = function(TABi_obj, covariate = colnames(TABi_obj$input.X)[2], CI = 1 - (0.05 / ncol(TABi_obj$input.y))){

	TABi_obj$posterior_df %>%

		# Filter covariate
		filter(covariate_idx==which(colnames(TABi_obj$input.X) == !!covariate)) %>%

		# Get statistics, Exclude gene for tidybayes bug
		select(-gene) %>%
		mean_qi( ) %>%

		# Add gene names
		ungroup() %>%
		left_join(
			TABi_obj$posterior_df %>% ungroup() %>% distinct(gene_idx, gene),
			by = "gene_idx"
		) %>%

		# Attribute significance
		mutate(sig = !(.lower <= 0 & .upper >=0 )) %>%

		# Plot
		ggplot(aes(x=gene, y=estimate, color=sig)) +
		geom_errorbar(aes(ymin=.lower, ymax=.upper), width = 0) +
		theme_bw()

}




