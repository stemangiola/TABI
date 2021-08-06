# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

#' Get the sigmoid model (intended for easy transition between development
#' and release)
#'
#' @return The model
#'
get_sigmoid_model = function() {
  #stanmodels$DE_sigmoid ##This should be the variant for release
  #rstan::stan_model(here::here("inst","stan","DE_sigmoid_norm_factor_log_space_fit_vert_trans.stan"))
  rstan::stan_model("./inst/stan/DE_sigmoid_hierarchical.stan", auto_write = T)
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
	model = stanmodels$DE_sigmoid_hierarchical,
	control=list(
		adapt_delta=0.8,
		stepsize = 0.1,
		max_treedepth =10
	),
	prior_only = 0,
	shards = 1,
	cores = cores,
	chains = chains
){

	#------------------------------#
	# Model 1
	#------------------------------#

	# Create dimensions for model
	G = ncol(y)
	T = nrow(y)
	R_1 = ncol(X)-1

	# Run model
	fit =
		sampling(
		  model,
			iter =   iter,
			warmup = warmup,
			chains = chains,
			cores = cores,
			control = control,
			init = init_fun <- function(...) list(inflection=as.array(c(0)))
		)

	#------------------------------#
	# Return
	#------------------------------#

	list(

		# All fits
		fit = fit,

		# Produce output table posterior
		posterior_df = 
		  fit %>%
		  tidybayes::gather_draws(beta[covariate_idx, gene_idx]) %>%
		  
		  # Drop bimodal
		  ungroup() %>% 
		  # parse_summary_check_divergence() %>% 
		  cluster_posterior_slopes() %>% 
		  
			left_join(
				tibble(
					gene_idx = 1:ncol(y),
					gene = colnames(y)
				), 
				by= "gene_idx"
			),

		# Generated quantities
		generated_quantities = 
		  fit %>%
			tidybayes::gather_draws(y_gen[sample_idx, gene_idx] ) %>%
			mean_qi()

	)
}

#@description Get which chain cluster is more opulated in case I have divergence
choose_chains_majority_roule = function(fit_parsed) {
  
  my_chains = 
    # Calculate modes
    fit_parsed %>%
    select(.chain, .value) %>%
    {
      main_cluster =
        (.) %>%
        pull(.value) %>%
        kmeans(centers = 2, nstart = 10) %>% {
          bind_cols(
            (.) %$% centers %>% as_tibble(.name_repair = "minimal") %>% setNames("center") ,
            (.) %$% size %>% as_tibble(.name_repair = "minimal")  %>% setNames("size")
          )
        } %>%
        arrange(size %>% desc) %>%
        slice(1)
      
      (.) %>%
        group_by(.chain) %>%
        summarise(
          .lower_chain = quantile(.value, probs = c(0.025)),
          .upper_chain = quantile(.value, probs = c(0.975))
        ) %>%
        ungroup %>%
        mutate(center = main_cluster %>% pull(center))
    } %>%
    
    # Filter cains
    rowwise() %>%
    filter(between(center, .lower_chain, .upper_chain)) %>%
    ungroup %>%
    distinct(.chain)
  
  if(nrow(my_chains)==0) my_chains = fit_parsed %>% distinct(.chain)
  
  fit_parsed %>%
    inner_join(my_chains ,  by = ".chain"  )
}

# @description Parse the stan fit object and check for divergences
parse_summary_check_divergence = function(draws) {
  draws %>%
    
    # If not converged choose the majority chains
    mutate(converged = diptest::dip.test(.value) %$%	`p.value` > 0.05) %>%
    
    # If some proportions have not converged chose the most populated one
    when(
      (.) %>% distinct(converged) %>% pull(1) %>% `!` ~ choose_chains_majority_roule(.),
      ~ (.)
    ) 
}

#' Get matrix from tibble
#'
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#' library(dplyr)
#'
#' tidybulk::se_mini %>% tidybulk() %>% select(feature, count) %>% head %>% as_matrix(rownames=feature)
#'
#' @export
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {
  rownames = enquo(rownames)
  tbl %>%
    
    # Through warning if data frame is not numerical beside the rownames column (if present)
    when(
      do_check &&
        tbl %>%
        # If rownames defined eliminate it from the data frame
        when(!quo_is_null(rownames) ~ (.)[,-1], ~ (.)) %>%
        dplyr::summarise_all(class) %>%
        tidyr::gather(variable, class) %>%
        pull(class) %>%
        unique() %>%
        `%in%`(c("numeric", "integer")) %>% not() %>% any()   ~ {
        warning("tidybulk says: there are NON-numerical columns, the matrix will NOT be numerical")
        (.)
      },
      ~ (.)
    ) %>%
    as.data.frame() %>%
    
    # Deal with rownames column if present
    when(
      !quo_is_null(rownames) ~ (.) %>%
        magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
        select(-1),
      ~ (.)
    ) %>%
    
    # Convert to matrix
    as.matrix()
}

#' @importFrom tidygraph tbl_graph
#' @keywords internal
#' 
#' @param .data A tibble
#' @param credible_interval A double
cluster_posterior_slopes = function(draws){
  
  my_chains = 
    draws %>% 
    group_by(covariate_idx, gene_idx, .chain, .variable) %>% 
    tidybayes::mean_qi() %>% 
    ungroup() %>% 
    select(.chain, .lower, .upper, .value) %>% 
    mutate_if(is.integer, as.character) %>% 
    nanny::combine_nest(
      .names_from = .chain,
      .values_from = c(.lower, .upper, .value)
    ) %>% 
    mutate(is_cluster = map_lgl(
      data,
      ~ {
        
          .x[1,]$.value %>% between(.x[2,]$.lower, .x[2,]$.upper) &
          .x[2,]$.value %>% between(.x[1,]$.lower, .x[1,]$.upper)
        
        }

    )) %>% 
    # Find communities based on cell type clusters
    {
      ct_levels = (.) %>% arrange(!is_cluster) %>% select(1:2) %>% as_matrix %>% t %>% as.character() %>% unique
      
      (.) %>%
        filter(is_cluster) %>% 
        select(1:2) %>%
        tbl_graph(
          edges = .,
          nodes = data.frame(name = ct_levels)
        )%>%
        mutate(community = as.factor(tidygraph::group_infomap())) 
    } %>% 
    suppressWarnings() %>% 
    as_tibble() %>% 
    nest(data = -community) %>% 
    mutate(n = map_int(data, ~ nrow(.x))) %>% 
    filter(n == max(n)) %>% 
    unnest(data) %>% 
    pull(name) %>% 
    as.integer()
  
  draws %>% 
    filter(.chain %in% my_chains)
  
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

