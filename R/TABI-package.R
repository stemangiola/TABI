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
#' @importFrom dplyr pull
#' @importFrom dplyr pull
#' @importFrom dplyr pull
#'
#' @importFrom tidyr spread
#'
#' @importFrom tidybayes gather_samples
#' @importFrom tidybayes mean_qi
#'
#'
#'
#'
NULL

#' Perform generalised linear model on RNA seq data
#'
#' @param formula A formula
#' @param data A tibble
#' @param link A character string
#' @return A tibble
#'
#' @export
#'
TABI_glm = function(
	formula,
	data,
	link = "sigmoid",
	house_keeping_genes = house_keeping_genes %>% pull(symbol)
){

	# Formula parser
	parse_formula <- function(fm) {
		if(attr(terms(fm), "response") == 1) stop("The formula must be of the kind \"~ covariates\" ")
		else as.character(attr(terms(fm), "variables"))[-1]
	}

	# Scale design matrix
	scale_design = function(df){
		df %>%
			setNames(c("sample_idx", "(Intercept)", parse_formula(formula))) %>%
			gather(cov, value, -sample_idx) %>%
			group_by(cov) %>%
			mutate( value = ifelse( !grepl("Intercept", cov) & length(union(c(0,1), value)) != 2, scale(value), value )) %>%
			ungroup() %>%
			spread(cov, value) %>%
			arrange(as.integer(sample_idx)) %>%
			select(`(Intercept)`, one_of(parse_formula(formula)))
	}

	# Create design matrix
	X =
		model.matrix(object = formula, data = data) %>%
		as_tibble(rownames="sample_idx") %>%
		scale_design()

	# Set up expression data frame
	y =
		data %>%
		mutate_if(is_numeric, as.integer) %>%
		select(-one_of(parse_formula(formula))) %>%
		select(one_of(house_keeping_genes), everything())

	switch(
		link,
		"sigmoid" = sigmoid_link(X, y, house_keeping_genes)
	)
}
