#' Formula parser
#'
#' @param fm A formula
#'
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
  pars = as.character(attr(terms(fm), "variables"))[-1]
  
  response = NULL
  if(attr(terms(fm), "response") == 1) response = pars[1]
  covariates = ifelse(attr(terms(fm), "response") == 1, pars[-1], pars)
  
  list(
    response = response,
    covariates = covariates
  )
}

#' Scale design matrix
#'
#' @param df A tibble
#' @return A tibble
#'
#'
scale_design = function(df, formula){
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

#' make a function to combine the results
stan.combine <- function(...) { return( sflist2stanfit( list(...) )  ) }

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import dplyr
#' @import tidyr
#'
#' @param input.df A tibble
#' @param condition A boolean
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
  switch(.p %>% `!` %>% sum(1),
         as_mapper(.f1)(.x),
         if (.f2 %>% is.null %>% `!`)
           as_mapper(.f2)(.x)
         else
           .x)
  
}

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import dplyr
#' @import tidyr
#'
#' @param .x A tibble
#' @param .p1 A boolean
#' @param .p2 ELSE IF condition
#' @param .f1 A function
#' @param .f2 A function
#' @param .f3 A function
#'
#' @return A tibble
ifelse2_pipe = function(.x, .p1, .p2, .f1, .f2, .f3 = NULL) {
  # Nested switch
  switch(# First condition
    .p1 %>% `!` %>% sum(1),
    
    # First outcome
    as_mapper(.f1)(.x),
    switch(
      # Second condition
      .p2 %>% `!` %>% sum(1),
      
      # Second outcome
      as_mapper(.f2)(.x),
      
      # Third outcome - if there is not .f3 just return the original data frame
      if (.f3 %>% is.null %>% `!`)
        as_mapper(.f3)(.x)
      else
        .x
    ))
}