# library(tidyverse)
# library(TABI)
# library(tidybulk)
# library(magrittr)


test_that("test slope 0",{
  
  TABI_TP <- 
    TABI::test_df %>%
    mutate(count = rnbinom(nrow(TABI::test_df), 
                            mu=50, 
                                 size = 30) %>% 
                                            as.integer) %>%
    TABI::TABI_glm(
      .data = .,
      ~CAPRA_S,
      .sample = sample, 
      .transcript = transcript, 
      .abundance = count,
      #model = rstan::stan_model("inst/stan/DE_sigmoid_hierarchical.stan"),
      control=list(adapt_delta=0.99, 
                            stepsize = 0.01,  
                                  max_treedepth =10),
      iter = 2000,
      warmup = 1000
    )
  
  
  testthat::expect_equal(
    TABI_TP$fit %>% 
      rstan::extract("beta") %$% 
      beta %>% 
      as.numeric %>%
      {x = (.); (x > 0) %>% 
                         which %>% 
                            length %>% 
                                   magrittr::divide_by(length(x)) 
                                                         } , 
    0.5 ,
    tolerance=0.05
  )
  
  # TABI_TP$fit %>% pairs(pars=c("beta", "inflection", "A", "od", "y_cross"))
  
  
})


test_that("test slope 2",{
  
  sigmoid_4_param<- function(scaled_CAPRA_S,
                             A,
                             y_0,
                             slope,
                             inflection) { 
    
    return((y_0 - A)*(
                       1+exp(slope*inflection)
                                              )/(
                                                1+exp(-scaled_CAPRA_S %*% 
                                                                 slope+inflection*slope %>% as.vector()) 
                                              )
                                                 )
    
    
    }
  
  TABI_TP <- 
    TABI::test_df %>%
    # Simulate with slope
    mutate(count =   rnbinom(nrow(TABI::test_df), 
                               mu =exp(
                                   sigmoid_4_param(
                                          TABI::test_df$CAPRA_S %>% scale,
                                           A = 2,
                                               y_0 = 3,
                                                      slope = 2,
                                                            inflection= 1
                                                                          )
                                                                             ), 
                                                     size=30) %>% 
                                                         as.integer()) %>%
    TABI::TABI_glm(
      .data = .,
      ~ CAPRA_S,
      .sample = sample, 
      .transcript = transcript, 
      .abundance = count,
      #model = rstan::stan_model("inst/stan/DE_sigmoid_hierarchical.stan"),
      control=list(adapt_delta=0.99,
                              stepsize = 0.01,  
                                      max_treedepth =10),
      iter = 2000,
      warmup = 1000
    )
  
  testthat::expect_equal(
    TABI_TP$fit %>% 
      rstan::extract("inflection") %$% 
      inflection %>% 
      as.numeric %>% mean() , 
    1, 
    tolerance=0.5 
  )
  
  
})

# # Divergence plot
# bayesplot::mcmc_parcoord(
#   as.array(TABI_slope_0$fit, pars = c("beta", "inflection", "A", "od", "y_cross")),
#   np = bayesplot::nuts_params(TABI_slope_0$fit),
#   transform = function(x) {(x - mean(x)) / sd(x)}
# ) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


