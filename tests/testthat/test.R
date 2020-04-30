library(tidyverse)
library(TABI)
library(tidybulk)
library(magrittr)

load("data/test_df.rda")

test_that("test slope 0",{
  
  TABI_TP <- 
    test_df %>%
    mutate(count = rnbinom(n(), mu=50, size = 30) %>% as.integer) %>%
    TABI_glm(
      ~ CAPRA_S,
      sample, transcript, count,
      #model = rstan::stan_model("inst/stan/DE_sigmoid_hierarchical.stan"),
      control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
      iter = 2000,
      warmup = 1000,
    )
  
  
  expect_gt(
    TABI_TP$fit %>% 
      rstan::extract("beta") %$% 
      beta %>% 
      as.numeric %>%
      shapiro.test %>%
      broom::tidy() %>%
      pull(p.value) , 
    0.1 
  )
  
  #TABI_TP$fit %>% pairs(pars=c("beta", "inflection", "A", "od", "y_cross"))
  
  
})


test_that("test slope 0",{
  
  
  
  TABI_TP <- 
    test_df %>%
    
    # Simulate with slope
    mutate(count =   rnbinom(nrow(test_df), mu =exp(sigmoid_4_param(
      test_df$CAPRA_S %>% scale,
      A = 2,
      y_cross = 0.5,
      slope = matrix(2, nrow = 1),
      inflection= 1
    )), size=30)) %>%
    mutate(CAPRA_S =  CAPRA_S %>% scale) %>%
    
    # Run model
    TABI_glm(
      ~ CAPRA_S,
      sample, transcript, count,
      #model = rstan::stan_model("inst/stan/DE_sigmoid_hierarchical.stan"),
      control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10 ),
      iter = 2000,
      warmup = 1000,
    )
  
  
  expect_equal(
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


