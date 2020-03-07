library(tidyverse)
library(TABI)
library(tidybulk)

TABI_TP <- 
  TABI::test_df %>%
  TABI_glm(
    ~ CAPRA_S,
    sample, transcript, count,
    model = rstan::stan_model("inst/stan/DE_sigmoid_hierarchical.stan")
  )
  