Test_TP <-
  data.frame(CAPRA_S = seq(from = 0, to = 8, by = 0.1),
             y_sim_gene = y_sim)

TABI_TP <- 
  TABI::test_df %>%
  TABI_glm(
    formula = ~ CAPRA_S,
    sample, transcript, count,
    model = rstan::stan_model("init/stan/DE_sigmoid_hierarchical.stann")
  )
  