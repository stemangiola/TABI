library(tidyverse)
library(tidybayes)
library(TABI)
library(rstan)
source("/stornext/Home/data/allstaff/m/mangiola.s/PhD/TABI/dev/Article_Sections/TABI_Article_Functions.R", echo=TRUE)

sim_df <- read_table2("dev/sim_df.csv")
sim_df = sim_df %>% setNames(colnames(sim_df) %>% gsub('"', "", .))
sim_df = sim_df %>% filter(slope>0)
sim_df= sim_df %>% filter(slope>0) %>% nest(data = -Gene_number) %>% slice(2) %>% unnest(data) 


my_m = rstan::stan_model("inst/stan/DE_sigmoid_hierarchical.stan")


fits = 
  
  # Decide parameters
  expand_grid(beta = seq(-2, 2, length.out = 5), K = 4, A = 2, alpha = seq(-6, 6, length.out = 5))  %>% 
  mutate(samples = pmap(
    list(beta, K, A, alpha),
    ~ sigmoidal_sim_df(n_true_tests = 1, #single gene
                       n_false_tests = 0, #No null tests
                       beta = ..1,  
                       k = ..2,
                       A = ..3,
                       sample_size = 63,
                       disp_size = 1, 
                       covar =  seq(from = -5, to =5, by =0.5), # set of X labels to simulate over
                       alpha = rep(..4, 2) * -1, # alpha = 0,  same as inflection = 0 
                       seed_n  = 151) 
  )) %>% 
  
  # Fit
  mutate(fit = map(
    samples, 
    ~  TABI::TABI_glm(
      .data = .x,
      ~CAPRA_S,
      .sample = sample_id, 
      .transcript = Gene_number, 
      .abundance = value,
      model = my_m,
      iter = 1000,
      warmup = 800
    )
  )) %>% 
  
  # Plot
  mutate(X = map(fit, ~ .x$input.X %>%  rowid_to_column("T"))) %>%
  
  # Converged chains
  mutate(chains_converged = map(
    fit,
    ~ .x$posterior_df %>% distinct(.chain) %>% pull(.chain)
  )) %>% 
  
  
  mutate(y = pmap(
    list(fit,  X, chains_converged),
    ~ ..1$fit %>% 
      gather_draws(log_y_hat[T, G]) %>% 
      ungroup() %>% 
      # filter(.iteration > 390) %>%
      left_join(..2, by =  "T") %>% 
      filter(.chain %in% ..3)
  )) %>%
  
  mutate(
    p = map2(
      samples,
      y,
      ~ .x %>% ungroup() %>%  mutate(CAPRA_S = scale(CAPRA_S) %>% as.numeric) %>%
        ggplot(aes(CAPRA_S, log(value))) +
        geom_point() +
        geom_smooth(method = "lm", color = "red") +
        geom_line(data = .y, aes(
          CAPRA_S, .value, group = .draw, color = factor(.chain)
        ), alpha = 0.2)
    )
  )

saveRDS(fits, "dev/fits.rds")

for(i in 1:100){
  fits %>% slice(i) %>% pull(p) %>% .[[1]] %>% plot()
  print(i)
  fits %>% slice(i) %>% select(beta, K, A, alpha) %>% print()
  readline(prompt="Press [enter] to continue")
}

fits %>% slice(11) %>% pull(fit) %>% .[[1]] %$% fit %>% traceplot( inc_warmup=T)
