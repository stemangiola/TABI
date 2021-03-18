

# Run TABI on simulated dataset with 100 genes, positive slope (1), 100 genes negative slope (-1)


library(furr)
library(TABI)
library(dplyr)
library(tictoc)


tic()
y<-TABI_df %>% 
  #filter(Gene_ref == "X6-1") %>% 
  rename(sample = Sample) %>% 
  split(.$Gene_ref) %>% # %$%
 # .[1:3]%>% 
future_map_dfr(~TABI::TABI_glm(
     .data = .x,
    ~ CAPRA_S,
    .sample = sample,
    .transcript = Gene_ref,
    chains = 4,
    cores = 4,
    .abundance = value,
    #control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
    iter = 1500,
    warmup = 300
  ) %$% 
  fit %>% 
   rstan::summary() %$%
   summary %>% 
   as.data.frame() %>% 
    tibble::rownames_to_column("Pars") %>% 
    filter(Pars == "inflection[1]")
  ) 
toc()
