




Fig_2_sim_data_sample_105 <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Presentation_30_July/ROC_FDR_plots/Fig_2_sim_data_sample_105.rds")



{
  
  set.seed(40)
  Fig_2_sim_data_sample_21<-inner_join(
    Fig_2_sim_data_sample_105,
    Fig_2_sim_data_sample_105 %>% 
      dplyr::select(CAPRA_S, sample) %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(CAPRA_S) %>% 
      sample_n(1))  
  
} 

library(tidyverse)
library(TABI)
library(tictoc)
tic()
x<-Fig_2_sim_data_sample_21 %>% 
  filter(Gene_ref == "X64_Beta_-1_k_3") %>% 
TABI_glm(
  ~ CAPRA_S,
  .sample = sample,
  .transcript = Gene_ref,
  chains = 3,
  cores = 1,
  .abundance = value,
 # control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
  iter = 900,
  warmup = 300,
)
toc()

tic()
Fig_2_sim_data_sample_21 %>% 
  filter(Gene_ref == "X1_Beta_-1_k_3") %>% 
  TABI_glm(
    ~ CAPRA_S,
    .sample = sample,
    .transcript = Gene_ref,
    chains = 4,
    cores = 4,
  .abundance = value,
  control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
  iter = 2000,
  warmup = 1000,
)
toc()


plan(multisession, workers = 5)

library(furrr)

x<-((Fig_2_sim_data_sample_21 %>% 
  split(.$Gene_ref))) %>% 
  .[1:50] %>% 
  future_map(~TABI::TABI_glm(
    .data = .x,
    ~ CAPRA_S,
    .sample = sample,
    .transcript = Gene_ref,
    chains = 3,
    cores = 1,
    .abundance = value,
   # control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
    iter = 900,
    warmup = 300,
  ) %$% 
    posterior_df, 
  .progress = TRUE) %>% 
  future_map_dfr(~as.data.frame(.))
