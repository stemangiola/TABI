

#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

# Script Intention


# (Section 1)
# Simulate sigmoidal based RNAseq datasets with varying sample sizes and varying slope sizes

# (Section 2)
# Use TABI to analyse them 

# (Section 3)
# Analyse them using common differential expression analysis software
# edgeR, DESeq2, Bayseq 

# (Section 4)
# Compare the mean error of TABI to edgeR, DESeq2, Bayseq



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


# (Section 1)

# Sample size = 21
# 
test_df_1= lapply(c(-1,1), function(x)
sigmoidal_sim(n_true_tests = 50,
              n_false_tests = 0,
              beta = x,
              k = 5,
              A = 3,
              sample_size = 21,
              disp_size = 1,
              x_cord = seq(from = -5, to = 5, by =0.5),
              alpha = "all") )

saveRDS(test_df_1, "test_df_1.rds")


library(furrr)

test_df_1 %>% 
  bind_rows() %>% 
  #mutate(sample_id = as.factor(sample_id)) %>% 
  #mutate(CAPRA_S = as.character(CAPRA_S)) %>% 
  split(.$Gene_ref) %>% 
  .[[1]] %>% 
  TABI::TABI_glm(
    ~CAPRA_S,
    .sample = sample_id,
    .transcript = Gene_ref,
    .abundance = value,
    chains = 4,
    cores = 1,
    # control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
    iter = 900,
    warmup = 300
  )



# Sample size = 42


# Sample size = 84




              
              
              
              
