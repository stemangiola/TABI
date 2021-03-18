

# Intention of Script 
library(dplyr)

{
  
 #C12orf49  

C12orf49_TABI= ed_PC_TCGA %>% 
  select(CAPRA_S = `CAPRA-S`, read_count_normalised = `read count normalised`, transcript, purity.score) %>% 
  filter(transcript == "C12orf49") %>% 
  mutate(read_count_normalised = as.integer(read_count_normalised)) %>% 
  #mutate(Sample = as.character(Sample)) %>% 
  filter(is.na(CAPRA_S) == F) %>% 
  arrange(CAPRA_S) %>% 
  mutate(sample = c(1:135)) %>% 
  TABI::TABI_glm(
    ~ CAPRA_S+purity.score,
    .sample = sample, .transcript = transcript, .abundance = read_count_normalised,
    control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
    iter = 5000,
    warmup = 2500,
  ) 

saveRDS(
  C12orf49_TABI,
  compress = "gzip", 
  file = "C12orf49_TABI.rds"
)


}


{
  # IQGAP3
IQGAP3_TABI= ed_PC_TCGA %>% 
  select(CAPRA_S = `CAPRA-S`, read_count_normalised = `read count normalised`, Sample = Sample.ID, transcript, purity.score) %>% 
  filter(transcript == "IQGAP3") %>% 
  mutate(read_count_normalised = as.integer(read_count_normalised)) %>% 
  #mutate(Sample = as.character(Sample)) %>% 
  filter(is.na(CAPRA_S) == F) %>% 
  arrange(CAPRA_S) %>% 
  mutate(sample = c(1:135)) %>% 
  TABI::TABI_glm(
    ~ CAPRA_S+purity.score,
    .sample = sample, .transcript = transcript, .abundance = read_count_normalised,
    control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
    iter = 5000,
    warmup = 2500,
  ) 


saveRDS(
  IQGAP3_TABI,
  compress = "gzip", 
  file = "IQGAP3_TABI.rds"
)

}


{
# CXCR5

CXCR5_TABI= ed_PC_TCGA %>% 
  select(CAPRA_S = `CAPRA-S`, read_count_normalised = `read count normalised`, Sample = Sample.ID, transcript, purity.score) %>% 
  filter(transcript == "CXCR5") %>% 
  mutate(read_count_normalised = as.integer(read_count_normalised)) %>% 
  #mutate(Sample = as.character(Sample)) %>% 
  filter(is.na(CAPRA_S) == F) %>% 
  arrange(CAPRA_S) %>% 
  mutate(sample = c(1:135)) %>% 
  TABI::TABI_glm(
    ~ CAPRA_S+purity.score,
    .sample = sample, .transcript = transcript, .abundance = read_count_normalised,
    control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
    iter = 5000,
    warmup = 2500,
  ) 


saveRDS(
  CXCR5_TABI,
  compress = "gzip", 
  file = "CXCR5_TABI.rds"
)

}



#PROK1

{

PROK1_TABI= ed_PC_TCGA %>% 
  select(CAPRA_S = `CAPRA-S`, read_count_normalised = `read count normalised`, Sample = Sample.ID, transcript, purity.score) %>% 
  filter(transcript == "PROK1") %>% 
  mutate(read_count_normalised = as.integer(read_count_normalised)) %>% 
  #mutate(Sample = as.character(Sample)) %>% 
  filter(is.na(CAPRA_S) == F) %>% 
  arrange(CAPRA_S) %>% 
  mutate(sample = c(1:135)) %>% 
  TABI::TABI_glm(
    ~ CAPRA_S+purity.score,
    .sample = sample, .transcript = transcript, .abundance = read_count_normalised,
    control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
    iter = 5000,
    warmup = 2500,
  ) 


saveRDS(
  PROK1_TABI,
  compress = "gzip", 
  file = "PROK1_TABI.rds"
)
}


# LGI3

{
  
  LGI3_TABI= ed_PC_TCGA %>% 
    select(CAPRA_S = `CAPRA-S`, read_count_normalised = `read count normalised`, Sample = Sample.ID, transcript, purity.score) %>% 
    filter(transcript == "LGI3") %>% 
    mutate(read_count_normalised = as.integer(read_count_normalised)) %>% 
    #mutate(Sample = as.character(Sample)) %>% 
    filter(is.na(CAPRA_S) == F) %>% 
    arrange(CAPRA_S) %>% 
    mutate(sample = c(1:135)) %>% 
    TABI::TABI_glm(
      ~ CAPRA_S+purity.score,
      .sample = sample, .transcript = transcript, .abundance = read_count_normalised,
      control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
      iter = 5000,
      warmup = 2500,
    ) 
  
  
  saveRDS(
    LGI3_TABI,
    compress = "gzip", 
    file = "LGI3_TABI.rds"
  )
  
  
}



# SPRY4

{
  
  SPRY4_TABI= ed_PC_TCGA %>% 
    select(CAPRA_S = `CAPRA-S`, read_count_normalised = `read count normalised`, Sample = Sample.ID, transcript, purity.score) %>% 
    filter(transcript == "SPRY4") %>% 
    mutate(read_count_normalised = as.integer(read_count_normalised)) %>% 
    #mutate(Sample = as.character(Sample)) %>% 
    filter(is.na(CAPRA_S) == F) %>% 
    arrange(CAPRA_S) %>% 
    mutate(sample = c(1:135)) %>% 
    TABI::TABI_glm(
      ~ CAPRA_S+purity.score,
      .sample = sample, .transcript = transcript, .abundance = read_count_normalised,
      control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
      iter = 5000,
      warmup = 2500,
    ) 
  
  
  saveRDS(
    SPRY4_TABI,
    compress = "gzip", 
    file = "SPRY4_TABI.rds"
  )
  
  
}




# SLC15A1

{
  
  SLC15A1_TABI= ed_PC_TCGA %>% 
    select(CAPRA_S = `CAPRA-S`, read_count_normalised = `read count normalised`, Sample = Sample.ID, transcript, purity.score) %>% 
    filter(transcript == "SLC15A1") %>% 
    mutate(read_count_normalised = as.integer(read_count_normalised)) %>% 
    #mutate(Sample = as.character(Sample)) %>% 
    filter(is.na(CAPRA_S) == F) %>% 
    arrange(CAPRA_S) %>% 
    mutate(sample = c(1:135)) %>% 
    TABI::TABI_glm(
      ~ CAPRA_S+purity.score,
      .sample = sample, .transcript = transcript, .abundance = read_count_normalised,
      control=list( adapt_delta=0.9,stepsize = 0.01,  max_treedepth =10  ),
      iter = 5000,
      warmup = 2500,
    ) 
  
  
  saveRDS(
    SLC15A1_TABI,
    compress = "gzip", 
    file = "SLC15A1_TABI.rds"
  )
  
  
}







