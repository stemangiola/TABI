#Parallel TABI using purity as a covariate 
#Getting fit for all values 

library(TABI)
library(tibble)
library(rstan)
library(dplyr)

#load compressed  RDS file for CAPRA_S moved from 8 to 7, and also with na removed
normalised_PC_TCGA <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/compress_normalised_narm_TCGA.rds")


#Function for saving fitting table 
#Creates tibble which can be combined in parallel
fit_purity<- function(gene_number) {
  #Name Gene to Test
  gene<- levels(normalised_PC_TCGA$transcript)[gene_number]
  
  
  #Required Packages
  
  require(TABI)
  require(dplyr)
  require(rstan)
  require(tibble)
  
  #Only evaluate if there is at least one non-zero read count 
  if (normalised_PC_TCGA %>% 
      filter(transcript == gene) %>% 
      filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>%
      summarise(non_zero_counts = sum(`read_count normalised`> 0)) > 0) {
    
    #Extract stanfit object from TABI testing and coherse into data frame
    
    Fit<-as.data.frame(summary( #Summarise stanfit object
      # Call TABI
      (TABI_glm(formula = ~CAPRA_S+purity.score, 
                data = normalised_PC_TCGA %>% filter(transcript == gene) %>% 
                  select(CAPRA_S, purity.score, `read_count normalised`),
                prop_DE =  0.1,
                scale_DE =5,
                model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_vert_log_space_od.stan")))
      #Extract the stanfit object
      $fit)$summary) %>% #Extract summary of those 
      rownames_to_column("parameters") %>% #Make row names (which are names of paramters) explicit column 
      as_tibble() %>% 
      filter(!grepl("y_hat", parameters), !grepl("y_gen", parameters), !grepl("phi", parameters)) %>% #Remove rows with y_hat and y_gen 
      mutate(Gene_name = gene) #Add gene name 
    
    return(Fit)
    
  }
}


#Function for saving entire stanfit file
fit_total_purity<- function(gene_number) {
  #Name Gene to Test
  gene<- levels(normalised_PC_TCGA$transcript)[gene_number]
  
  
  #Required Packages
  
  require(TABI)
  require(dplyr)
  require(rstan)
  require(tibble)
  
  #Only evaluate if there is at least one non-zero read count 
  if (normalised_PC_TCGA %>% 
      filter(transcript == gene) %>% 
      filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>%
      summarise(non_zero_counts = sum(`read_count normalised`> 0)) > 0) {
    
    #Extract stanfit object from TABI testing and coherse into data frame
    
    #Summarise stanfit object
    # Call TABI
    TABI_glm(formula = ~purity.score +CAPRA_S, 
             data = normalised_PC_TCGA %>% filter(transcript == gene) %>% select(CAPRA_S, purity.score, `read_count normalised`),
             prop_DE =  0.1,
             scale_DE =5,
             model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_log_space_od.stan"))$fit
    
    
    
  }
}




#Function for saving entire TABI output



#Setting up parallelisation
library(doParallel)
library(bigstatsr)

cl <- parallel::makeCluster(bigstatsr::nb_cores())
doParallel::registerDoParallel(cl)

library(foreach)

#To create important fit table for plotting etc. 
#First 1000 genes

fit_table_purity_1000<-foreach(i=1:1000,.packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_purity(i) }

test_100<-foreach(i=1:100,.packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_purity(i) }





