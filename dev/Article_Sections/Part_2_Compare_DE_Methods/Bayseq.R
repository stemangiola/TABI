#Bayseq Analysis 


#BiocManager::install("baySeq")

{library(baySeq)
library(dplyr)
library(tidyr)
library(tibble)}


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
Sig_multiple_value_data<-test %>% mutate(Sample = rep(1:105, 15000))

#Set up Read Count Data
# Convert simulated data in form for baySeq analysis 
# Data must be in the form – each column is a replicate / condition / CAPRA_S value
# Each row is a transcript / gene 
simdata <-Sig_multiple_value_data %>% 
  select(Gene_ref, Sample, value) %>% 
  tidyr::pivot_wider(names_from = Sample, values_from = value) %>%
  tibble::column_to_rownames(var = "Gene_ref")


sample_size<-Sig_multiple_value_data %>% 
  select(Sample) %>% 
  distinct() %>% 
  max()


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

         

#Define replicate names 
replicates <-sapply(seq(from = -5, to=5, by = 0.5), function(x) paste0("CAPRA_S_", x))  %>% 
           c() %>% 
  sapply(., function(x) rep(x, sample_size/21))


         
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




# Models / Groups
{         
#Null model (all groups are the same)
NDE<-rep(1, 
         21)
                  
#2 Group Model, there is a difference between early and late expression
# Early defined as CAPRA_S<0, Late >=0
DE_2_group<- sapply(seq(from = -5, to=5, by = 0.5),  
                    function(x) ifelse(x<0, 1, 2))

                  
#3 Group Model: Early / Middle / Late
# Early <= CAPRA_S  -2
# Middle  -2 < CAPRA_S < 2
# Late >= CAPRA_S 2
                  
                  
#3 Groups Model, there is a difference between early, middle and late expression
DE_3_group_all_diff<- sapply(seq(from = -5, to=5, by = 0.5),  
                             function(x) ifelse(x<=2*-1, 1, ifelse(x>=2, 3, 2)))
                  
#3 Groups Model, There is a difference between early / middle but not late / middle
                  DE_3_group_early_change<- sapply(seq(from = -5, to=5, by = 0.5),  
                                                   function(x) ifelse(x<=2*-1, 1, 2))
                  
                  # 3 Groups Model, There is a difference between late / middle but not middle / early
                  DE_3_group_late_change<-sapply(seq(from = -5, to=5, by = 0.5),  
                                                 function(x) ifelse(x<2, 1, 2))
                  
                  
                  
#4 Group Model 
# Group A (Early) 
# Group B (Early Middle)
# Group C (Late Middle)
# Group D (Late)
                  
# Difference between all groups (A,B,C,D)
DE_4_group_all_diff<- sapply(seq(from = -5, to=5, by = 0.5),  
                                               function(x) ifelse(x<2.5*-1, 1, ifelse(x<0, 2, ifelse(x<2.5, 3, 4))))
                  

# Difference between Group A and all others 
DE_4_group_A_diff<- sapply(seq(from = -5, to=5, by = 0.5),  
                                             function(x) ifelse(x<2.5*-1, 1, 2)) 
                            
                 
# Difference between Group D and all others 
DE_4_group_D_diff<- sapply(seq(from = -5, to=5, by = 0.5),  
                           function(x) ifelse(x<2.5, 1,2))

# A / B same, C different, D different 
DE_4_group_AB_same <- sapply(seq(from = -5, to=5, by = 0.5),  
                             function(x) ifelse(x<2.5*-1, 1, ifelse(x<0, 1, ifelse(x<2.5, 2, 3))))
                                             
}


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



# Convert to having more than one sample per group 

convert_sample<- function(total_sample_size) {
  require(dplyr)
  
  # Groups 
  groups<-list(NDE, 
               DE_2_group,  
               DE_3_group_all_diff, 
               DE_3_group_early_change, 
               DE_3_group_late_change, 
               DE_4_group_all_diff, 
               DE_4_group_A_diff, 
               DE_4_group_D_diff, 
               DE_4_group_AB_same) 
  
 convert<-lapply(groups, 
                function(x) {
                sapply(x, 
                       function(y) rep(y, total_sample_size/21) 
                       ) %>%  c()}
                )
 
 return(convert)
}

groups<-convert_sample(sample_size)

groups[[1]]<-rep(1,sample_size)
 
 
 
 
 
                                             
 #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
                                            
 
 
 
# Count Data 

 
CD <- new("countData", 
          data = simdata, 
          replicates = c(replicates),
          groups = groups)
                                             

libsizes(CD) <- getLibsizes(CD)
                                             
# Set up Parallelisation

{library(doParallel)
library(bigstatsr)
  
                                               
cl <- parallel::makeCluster(bigstatsr::nb_cores())
doParallel::registerDoParallel(cl)
                                      
         
library(foreach)
library(doSNOW)
                                               
 registerDoSNOW(cl)} 
                            

                 
                                             
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


                                             
# Negative Binomial Model
CD <- getPriors.NB(CD,
                   samplesize = 10000, #Recommended sample size of 10000
                   cl = cl)
                                             
CD <- getLikelihoods(CD, 
                     cl = cl, 
                     bootStraps = 3, 
                     verbose = TRUE)


CD@estProps

test<-CP

PEP_table<-exp(test)[,2] %>% 
  as_tibble() %>% 
  mutate(PIP = value) %>% 
  mutate(PEP = 1-PIP) %>% 
  arrange(PEP) %>% 
  mutate(qvalue = cummean(PEP))# %>% 
  ggplot(aes(x=PEP)) + 
  geom_histogram()

N_sig_expFDR<-sapply(seq(from =0, to = 1, by = 1/1000), # Each FDR 
                     function(x) PEP_table %>% 
                       filter(qvalue<x) %>% 
                       nrow() 
)

topCounts(CD, group = 2) %>% dim()


cbind(sim_data %>% 
  dplyr::select(Gene_number) %>% 
  dplyr::distinct(), PEP_table)

data.frame(N_sig = N_sig_expFDR,
           FDR_pred = seq(from =0, to = 1, by = 1/1000)) %>% 
  ggplot(aes(x=FDR_pred, N_sig)) + 
  geom_line()

Bayseq_3_group=foreach(i=c(seq(from = 1, to = 15000, by =1000)),
                                .verbose = T) %do%  {
                                  
                                  CD <- new("countData", 
                                            data = simdata[i:(999+i),] %>% as.matrix(), 
                                            replicates = replicates %>% c(),
                                            groups = list(groups[[1]], groups[[3]], groups[[4]], groups[[5]]))
                                  
                                  
                                  
                                  
                                  libsizes(CD) <- getLibsizes(CD)
                                  
                                  # Set up Parallelisation
                                  
                                  
                                  
                                  #Bayseq_2_group_pval
                                  
                                  
                                  
                                  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
                                  
                                  
                                  
                                  # Negative Binomial Model
                                  CD <- getPriors.NB(CD,
                                                     samplesize = 10000, #Recommended sample size of 10000
                                                     cl = cl)
                                  
                                  CD <- getLikelihoods(CD, 
                                                       cl = cl, 
                                                       bootStraps = 3, 
                                                       verbose = TRUE)
                                  
                                  #saveRDS(CD@posteriors, compress = "gzip", file = paste0("Bayseq", i, "2_group_2000_2.rds"))
                                  
                                  CD@posteriors
                                  
                                }
saveRDS(Bayseq_2_group, compress = "gzip", file = "Bayseq_2_group.rds")



Bayseq_2_group=foreach(i=c(seq(from = 1, to = 15000, by =1000)),
                       .verbose = T) %do%  {
                         
                         CD <- new("countData", 
                                   data = simdata[i:(999+i),] %>% as.matrix(), 
                                   replicates = replicates %>% c(),
                                   groups = list(groups[[1]], groups[[2]]))
                         
                         
                         
                         
                         libsizes(CD) <- getLibsizes(CD)
                         
                         # Set up Parallelisation
                         
                         
                         
                         #Bayseq_2_group_pval
                         
                         
                         
                         #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
                         
                         
                         
                         # Negative Binomial Model
                         CD <- getPriors.NB(CD,
                                            samplesize = 10000, #Recommended sample size of 10000
                                            cl = cl)
                         
                         CD <- getLikelihoods(CD, 
                                              cl = cl, 
                                              bootStraps = 3, 
                                              verbose = TRUE)
                         
                         #saveRDS(CD@posteriors, compress = "gzip", file = paste0("Bayseq", i, "2_group_2000_2.rds"))
                         
                         CD@posteriors
                         
                       }








Bayseq_test<-lapply(Bayseq_test, function(x) {x %>% 
                       as_tibble()}) %>% 
  bind_rows()

Prob_table<-cbind(Bayseq_test %>% 
  exp(),
  Gene_ref = simdata %>% rownames()) 

saveRDS(Prob_table, 
        compress = "gzip", 
        file = paste0("Bayseq_2_group_samplesize", k, "_boostrap", j, ".rds"))
} 
  
  
  
  
  ROC_Bayseq<-function(table, simulated_data) {
    
    simdata <-simulated_data %>% 
      select(Gene_ref, Sample, value) %>% 
      tidyr::pivot_wider(names_from = Sample, values_from = value) %>%
      tibble::column_to_rownames(var = "Gene_ref")
    
    N_sig<-sapply(seq(from=0, to=1, by=1/1000), 
                  function(x) {
                    full_table<-cbind(table %>% 
                                        exp(), names = simdata %>% rownames()) %>% 
                      as_tibble() %>% 
                      mutate(V1 = 1-V1) %>% 
                      filter(V1>=(1-x))
                    
                    #False<-full_table %>% filter(grepl("V",names)) %>% nrow()
                    
                    Total<-full_table %>% nrow()
                    
                    return(Total)
                    
                  })
    
    FP<-sapply(seq(from=0, to=1, by=1/1000), 
               function(x) {
                 full_table<-cbind(table %>% 
                                     exp(), names = simdata %>% rownames()) %>% 
                   as_tibble() %>% 
                   mutate(V1 = 1-V1) %>% 
                   filter(V1>=(1-x))
                 
                 False<-full_table %>% filter(grepl("V",names)) %>% nrow()
                 
                 #Total<-full_table %>% nrow()
                 
                 return(False)
                 
               })
    
    
    
    
    return(data.frame(FDR_true = FP/N_sig,
                      N_sig = N_sig))
    
  } 
  
  
  FDR_ROC_Bayseq<-function(table, simulated_data) {
    
    simdata <-simulated_data %>% 
      select(Gene_ref, Sample, value) %>% 
      tidyr::pivot_wider(names_from = Sample, values_from = value) %>%
      tibble::column_to_rownames(var = "Gene_ref")
    
    N_sig<-sapply(seq(from=0, to=1, by=1/1000), 
                  function(x) {
                    full_table<-cbind(table %>% 
                                        exp(), names = simdata %>% rownames()) %>% 
                      as_tibble() %>% 
                      mutate(V1 = 1-V1) %>% 
                      filter(V1>=(1-x))
                    
                    #False<-full_table %>% filter(grepl("V",names)) %>% nrow()
                    
                    Total<-full_table %>% nrow()
                    
                    return(Total)
                    
                  })
    
    FP<-sapply(seq(from=0, to=1, by=1/1000), 
               function(x) {
                 full_table<-cbind(table %>% 
                                     exp(), names = simdata %>% rownames()) %>% 
                   as_tibble() %>% 
                   mutate(V1 = 1-V1) %>% 
                   filter(V1>=(1-x))
                 
                 False<-full_table %>% filter(grepl("V",names)) %>% nrow()
                 
                 #Total<-full_table %>% nrow()
                 
                 return(False)
                 
               })
    
    
    return(data.frame(FDR_true = FP/N_sig,
                      N_sig = N_sig))
    
  } 
} 
                                             
