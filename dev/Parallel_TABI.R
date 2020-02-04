#Set Up Clusters For Parallelisation

library(doParallel)
library(bigstatsr)

cl <- parallel::makeCluster(bigstatsr::nb_cores())
doParallel::registerDoParallel(cl)

#Change all CAPRA_S 8 to 7 

library(dplyr)

library(foreach)
normalised_PC_TCGA <- normalised_PC_TCGA %>% filter(is.na(CAPRA_S)==F, is.na(transcript)==F)

normalised_PC_TCGA <- normalised_TCGA


saveRDS(normalised_PC_TCGA, compress="gzip", file="compress_normalised_narm_TCGA.rds")

for (j in 1:length(normalised_PC_TCGA$CAPRA_S)) {
  if (normalised_PC_TCGA$CAPRA_S[j]>7) {
    normalised_PC_TCGA$CAPRA_S[j]=7
    print(j)
  }
  else {
  print(j) }
}

saveRDS(normalised_PC_TCGA, compress="gip", file="compress_normalised_narm_TCGA")

#Check is.na

summary(normalised_PC_TCGA$CAPRA_S)

saveRDS(normalised_PC_TCGA, compress = "gzip", file="compress_normalised_na_rm_TCGA")


#Create Fitting Function 
#Compresssed Fitting function (just )

fit<- function(gene_number) {
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
    (TABI_glm(formula = ~CAPRA_S, 
              data = normalised_PC_TCGA %>% filter(transcript == gene) %>% filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>% select(CAPRA_S, `read_count normalised`),
              prop_DE =  0.1,
              scale_DE =5,
              model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_log_space_od.stan")))
    #Extract the stanfit object
    $fit)$summary) %>% #Extract summary of those 
    rownames_to_column("parameters") %>% #Make row names (which are names of paramters) explicit column 
    as_tibble() %>% 
    #filter(!grepl("y_hat", parameters), !grepl("y_gen", parameters), !grepl("phi", parameters)) %>% #Remove rows with y_hat and y_gen 
    mutate(Gene_name = gene) #Add gene name 
  
  return(Fit)
  
  }
}

normalised_PC_TCGA %>% names()

tab<-fit(8)
tab %>% filter(grepl("log_y_hat", parameters))
#phi 0.03 to 4
#log_y_hat -6 to -2.70 
exp(-6)
exp(-2.70)
rnbinom(n=20, mu= exp(-2.7), size = exp(-6))



(TABI_glm(formula = ~CAPRA_S, 
          data = normalised_PC_TCGA %>% filter(transcript == "A2ML1-AS2") %>% filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>% select(CAPRA_S, `read_count normalised`),
          prop_DE =  0.1,
          scale_DE =5,
          model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_log_space_od.stan")))

#Running on problem simulated data

library(dplyr)
library(TABI)
library(ggplot2)

data_disp <- data.frame(CAPRA_S = seq(from=0, to=7, by=0.5),  #Simulate random CAPRA_S score
                        read =  rnbinom(n = length(seq(from=0, to=7, by=0.5)), size =  0.005, mu = (seq(from=0, to=7, by=0.5))))  #with negbinom distribution
data_disp %>% 
  ggplot(aes(y=read, x=CAPRA_S)) +geom_point() #Plot simulated data 

#Run TABI

TABI_glm(formula = ~CAPRA_S, 
          data = data_disp,
          prop_DE =  0.1,
          scale_DE =5,
          model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_vert_log_space_od.stan"))







#Complete $fit Save function 
fit_total<- function(gene_number) {
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
      TABI_glm(formula = ~CAPRA_S, 
                data = normalised_PC_TCGA %>% filter(transcript == gene) %>% filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>% select(-transcript, -TMM, -purity.score),
                prop_DE =  0.1,
                scale_DE =5,
                model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_vert_log_space_od.stan"))$fit
    
    
    
  }
}

#Fit for CAPRA_S + purity.score (purity as a covariate)
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
                data = normalised_PC_TCGA %>% filter(transcript == gene) %>% filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>% select(-transcript, -TMM, -purity.score),
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

test_disp$generated_quantities
#Complete Fit Save
fit_total<- function(gene_number) {
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
    TABI_glm(formula = ~CAPRA_S, 
             data = normalised_PC_TCGA %>% filter(transcript == gene) %>% filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>% select(-transcript, -TMM, -purity.score),
             prop_DE =  0.1,
             scale_DE =5,
             model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_log_space_od.stan"))$fit
    
    
    
  }
}


#Create parallelisation loop

library(foreach)

#Number of 
length(levels(normalised_PC_TCGA$transcript))

fit_purity_p1 <- foreach(i=sample.int(37318, 100), .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit(i)
}

test_parameters <- fit_purity_p1 

saveRDS(test_parameters, compress='gzip', file="test_parameters")

saveRDS(fit_purity_p1, compress="gzip", file="fit_purity_p1.rds")


library(TABI)

TABI_glm(formula = ~CAPRA_S+purity_score, 
         data = normalised_PC_TCGA %>% filter(transcript == "AMBN") %>% filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>% select(-transcript, -TMM),
         prop_DE =  0.1,
         scale_DE =5,
         model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_log_space_od.stan"))


normalised_PC_TCGA <- normalised_PC_TCGA %>% rename(purity_score = purity.score)


fit_ed_p1<-foreach(i=1:5000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit(i)
}

#Rerun TABI fitting function on genes for

#altered normalised_PC_TCGA with CAPRA_S 8 -> 7 

ed_normalised_PC_TCGA <-normalised_PC_TCGA

saveRDS(ed_normalised_PC_TCGA, file="ed_normalised_PC_TCGA.rds", compress="gzip")

save(ed_normalised_PC_TCGA, file="ed_normalised_PC_TCGA.rda")

fit_ed_p1<-foreach(i=1:5000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit(i)
}

saveRDS(fit_ed_p1, compress="gzip", file="fit_ed_p1.rds")

save(fit_ed_p1, file="fit_ed_p1.rda")

fit_ed_p2<-foreach(i=5001:10000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit(i)
}

saveRDS(fit_ed_p2, compress="gzip", file="fit_ed_p2.rds")

fit_ed_total_p1<-foreach(i=1:5000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble')) %dopar% {
  fit_total(i)
}

saveRDS(fit_ed_total_p1, compress="gzip", file="fit_ed_total_p1.rds")


fit_ed_p3<-foreach(i=10001:15000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit(i)
}

saveRDS(fit_ed_p3, compress="gzip", file="fit_ed_p3.rds")

fit_ed_p4<-foreach(i=15001:20000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit(i)
}

saveRDS(fit_ed_p4, compress="gzip", file="fit_ed_p4.rds")



Greater_total_1000<-foreach(i=1:1000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_total(i)
}

# Simulate Data with mu = exp(CAPRA_S) and low NB
# Just dispersion - no differential change 


data_disp <- data.frame(CAPRA_S = seq(from=0, to=7, by=0.005),  
                        read =  rnbinom(n = length(seq(from=0, to=7, by=0.005)), size = 0.01, mu = exp(seq(from=0, to=7, by=0.005)))) 

data_disp %>% 
  ggplot(aes(x = CAPRA_S, y = read)) +geom_point()

test_disp<-TABI_glm(
  formula = ~CAPRA_S, 
  data  = data_disp, 
  prop_DE =  0.1,
  scale_DE =5,
  model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_vert_log_space_od.stan", auto_write = F)
  )

#Comparing TABI generated data with + simulated data

test_disp$generated_quantities %>% 
  ggplot(aes(x= seq(from=0, to=7, by=0.05), y=estimate)) + 
  geom_point(col="blue") +
  geom_point(aes(x=CAPRA_S, y=read), data=data-disp)
  scale_y_log10()
  
tp<-TABI_glm(
  formula =~CAPRA_S, 
  data = normalised_PC_TCGA %>% filter(transcript =="A1BG") %>% filter(is.na(CAPRA_S)==F, is.na(`read_count normalised`==F)) %>% select(-TMM, - purity.score, -transcript),
  prop_DE =  0.1,
  scale_DE =5,
  model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_vert_log_space_od.stan"))

TABI_glm(
  formula=~purity.score +CAPRA_S, 
  data = data, 
  prop_DE = 0.1,
  scale_DE =5,
  model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_log_space_od.stan")
)

#Real Data that is causing problems
#Problems 8, 23627

normalised_PC_TCGA %>% 
  filter(transcript== "A2ML1-AS2") %>% 
  ggplot(aes(x=CAPRA_S, y=`read_count normalised`+1)) + 
  geom_point() +
  scale_y_log10()


normalised_PC_TCGA %>% 
  filter(transcript==  "AMBN") %>% 
  ggplot(aes(x=CAPRA_S, y=`read_count normalised`)) + 
  geom_point() +
  scale_y_log10()

reformate_PC_CAPRA_S <- normalised_PC_TCGA %>%
  filter(transcript = c("A1BG", "A1BG-AS1", "A1CF", "A2M")) %>% 
  group_by(transcript) %>%
  mutate(row = row_number()) %>% 
  filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>% 
  select(-TMM, -purity.score) %>% 
  pivot_wider(names_from=transcript, values_from =`read_count normalised`)




test_A2<-TABI_glm(
  formula = ~CAPRA_S,
  data = normalised_PC_TCGA %>% 
    filter(transcript== "A1BG") %>% select(-TMM, -transcript, -purity.score) %>% filter(is.na(CAPRA_S) == F, is.na(`read_count normalised`)==F),
  prop_DE =  0.1,
  scale_DE =5,
  model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_log_space_od.stan")
  
)


normalised_PC_TCGA %>% 
  filter(transcript== "A1BG") %>% 
  select(-TMM, -transcript, -purity.score) %>% 
  filter(is.na(CAPRA_S) == F, is.na(`read_count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y= `read_count normalised`)) +geom_point(col="blue") + geom_point(aes(x=CAPRA_S, y= estimate), data =
                                                                                  test_A2$generated_quantities %>% ungroup %>% mutate(CAPRA_S = 
                                                                                                                            (normalised_PC_TCGA %>% 
                                                                                                                             filter(transcript== "A1BG") %>% 
                                                                                                                               filter(is.na(CAPRA_S) == F, is.na(`read_count normalised`)==F))$CAPRA_S))
# Differential Change + Dispersion


#Create Fitting Function 
#Compresssed Fitting function (just )

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
      (TABI_glm(formula =~purity.score +CAPRA_S, 
                data = normalised_PC_TCGA %>% filter(transcript == gene) %>% filter(is.na(transcript)==F, is.na(CAPRA_S)==F) %>% select(CAPRA_S, purity.score, `read_count normalised`),
                prop_DE =  0.1,
                scale_DE =5,
                model = rstan::stan_model("~/TABI/inst/stan/DE_sigmoid_one_gene_log_space_od.stan")))
      #Extract the stanfit object
      $fit)$summary) %>% #Extract summary of those 
      rownames_to_column("parameters") %>% #Make row names (which are names of paramters) explicit column 
      as_tibble() %>% 
      #filter(!grepl("y_hat", parameters), !grepl("y_gen", parameters), !grepl("phi", parameters)) %>% #Remove rows with y_hat and y_gen 
      mutate(Gene_name = gene) #Add gene name 
    
    return(Fit)
    
  }
}

purity_score_fit<-foreach(i=1:100, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_purity(i)
}

purity_score_fit_p1<-foreach(i=1:1000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_purity(i)
}

saveRDS(purity_score_fit_p1, compress="gzip", file="purity_score_fit_1.rds")

fit_non_purity_p1<-foreach(i=1:1000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit(i) }

saveRDS(fit_non_purity_p1, compress="gzip", file="fit_non_purity_p1.rds" )

purity_score_fit_p2<-foreach(i=1001:2000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_purity(i)
}

saveRDS(purity_score_fit_p2, compress="gzip", file="purity_score_fit_2.rds")

purity_score_fit_p3<-foreach(i=2001:3000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_purity(i)
}


saveRDS(purity_score_fit_p3, compress="gzip", file="purity_score_fit_3.rds")

purity_score_fit_p4<-foreach(i=3001:4000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_purity(i)
}

saveRDS(purity_score_fit_p4, compress="gzip", file="purity_score_fit_4.rds")

purity_score_fit_p5<-foreach(i=4001:5000, .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_purity(i)
}

saveRDS(purity_score_fit_p5, compress="gzip", file="purity_score_fit_5.rds")

purity_score_random_1000<-foreach(i=sample.int(size=1000, n = length(levels(normalised_PC_TCGA$transcript))), .packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit_purity(i)
}

saveRDS(purity_score_random_1000, compress="gzip", file="purity_score_random_5000.rds")

fit_ed_random_1000<-foreach(i=sample.int(size=1000, n = length(levels(normalised_PC_TCGA$transcript))),.packages=c('TABI', 'dplyr', 'rstan', 'tibble'), .combine="rbind") %dopar% {
  fit(i) }

saveRDS(fit_ed_random_1000, compress="gzip", file="fit_ed_random_5000.rds")




