

# Section 3 

# Function which performs differential expression using Bayseq

for(i in 2:n_group) {
  
  DE_result = result_list[[i-1]] %>% 
    as.data.frame() %>%
    select(names(result_list)[i-1] = DE)
  
}


Bayseq_DE<-function(simulated_data,
                          n_group,
                          x_cord = seq(from = -5, to =5, by = 0.5)) {
  # Set up
  
  {require(baySeq)
    require(dplyr)
    require(tidyr)
    require(tibble)
    library(doParallel)
    library(bigstatsr)
    library(foreach)
    library(doSNOW)}
  
  
  
  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
  
  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
  
  
  #Set up Read Count Data
  # Convert simulated data in form for baySeq analysis 
  # Data must be in the form – each column is a replicate / condition / CAPRA_S value
  # Each row is a transcript / gene 
  simdata <-simulated_data %>% 
    select(Gene_ref, Sample, value) %>% 
    tidyr::pivot_wider(names_from = Sample, values_from = value) %>%
    tibble::column_to_rownames(var = "Gene_ref")
  
  #Calculate sample size
  sample_size<-simulated_data %>% 
    select(Sample) %>% 
    dplyr::distinct() %>% 
    max()
  
  
  #Calculate total number of simulated genes
  n_sim_genes<-simulated_data %>% 
    select(Gene_ref) %>% 
    dplyr::distinct() %>% 
    nrow()
  
  
  # Number of distinct x-coordinate values / CAPRA - S values
  n_cord<- length(x_cord)
  
  #Define replicate names 
  #For each sample 
  replicates <-sapply(x_cord, 
                      function(x) paste0("CAPRA_S_", x))  %>% 
    c() %>% 
    sapply(., 
           function(x) rep(x, sample_size/n_cord))
  
  
  
  
  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
  
  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
  
  
  
  # Models / Groups
  {         
    #Null model (all groups are the same)
    NDE<-rep(1, 
             sample_size)
    
    # DE Model (using defined groups)
    DE<-match(make_groups(simulated_data,
                    n_group,
                    x_cord) %$% 
      group %$% 
      group, 
      LETTERS)-1
    
    models<-list(NDE, DE)
    } 
  
  
  # Set up parallelisation 
  {
    
    cl <- parallel::makeCluster(bigstatsr::nb_cores())
    doParallel::registerDoParallel(cl)
    
    registerDoSNOW(cl)} 

  # Breaking up dataset into genes sets of ~1000 
  if (sample_size<1001) {
    par_seq<-1000
  } if else (sample_size %% 1000 == 0) {
    par_seq<-c(1:(sample_size %/% 1000))
  } else {
    par_seq<-c(1:(sample_size %/% 1000))
    par_seq<-c(par_seq, 
               sum(par_seq) + sample_size %% 1000)
  }
  
  
  
  # Run analysis comparing Bayseq likelihood of DE groups to null (of no DE)
    
  Bayseq_DE_df=foreach(i=length(par_seq),
          .verbose = T) %do% {
            
            CD <- new("countData", 
                      data = simdata[parseq[i]:parseq[i+1],] %>% as.matrix(), 
                      replicates = replicates %>% c(),
                      groups = models)
           

            libsizes(CD) <- getLibsizes(CD)
            
            
 #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
            
            
            
            # Negative Binomial Model
            CD <- getPriors.NB(CD,
                               samplesize = 10000, #Recommended sample size of 10000
                               cl = cl)
            
            CD <- getLikelihoods(CD, 
                                 cl = cl, 
                                 bootStraps = 3, #Recommonded bootstrap number of 3
                                 verbose = TRUE)
          }
  
  
  return(Bayseq_DE_df)
}


# Take Bayseq results and compute FDR / ROC analysis 

FDR_ROC_Bayseq<-function(simulated_data, 
                         .sample = sample_id, #Sample ID column 
                         .transcript = Gene_ref, #transcipt / gene ID colum
                         .abundance = value, 
                         n_group,
                         x_cord = seq(from = -5, to =5, by = 0.5)) {
  
  # Set up columns
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  
  results_list = Bayseq_DE(simulated_data = simulated_data, 
                            .sample = .sample,
                            .transcript = .transcript,
                            .abundance = .abundance, 
                            n_group = n_group,
                            x_cord = x_cord
  )
  
  DE_Bayseq_result = list()
  
  for(i in 2:n_group) {
    
    DE_Bayseq_result = result_list[[i-1]] %>% 
      as.data.frame() %>%
      select(names(result_list)[i-1] = DE)
    
  }
  
  DE_result_table = do.call(cbind, 
                            DE_Bayseq_result) %>% 
    tibble::rowid_to_column("ID") 
  
  
  # Inflection results 
  
  pivot_longer()
  
  contrast_names=sapply(2:n_group, function(x) paste0(LETTERS[x-1], LETTERS[x]))
  
  contrast_values=make_groups(table = simulated_data,
                              .sample= !!.sample,
                              .transcript = !!.transcript,
                              n_group = nn_group,
                              x_cord = x_cord) %$% 
    contrast_values
  
  
  
  
  
  
  prob_range<-seq(from=0, to=1, by=1/1000)
  
  #Calculate number of Genes declared significant at each Probability value 
  
  N_sig<-sapply(prob_range, 
                function(x) {
                  full_table = results_table %>% 
                    exp() %>% 
                    as_tibble() %>% 
                    mutate(V1 = 1-V1) %>% 
                    filter(V1>=(1-x))
                  
                  
                  Total_n<-full_table %>% 
                    nrow()
                  
                  return(Total_n)
                  
                })
  
  # Calculate number of Genes incorrectly declared significant at each Probability value (i.e. False positives)
  
  N_sig_F<-sapply(prob_range, 
                  function(x) {
                    full_table = results_table %>% 
                      exp() %>% 
                      as_tibble() %>% 
                      mutate(V1 = 1-V1) %>% 
                      filter(V1>=(1-x))
                    
                    
                    False_n<-full_table %>% 
                      filter(grepl("V",
                                   names)) 
                    nrow()
                    
                    return(False_n)
                    
                  })
  
  True_FDR=N_sig_F/N_sig
  
  t_FDR_df<-data.frame(
    True_FDR = True_FDR
  )
  
  
  # ROC Statistics 
  
  # True Positives 
  
  P<-simulated_data %>% 
    select(!!.transcript) %>% 
    dplyr::distinct() %>% 
    filter(grepl("X", !!.transcript)) %>% 
    nrow()
  
  # True Negatives
  
  N<-simulated_data %>% 
    select(.!!transcript) %>% 
    dplyr::distinct() %>% 
    filter(grepl("V", !!.transcript)) %>% 
    nrow()
  
  #Calculate number of True Positives at each Probability Value
  
  TP<-sapply(pval_range, 
             function(x) {
               full_table<-cbind(results_table %>% 
                                   exp(), 
                                 names = simdata %>% rownames()) %>% 
                 as_tibble() %>% 
                 mutate(V1 = 1-V1) %>% 
                 filter(V1>=(1-x))
               
               True<-full_table %>% 
                 filter(grepl("X",names)) %>% 
                 nrow()
               
               return(True)
               
             })
  
  
  # For each probability threshold
  # how many false negatives
  # = number of positives - number of true positives
  FN<-P-TP
  
  
  # True positive rate
  # = Number of true positives / number of positives
  TPR<-TP/P
  
  # False negative rate
  # = Number of false negatives / number of positives
  FNR<-FN/P
  
  ROC_df<-data.frame(
    False_Neg = FN,
    True_Pos = TP,
    Pos = P,
    TPR = TPR,
    FNR = FNR
  )  
  
  
  # Estimation of inflection 
  
  
  
  return(list(
    #pred_FDR = FDR_table,
    true_FDR = t_FDR_df,
    ROC = ROC_df,
    Inflect_est = Inflect_df
  ))
} 
  