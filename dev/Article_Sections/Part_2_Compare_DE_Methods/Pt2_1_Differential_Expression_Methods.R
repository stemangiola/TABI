


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #  #




#Intention of Script 
#Create Functions which run tidybulk differential expression analysis (e.g. edgeR, DESeq2, limma + voom)
#on simulated data and then calculate FDR, ROC and estimate inflection values (so can be compared to TABI)


# (Section 1) make_groups function
# Function which is able to make a data frame to match CAPRA_S values to group labels
# Given a desired number of groups to sort the simulated samples into 


# (Section 2) 
# Function which performs differential expression using tidybulk
# Can be edgeR, DESEq2, limma + voom


# (Section 3)
# Function which performs differential expression using Bayseq


# (Section 4)
# Function which takes table of differential expression from tidybulk
# Computes FDR, ROC and inflection estimate tables


# (Section 5)
# Function which takes table of differential expression from Bayseq
# Computes FDR, ROC and inflection estimate tables





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


# (Section 1)

# Make groups for using tidybulk differential expression


# Input simulated df table from function sigmoidal_sim_df 

# Output: list consisting of $group (df matching CAPRA-S / x coordinate value to group letter)
# and $contrast_values (vector of values which correspond to the %>% ddle value between each group  - 
# i.e. if highest CAPRA-S for group A is 4 and lowest for group B is a 5 then the contrast value is 4.5)


make_groups<-function(table, #Table of simulated RNA seq data counts
                      .sample = sample_id, #Sample ID column 
                      .transcript = Gene_ref, #transcript / gene ID column
                      .abundance = value, # read count columns
                      .x = CAPRA_S,# x - coordinate column / factor of interest, in this case CAPRA_S
                      n_group, # number of groups 
                      x_cord = seq(from = -5, to =5, by = 0.5) #For grouping, vector of discrete x-coordinates values / simulated CAPRA values
                        ){

  
  # Set up columns
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  .x = enquo(.x)
  
  # Number of distinct simulated x-coordinate values 
  n_cord = length(x_cord)

  # Reduce df of simulated data to only required samples
  # sample id column (sample)
  # name of transcript 
  # RNAseq count column
  table = table %>% 
    select(!!.sample,
           !!.transcript,
           !!.abundance)
   

  
  if((n_cord %% n_group) ==0){ #If number of groups to make is a multiple of possible x_cord values
    
    # Then number of different CAPRA S value or x-coordinate values with per groups is : 
    n_per_group<-n_cord/n_group # Number of different x-cord per group 
    
    # Vector of number of samples in each group (assuming one sample per CAPRA-S value)
    groups<-rep(n_per_group, 
                n_group)
    
    x_cord_id<-cumsum(groups)[-length(groups)]

    # Contrast values are the point estimate of the 'jump' 
    # i.e. if there is a change detected between the two groups, what is the estimate
    contrast_values<-(x_cord[x_cord_id]+x_cord[x_cord_id+1])/2
    # # Assign each sample a group letter as a group id
  
    
    # Assign each x-coordinate value a group letter based on how many x-coordinate values there are per group
    group_let<-sapply(LETTERS[1:n_group], 
                      function(x) rep(x, n_per_group)) %>% 
                       c()
    

  
    # Combine group id and corresponding x-cord value into a data frame
    df<-data.frame(group = group_let,
                   CAPRA_S = x_cord)

    
  } 
  else{ #if the number of groups is not a multiple of possible x_cord values
    # Then group 'lower' x-coordinate values together, such that there are bigger groups
    # with earlier CAPRA-S values
    
    n_per_group<-n_cord %/% n_group # How many samples per group, if only accept whole numbers
    
    n_left_over<-n_cord %% n_group # How many samples would be left over if we used the above number

    # Vector of number of samples in each group (assuming one sample per CAPRA-S value)
    # Give early / left most groups add an additional sample
    # than late / right most groups
    groups<-c(rep(n_per_group+1,n_left_over), 
              rep(n_per_group, n_group-n_left_over))
    
    

    group_let<-sapply(1:n_group, 
                      function(x) rep(LETTERS[x], groups[x])) %>% 
               Reduce(c, .)
    
   
    x_cord_id<-cumsum(groups)[-length(groups)]
  
    
    # Contrast values are the point estimate of the 'jump' 
    # i.e. if there is a change detected between the two groups, what is the estimate
    contrast_values<-(x_cord[x_cord_id]+x_cord[x_cord_id+1])/2 
    
    
    # Combine group id and corresponding x-cord value into a data frame
    df<-data.frame(group = group_let,
                   CAPRA_S = x_cord)
    
  } 
  
  
# Return data frame of group id and corresponding x-cord value
  # And point estimate the point estimate of the 'jump' between two groups
  
  return(list(group = df,
              contrast_values = contrast_values)) 
  
  
}



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




# Section 2 

# Function which performs differential expression using tidybulk
# Can be edgeR, DESEq2, limma + voom


tidybulk_DE<-function(table, #Table of simulated RNA seq data counts
                      .sample = sample_id, #Sample ID column 
                      .transcript = Gene_ref, #transcript / gene ID colum
                      .abundance = value, 
                      .x = CAPRA_S,
                      x_cord = seq(from = -5, to = 5, by = 0.5),
                     n_group, # Number of groups for DE 
                     method = "edgeR_quasi_likelihood"){
  
  
  require(tidybulk)
  require(dplyr)
  require(tidyr)
  require(tibble)
  
  
  # Set up columns
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  .x = enquo(.x)
  

  # Combine group names created using make_groups function
  # with table of simulated RNA seq counts
  tidybulk_df<-dplyr::left_join(table, 
                  make_groups(table %>% dplyr::rename(CAPRA_S = !!.x),
                              .sample = !!.sample,
                              .transcript = !!.transcript,
                              .abundance = !!.abundance,
                              .x = !!.x, 
                              x_cord = x_cord,
                              n_group = n_group) %$% group) %>% 
    as_tibble()   %>%
    mutate(alpha = as.numeric(alpha))  %>%
    mutate(Sample = as.factor(!!.sample))  %>%
    distinct()  %>%
    tidybulk::tidybulk(.sample = Sample,
                       .transcript = !!.transcript,
                       .abundance = !!.abundance)

  # Create contrasts to feed to tidybulk::test_differential_abundance function
  
  if(method == "DESeq2"){
    
    contrasts = list() #Then the contrasts should actually be a list of contrasts, in the form
    
    for (i in 2:(n_group)) {
      
      contrasts[[i-1]]= c("group", LETTERS[i], LETTERS[i-1])
    }
    
  }
  
  else{# Else the contrasts should be a vector to send to test_differential_abundance
  contrasts<-sapply(2:n_group, 
                    function(x) 
                      paste0("group", LETTERS[x], 
                             "-", 
                             "group", LETTERS[x-1]))
  
  }
  
  # Run DE analysis using specified parameters
  
  tidybulk_analysis<-tidybulk::test_differential_abundance(
    tidybulk_df,
    ~0 + group,
    .contrasts = contrasts,
    method = method,
    scaling_method = "none",
    omit_contrast_in_colnames = FALSE,
    #prefix = "",
    action = "get",
    significance_threshold = 1, # Return all results 
    fill_missing_values = FALSE
  )
  
  return(tidybulk_analysis)
} 


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #





# Section 3 

# Function which performs differential expression using Bayseq


Bayseq_DE<-function(simulated_data,
                    .sample = sample_id, #Sample ID column 
                    .transcript = Gene_ref, #transcipt / gene ID colum
                    .abundance = value, 
                    n_group,
                    x_cord = seq(from = -5, to =5, by = 0.5)) {
  
  
  # Set up columns
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  
  
  
  # Set up
  
  {require(baySeq)
    require(dplyr)
    require(tidyr)
    require(tibble)
    require(magrittr)
    #require(rlang)
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
    select(!!.transcript, !!.sample, !!.abundance) %>% 
    tidyr::pivot_wider(names_from = !!.sample, values_from = !!.abundance) %>%
    tibble::column_to_rownames(var = rlang::as_name(.transcript))


    #Calculate sample size
  sample_size<-simulated_data %>% 
    select(!!.sample) %>% 
    dplyr::distinct() %>% 
    max()
  
  #Calculate total number of simulated genes
  n_sim_genes<-simulated_data %>% 
    select(!!.transcript) %>% 
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
    
    # DE Model (using defined groups)
    all_DE<-match(make_groups(simulated_data,
                          .sample = !!.sample,
                          .transcript = !!.transcript,
                          .abundance = !!.abundance,
                          n_group = n_group,
                          x_cord = x_cord) %$% 
                group %$% 
                group, 
              LETTERS)
    } 
  
    DE_list = list()
    
    sim_df_list = list()
    

      groups_names = make_groups(
      table = simulated_data,
      .sample = !!.sample,
      .transcript = !!.transcript,
      .abundance = !!.abundance,
      n_group = n_group,
      x_cord = x_cord
    ) %$% 
      group %$% 
      group
    
    for(i in 2:n_group) {
    
    DE = replace(all_DE,
            which(all_DE == i-1),
            LETTERS[i-1])

    DE = replace(DE,
                 which(all_DE == i),
                 LETTERS[i])
    DE = replace(DE,
                 which(DE != LETTERS[i-1] & DE != LETTERS[i]),
                 NA) %>%
      na.omit()

    sample_id_to_keep = which(groups_names == LETTERS[i-1]|groups_names == LETTERS[i]) 
    
    #assign(paste0(LETTERS[i-1], LETTERS[i]), DE) 
    
    DE_list[[i-1]]<-match(DE, LETTERS)
    
    sim_df_list[[i-1]] = sample_id_to_keep
        
    }
    
    names(DE_list)<- sapply(2:n_group, 
           function(x) paste0(LETTERS[x-1], LETTERS[x]))
    
    #models<-list(NDE, DE_list)
  
  # Set up parallelisation 
  {
    
    cl <- parallel::makeCluster(bigstatsr::nb_cores())
    doParallel::registerDoParallel(cl)
    
    registerDoSNOW(cl)} 
  
  # Breaking up dataset into genes sets of ~1000 
  if(n_sim_genes<1001) {
    par_seq<-c(1, n_sim_genes)
  } else if (n_sim_genes %% 1000 == 0) {
    par_seq<-c(1, 1:(n_sim_genes %/% 1000)*1000)
  } else {
    par_seq<-c(1, 1:(n_sim_genes %/% 1000)*1000)
    par_seq<-c(par_seq, 
               sum(par_seq) + n_sim_genes %% 1000)
  }
  
  
  # Run analysis comparing Bayseq likelihood of DE groups to null (of no DE)
  
  #return(list(models, replicates, simdata)) } 
   
  #return(par_seq)} 

    # DE = rep(1, length(DE_list[[1]]))
    # 
    # return(list(simdata[, sim_df_list[[1]]] ,
    #             (replicates %>%
    #                c())[sim_df_list[[1]]],
    #             list(DE = DE, DE_list[[1]]) )) }

    Bayseq_DE_df=foreach(j = 1:length(DE_list)) %do% {
    
      
    DE = rep(1, length(DE_list[[j]]))
    
    models= list(DE = DE, 
                 DE_list[[j]] %>% as.vector()) 
    

    filt_simdata = simdata[, sim_df_list[[j]]] 
    
    
    replicates_sub = (replicates %>% 
                        c())[sim_df_list[[j]]] 
    
  #   return(list(filt_simdata, replicates_sub, models))
  #   } 
  # return(Bayseq_DE_df)  
# } 
           foreach(i=1:(length(par_seq)-1),
                       .verbose = T, 
                      .combine = 'rbind') %do% {
                         
                        if(i==1) {
                        
                         CD <- new("countData", 
                                   data = filt_simdata[par_seq[i]:par_seq[i+1],], 
                                   replicates = replicates_sub,
                                   groups = models)
                        } 
                        
                        else {
                          CD <- new("countData", 
                                    data = filt_simdata[(par_seq[i-1]+1):par_seq[i+1],], 
                                    replicates = replicates_sub,
                                    groups = models)
                          
                        }
                         
                         libsizes(CD) <- getLibsizes(CD)
                         
  
                         #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
                         
                         
                         
                         # Negative Binomial Model
                         CD_priors <- getPriors.NB(CD,
                                            samplesize = 10000, #Recommended sample size of 10000
                                            cl = cl)
                         
                         CD_likelihoods <- getLikelihoods(CD_priors, 
                                              cl = cl, 
                                              bootStraps = 3, #Recommonded bootstrap number of 3
                                              verbose = TRUE)
                         
                         return(CD_likelihoods@posteriors) 
                       }
  
  
  }
  
  names(Bayseq_DE_df)<- names(DE_list)
  
  return(Bayseq_DE_df)
} 






#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




# Section 4 


###################################################### 

# 

######################################################

FDR_ROC_Inflect_table<-function(table, 
                        .sample = sample_id,
                        .transcript = Gene_ref,
                        .abundance = value,
                        .x = CAPRA_S, 
                        x_cord = seq(from = -5, to = 5, by = 0.5),
                        n_group,
                        method = "edgeR_quasi_likelihood"
                        ) {
  
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  .x = enquo(.x)
  
  require(dplyr)
  require(tibble)
  require(tidybulk)
  
  
  # Run tidybulk analysis - use tidybulk DE function
  DE_table<-tidybulk_DE(table = table,
                        .sample = !!.sample,
                        .transcript = !!.transcript,
                        .abundance = !!.abundance,    
                        .x = !!.x,
                        x_cord = x_cord,
                        n_group = n_group,
                        method = method)
  
  
  #################################################################
  
  ### Set up for converting alpha to inflection value 
  scale_design = function(df, formula){
    df %>%
      setNames(c("sample_idx", "(Intercept)", parse_formula(formula))) %>%
      tidyr::gather(cov, !!.abundance, -sample_idx) %>%
      group_by(cov) %>%
      mutate( value = ifelse( !grepl("Intercept", cov) & length(union(c(0,1), value)) != 2, scale(value), value )) %>%
      ungroup() %>%
      tidyr::spread(cov, value) %>%
      arrange(as.integer(sample_idx)) %>%
      select(`(Intercept)`, 
             one_of(parse_formula(formula)))
  }
  
  parse_formula <- function(fm) {
    if(attr(terms(fm), "response") == 1) stop("The formula must be of the kind \"~ covariates\" ")
    else as.character(attr(terms(fm), 
                           "variables"))[-1]
  }
  
  ####################################################################
  
  # Create tibble of covariate data 
  set.seed(100)
  alpha_X  = table %>% 
    select(!!.transcript) %>% 
    ungroup() %>% 
    distinct() %>% 
    sample_n(1)
  
  covariate_data = 
    inner_join(table, alpha_X) %>%
    select(one_of(parse_formula(~CAPRA_S)), !!.sample) %>%
    distinct()

    X =
    model.matrix(object = ~CAPRA_S, 
                 data = covariate_data) %>%
    as_tibble(rownames="sample_idx") %>%
    scale_design(~CAPRA_S) %>% 
    select(CAPRA_S)
  
    
    #########################################################################
  
  #Calculate estimated FDR using returned table from tidybulk DE
    # Because tidybulk returns a slightly different table depending on the method
    # estimating/ extracting FDR estimates is different for each method
  
  if (method == "limma_voom") {
    
    FDR_table<-DE_table %>% 
      select(!!.transcript,
             contains("P.Value"))  %>% 
      pivot_longer(!(!!.transcript),
                   names_to= "Contrast",
                   values_to = "PValue")%>%
      group_by(!!.transcript)  %>%
      filter(rank(-PValue) == 1) %>% 
      separate(Contrast, 
               c("name", "Contrast"), 
               sep = "___") 
    
  }
  
  else if(method == "DESeq2") {
    FDR_table<-DE_table %>% 
      select(!!.transcript,
             contains("pvalue"))  %>% 
      pivot_longer(!(!!.transcript),
                   names_to= "Contrast",
                   values_to = "PValue")%>%
      group_by(!!.transcript)  %>%
      filter(rank(-PValue) == 1) %>% 
      separate(Contrast, 
               c("name", "Contrast"), 
               sep = "___") 
    
  }
  else{
  FDR_table<-DE_table %>% 
    select(!!.transcript,
           contains("FDR"))  %>% 
    pivot_longer(!(!!.transcript),
                 names_to= "Contrast",
                 values_to = "FDR")%>%
    group_by(!!.transcript)  %>%
    filter(rank(-FDR) == 1) %>% 
    separate(Contrast, 
             c("name", "Contrast"), 
             sep = "___") 
  } 
  
    #############################################
    
  # Extract Pvalues to calculate True FDR
    
    #######################################
  if (method == "limma_voom") {
    
    Pval_table<-DE_table %>% 
      select(!!.transcript,
             contains("P.Value")) %>% 
      pivot_longer(!(!!.transcript), 
                   names_to= "Contrast", 
                   values_to = "PValue") %>% 
      group_by(!!.transcript) %>%
      filter(rank(-PValue) == 1) %>% 
      separate(Contrast, 
               c("name", "Contrast"), 
               sep = "___") 
    
  }
  
  else if (method == "DESeq2") { 
    Pval_table<-DE_table %>% 
      select(!!.transcript,
             contains("pvalue")) %>% 
      pivot_longer(!(!!.transcript), 
                   names_to= "Contrast", 
                   values_to = "PValue") %>% 
      group_by(!!.transcript) %>%
      filter(rank(-PValue) == 1) %>% 
      separate(Contrast, 
               c("name", "Contrast"), 
               sep = "___") 
  }
  
  else{
  Pval_table<-DE_table %>% 
    select(!!.transcript,
           contains("PValue")) %>% 
    pivot_longer(!(!!.transcript), 
                 names_to= "Contrast", 
                 values_to = "PValue") %>% 
    group_by(!!.transcript) %>%
    filter(rank(-PValue) == 1) %>% 
    separate(Contrast, 
             c("name", "Contrast"), 
             sep = "___") 
  }
  
    

  pval_range<-seq(from = 0, #Pval_table$PValue %>% min(), 
                  to = 1, #Pval_table$PValue %>% min(),
                  length.out = 1000)
  
  N_sig<-sapply(pval_range, 
                function(x) Pval_table %>% 
                            filter(PValue<=x) %>% 
                            nrow())
  
  N_sig_F<-sapply(pval_range, 
                  function(x) Pval_table %>% 
                    filter(grepl("V",!!.transcript)) %>% 
                    filter(PValue<=x) %>% 
                    nrow())
  
  True_FDR=N_sig_F/N_sig
  
  t_FDR_df<-data.frame(
    True_FDR = True_FDR
  )
  
  
  
  # ROC Curve
  
  # For each pvalue threshold
  # how many true positives
  TP<-sapply(pval_range, 
             function(x) Pval_table %>% 
               filter(grepl("X",!!.transcript)) %>% 
               filter(PValue<=x) %>% 
               nrow())
  
  # Total number of true positives
  P<-Pval_table %>% 
    filter(grepl("X",!!.transcript)) %>% 
    nrow()
  
  # For each pvalue threshold
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
  # Rank using logFC 
  
  if (method == "DESeq2") { 
    logFC<-DE_table %>% 
      select(!!.transcript, 
             contains("FoldChange"))  %>% 
      mutate_if(is.numeric, abs) %>%
      pivot_longer(!(!!.transcript), 
                   names_to= "Contrast", 
                   values_to = "FC") %>%
      group_by(!!.transcript) %>%
      filter(rank(-FC) == 1) %>% 
      separate(Contrast, 
               c("name", "contrast"), 
               sep = "___")
    
    }
  

  
  else {
  logFC<-DE_table %>% 
    select(!!.transcript, 
           contains("FC"))  %>% 
    mutate_if(is.numeric, abs) %>%
    pivot_longer(!(!!.transcript), 
                 names_to= "Contrast", 
                 values_to = "FC") %>%
    group_by(!!.transcript) %>%
    filter(rank(-FC) == 1) %>% 
    separate(Contrast, 
             c("name", "contrast"), 
             sep = "___")
  
  } 
  
  
  if (method == "DESeq2") {
    
    contrasts= sapply(2:n_group, 
             function(x) paste0("group", " ", LETTERS[x], "-", LETTERS[x-1])) } 
    
      
  
  else{
    
    contrasts<-sapply(2:n_group, 
                    function(x) 
                      paste0("group", LETTERS[x], "-", "group", LETTERS[x-1]))} 
  
  
  contrasts_names<-sapply(2:n_group, 
                          function(x) 
                            paste0(LETTERS[x-1], LETTERS[x]))
  
  
  jump<-make_groups(table, 
                    .sample = !!.sample,
                    .transcript = !!.transcript,
                    .abundance = !!.abundance,   
                    .x = !!.x,
                    x_cord = x_cord,
                    n_group = n_group)$contrast_values
  
  con_df<-data.frame(contrast = contrasts,
                     group = contrasts_names,
                     jump_est = jump)
  

  
  alpha<- table %>% 
    select(!!.transcript, slope, alpha) %>% 
    distinct() %>% 
    mutate(alpha = ifelse(slope>0,
                          -1*as.numeric(alpha), 
                          as.numeric(alpha))) 
  # set.seed(100)
  # alpha_X  = table %>% 
  #   select(!!.transcript) %>% 
  #   ungroup() %>% 
  #   distinct() %>% 
  #   sample_n(1)
  # 
  # alpha_X  = inner_join(table, alpha_X) %>%  
  #          select(!!.transcript, CAPRA_S) %>% 
  #          distinct() %>% 
  #   select(alpha = CAPRA_S)
  # 
  # 
  # 
  # alpha_X = data.frame(alpha_X %>% rename(CARPA_S = alpha), 
  #                      X %>% 
  #                        rename(inflect = CAPRA_S))
  # 
  # 
  # alpha = data.frame(alpha, 
  #                   alpha_X) %>% 
  #   select(!!.transcript, 
  #          CAPRA_S, 
  #          inflect) %>% 
  #   distinct()
  # 
  # return(alpha) } 
  
  sample_size = table %>% 
    select(sample_size) %>% 
    distinct() %$%
    sample_size
  
  Inflect_df<-left_join(alpha,
                left_join(logFC,
                          con_df)) %>% 
              mutate(n_groups = n_group) %>% 
              mutate(sample_size = sample_size) 
    
  
  return(list(
    pred_FDR = FDR_table,
    true_FDR = t_FDR_df,
    ROC = ROC_df,
    Inflect_est = Inflect_df
  ))
}

# table_run = list()
#   
# for (i in 3:21) {
#   table_run[[i-2]]=FDR_ROC_Inflect_table(rbind( a %>% rename(True_Slope = slope, Sample = sample_id),
#                                             a_1 %>% 
#                                  mutate(Gene_number = paste0(Gene_number, "neg")) %>% 
#                                                              rename(True_Slope = slope, Sample = sample_id) ),
#                         i)$Inflect_est 
#     
# 
#     
# }


# Load simulation data

sim_df_slope_sample_size <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Article_Sections/Part_1_TABI_Analysis/sim_df_slope_sample_size.rds")


library(purr)

# Run edgeR analysis



# If package box is not installed, then install package 
if (
  !("box" %in% installed.packages()[ , "Package"])) {
  install.packages("box")
}

# Load functions from Functions_Simulate_RNA_seq_df using box 
options(box.path = setwd("~/TABI/dev/Article_Sections"))   


#' @export
box::use(./Functions_Simulate_RNAseq_df) 



par_vector= expand.grid(c(0.25, -0.25, 0.5, -0.5, 1, -1, 3, -3),
                        c(11,22,55, 77))


library(furrr)

plan(multisession, workers = 15)

sim_df_slopes_sample_size = purrr::pmap_dfr(.l = par_vector, 
                                              ~Functions_Simulate_RNAseq_df$sigmoidal_sim_df(
                                              n_true_tests = 5,
                                              n_false_tests = 0,
                                              beta = .x, 
                                              k = 8,
                                              A = 3,
                                              sample_size = .y,
                                              disp_size = 0.85, 
                                              alpha = "all",
                                              covar =  seq(from = -2.5, to = 2.5, by = 0.5)
                                                                                             )
                                                                                              ) 



by_Gene_ref = sim_df_slopes_sample_size %>% 
  mutate(sample_id = as.factor(sample_id)) %>% 
  split(.$Gene_ref)

n<- 15

# Set up parallelisation clusters
plan(multisession, 
     workers = 15)


TABI_sample_size_slope = by_Gene_ref %>% 
  future_map(
    ~TABI::TABI_glm(
      .data = .,
      ~CAPRA_S,
      .abundance = value, # simulated RNA seq value
      .sample = sample_id, #sample column is sample_id
      .transcript = Gene_ref, # column of transcript / gene ids 
      chains = 4,
      cores = 1,
      iter = 2000, # total number of iterations
      warmup = 1000
    ) %$% 
      fit, .progress = TRUE) 

gene_names = names(TABI_sample_size_slope)

TABI_list = lapply(1:length(gene_names),
       function(x) TABI_sample_size_slope[[x]] %>% 
                   rstan::summary() %$% 
                   summary %>% 
                   as.data.frame() %>% 
                   tibble::rownames_to_column("Pars") %>% 
                   filter(Pars == "inflection[1]") %>% 
                   mutate(Gene_ref = gene_names[x]))

TABI_df = data.table::rbindlist(TABI_list)

inner_join(TABI_df, 
           sim_df_slopes_sample_size %>% 
             select(Gene_ref, inflect, slope, sample_size) %>% 
             distinct()) %>% 
  rowwise() %>% 
  mutate(dif = abs(mean - inflect)) %>% 
  select(dif, slope, sample_size) %>% 
  ggplot(aes(y=dif)) + 
  facet_grid(rows = vars(sample_size), cols = vars(slope)) + 
  geom_boxplot()

plot_rsim_df()
inflect_list = list()

method = c("edgeR_quasi_likelihood", "DESeq2", "limma+voom")


sim_df_slopes_sample_size %>% 
  filter(sample_size == 77) %>% 
  tidybulk_DE(., n_group = 3, 
                        x_cord = seq(from = -2.5, to = 2.5, by = 0.5))

run_edgeR =lapply(
                  c(11,22,55, 77), 
                       function(y) lapply(3:11, 
                                          function(x) { if(x == 11 & y == 11) {}
                                              else{
                                              print(paste(x,y))
                                          
                                              sim_df_slopes_sample_size %>% 
                                              filter(sample_size == y) %>% 
                                              FDR_ROC_Inflect_table(., 
                                                                       n_group = x, 
                                                                    x_cord = seq(from = -2.5, to = 2.5, by = 0.5)
                                                                    ) %$% 
                                                                         Inflect_est %>% 
                                                                          mutate(N_groups = x) } } 
                                                                                      ) 
                                                                                       )
run_DESeq2 = lapply(c(11,22,55, 77), 
function(y) lapply(3:11, 
                   function(x) { if(x == 11 & y == 11) {}
                     else{
                       print(paste(x,y))
                       
                       sim_df_slopes_sample_size %>% 
                         filter(sample_size == y) %>% 
                         FDR_ROC_Inflect_table(., 
                                               n_group = x, 
                                               method = "DESeq2",
                                               x_cord = seq(from = -2.5, to = 2.5, by = 0.5)
                         ) %$% 
                         Inflect_est %>% 
                         mutate(N_groups = x) } } 
) 
)


run_limma_voom = lapply(c(11,22,55, 77), 
                    function(y) lapply(3:11, 
                                       function(x) { if(x == 11 & y == 11) {}
                                         else{
                                           print(paste(x,y))
                                           
                                           sim_df_slopes_sample_size %>% 
                                             filter(sample_size == y) %>% 
                                             FDR_ROC_Inflect_table(., 
                                                                   n_group = x, 
                                                                   method = "limma_voom",
                                                                   x_cord = seq(from = -2.5, to = 2.5, by = 0.5)
                                             ) %$% 
                                             Inflect_est %>% 
                                             mutate(N_groups = x) } } 
                    ) 
)


run_DESeq2 =lapply(3:11, 
            function(x)
              sim_df_slopes_sample_size %>% 
              filter(sample_size == 77) %>% 
              FDR_ROC_Inflect_table(., n_group = x, method = "DESeq2",
                                    x_cord = seq(from = -2.5, to = 2.5, by = 0.5)) %$% Inflect_est) 

run_limma_voom = lapply(3:11, 
                        function(x)
                          sim_df_slopes_sample_size %>% 
                          filter(sample_size == 77) %>% 
                          FDR_ROC_Inflect_table(., n_group = x, method = "limma_voom",
                                                x_cord = seq(from = -2.5, to = 2.5, by = 0.5)) %$% Inflect_est) 


bind_rows(run_edgeR) %>% mutate(method = "edgeR")
bind_rows(run_DESeq2) %>% mutate(method = "DESeq2")
bind_rows(run_limma_voom) %>% mutate(method = "limma+voom")

%>% 
sim_df_slopes_sample_size %>% 
  filter(sample_size == 11) %>% 
  FDR_ROC_Inflect_table(., n_group = 10, 
              x_cord = seq(from = -2.5, to = 2.5, by = 0.5)) %$% Inflect_est


sim_df_slopes_sample_size %>% 
  filter(sample_size == 77) %>% 
  FDR_ROC_Inflect_table(., n_group = 4, 
              x_cord = seq(from = -2.5, to = 2.5, by = 0.5))

sim_df_slopes_sample_size %>% 
  filter(sample_size == 77) %>% 
  tidybulk_DE(., n_group = 11, method = "DESeq2",
                        x_cord = seq(from = -2.5, to = 2.5, by = 0.5))
  

library(magrittr)

sim_df_slopes_sample_size %>% 
  filter(sample_size == 11) %>% 
  FDR_ROC_Inflect_table(., n_group = 3, 
                        x_cord = seq(from = -2.5, to = 2.5, by = 0.5)) %$% Inflect_est

for(i in c(3:21)) {

inflect_list[[i]] = sim_df_slopes_sample_size %>% 
                    split(.$sample_size) %>% 
                    purrr::map(~FDR_ROC_Inflect_table(., 
                                                      n_group = i, 
                                                      x_cord = seq(from = -2.5, to = 2.5, by = 0.5))$Inflect_est) 

} 
  



sim_df_slopes_sample_size  %>% 
  filter(sample_size == 11) %>% 
  make_groups(., 
                          n_group = 3, 
                          x_cord = seq(from = -2.5, to = 2.5, by = 0.5))



sim_df %>% 
  split(.$sample_size) %>% 
  purrr::map(~tidybulk_DE(., 
                                    n_group = 3)$Inflect_est) 















# Plotting of errors 

# 
# tidybulk_DE(sim_slope_sample_size_df, n_group = 3)
# 
# library(foreach)
# plot_table=foreach(i=c(1:19)) %do% {
#   table_run[[i]] %>% mutate(N_groups = i+2)
#   
# }
# test = rbind(sigmoidal_sim_df(100, 0, 1, 3, 5, 21*3, 1),
#              sigmoidal_sim_df(100, 0, -1, 3, 5, 21*3, 1))
# 
# 
# edgeR_inflection=lapply(2:21, 
#        function(x) (FDR_ROC_Inflect_table(test, 
#                                          n_group =x, 
#                                          method = "edgeR_quasi_likelihood"))$Inflect_est %>% 
#          mutate(N_groups = x)) 
# 
# DESeq2_inflection=lapply(2:21, 
#                         function(x) (FDR_ROC_Inflect_table(test, 
#                                                            n_group =x, 
#                                                            method = "DESeq2"))$Inflect_est %>% 
#                           mutate(N_groups = x)) 
# 
# 
# limma_voom_inflection=lapply(2:21, 
#                          function(x) (FDR_ROC_Inflect_table(test, 
#                                                             n_group =x, 
#                                                             method = "limma_voom"))$Inflect_est %>% 
#                            mutate(N_groups = x)) 
# 
# library(gridExtra)

data = rbind(
  bind_rows(run_edgeR) %>%
                 mutate(Method = "edgeR"),
  bind_rows(run_limma_voom) %>%
    mutate(Method = "limma + voom"),
  bind_rows(run_DESeq2) %>%
  mutate(Method = "DESeq2")) %>%
  rowwise() %>% 
  mutate(Value = paste0(Method, n_groups)) %>% 
  mutate(alpha = alpha/slope) %>% 
  mutate(error = abs(jump_est-alpha))


TABI_error= inner_join(TABI_df, 
           sim_df_slopes_sample_size %>% 
             select(Gene_ref, inflect, slope, sample_size) %>% 
             distinct()) %>% 
  rowwise() %>% 
  mutate(error = abs(mean - inflect)) %>% 
  mutate(Method = "TABI") %>% 
  mutate(n_groups = 0)

 
data = full_join(data, TABI_error)
# 
# 
# 
# 
# 
# 


data_bigger = data %>% ungroup() %>% filter(slope > 0)

data_smaller = data %>% ungroup() %>% filter(slope < 0)

# 

  data = bocplot_eror_df 
data_bigger %>% 
  filter(Method == "TABI") %>% 
  mutate(n_groups = 0) %>% 
  ggplot(aes(y=error,
             x = n_groups,
             group = Value,
             fill = as.factor(Method))) + 
  geom_boxplot() + 
  facet_grid(rows = vars(sample_size), cols = vars(slope))
  
data_bigger %>% 
  group_by(Value, Method, n_groups, sample_size, slope) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(slope = as.factor(slope)) %>% 
  filter(Method == "TABI") %>% 
  ggplot(aes(y=mean_error, x=sample_size, group = slope, col = slope)) + 
  geom_line()

p1 = data_bigger %>% 
  group_by(Value, Method, sample_size, slope) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  filter(Method == "TABI") %>% 
  ggplot(aes(y=mean_error, x=slope, group = sample_size, col = sample_size)) + 
  geom_line() + 
  labs(y = "", x = "", title = "TABI", col = "Sample Size") + 
  ylim(0,2.5) + 
  custom_theme_2


data_bigger %>% 
  group_by(Method, n_groups, sample_size) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  filter(Method == "limma + voom") %>% 
  ggplot(aes(y=mean_error, x=n_groups, group = sample_size, col = sample_size)) + 
  geom_line() + 
  labs(y = "", x = "", title = "TABI", col = "Sample Size") + 
  ylim(0,2.5) + 
  custom_theme_2

p11 = data_bigger %>% 
  group_by(Method, n_groups, slope) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(slope = as.factor(slope)) %>% 
  filter(Method == "edgeR") %>% 
  ggplot(aes(y=mean_error, x=n_groups, group = slope, col = slope)) + 
  geom_line() + 
  labs(y = "", x = "", title = "edgeR", col = "Slope") + 
  ylim(0,2.5) + 
  custom_theme_2

p21 = data_bigger %>% 
  group_by(Method, n_groups, sample_size) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  filter(Method == "edgeR") %>% 
  ggplot(aes(y=mean_error, x=n_groups, group = sample_size, col = sample_size)) + 
  geom_line() + 
  labs(y = "", x = "", title = "edgeR", col = "Sample Size") + 
  ylim(0,2.5) + 
  custom_theme_2


p12 = data_bigger %>% 
  group_by(Method, n_groups, slope) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(slope = as.factor(slope)) %>% 
  filter(Method == "DESeq2") %>% 
  ggplot(aes(y=mean_error, x=n_groups, group = slope, col = slope)) + 
  geom_line() + 
  labs(y = "", x = "", title = "DESeq2", col = "Slope") + 
  ylim(0,2.5) + 
  custom_theme_2

p22 = data_bigger %>% 
  group_by(Method, n_groups, sample_size) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  filter(Method == "DESeq2") %>% 
  ggplot(aes(y=mean_error, x=n_groups, group = sample_size, col = sample_size)) + 
  geom_line() + 
  labs(y = "", x = "", title = "DESeq2", col = "Sample Size") + 
  ylim(0,2.5) + 
  custom_theme_2

p13 = data_bigger %>% 
  group_by(Method, n_groups, slope) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(slope = as.factor(slope)) %>% 
  filter(Method == "limma + voom") %>% 
  ggplot(aes(y=mean_error, x=n_groups, group = slope, col = slope)) + 
  geom_line() + 
  labs(y = "", x = "", title = "limma + voom", col = "Slope") + 
  ylim(0,2.5) + 
  custom_theme_2

p23 = data_bigger %>% 
  group_by(Method, n_groups, sample_size) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  filter(Method == "limma + voom") %>% 
  ggplot(aes(y=mean_error, x=n_groups, group = sample_size, col = sample_size)) + 
  geom_line() + 
  labs(y = "", x = "", title = "limma + voom", col = "Sample Size") + 
  ylim(0,2.5) + 
  custom_theme_2

ggarrange(p11, p12, p13, p21, p22, p23, nrow =2, ncol = 3)


p2 = data_bigger %>% 
  group_by(Method, sample_size, slope) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  filter(Method == "DESeq2") %>% 
  ggplot(aes(y=mean_error, x=slope, group = sample_size, col = sample_size)) + 
  geom_line() +
  labs(y = "", x = "", title = "DESeq2", col = "Sample Size") + 
  ylim(0,2.5) + 
  custom_theme_2


p3 = data_bigger %>% 
  group_by(Method,sample_size, slope) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  filter(Method == "edgeR") %>% 
  ggplot(aes(y=mean_error, x=slope, group = sample_size, col = sample_size)) + 
  geom_line() + 
  labs(y = "", x = "", title = "edgeR", col = "Sample Size") + 
  ylim(0,2.5) + 
  custom_theme_2

p4 = data_bigger %>% 
  group_by(Method, sample_size, slope) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  filter(Method == "limma + voom") %>% 
  ggplot(aes(y=mean_error, x=slope, group = sample_size, col = sample_size)) + 
  geom_line() + 
  labs(y = "", x = "", title = "limma + voom", col = "Sample Size") + 
  ylim(0,2.5) + 
  custom_theme_2




bottom_plot = ggpubr::annotate_figure(ggpubr::ggarrange(p1, p2,p3, p4, 
                                          ncol = 4, 
                                          common.legend = TRUE, 
                                          legend = "bottom"), 
bottom = text_grob("Slope", vjust = -4, size = 10),
                                left = text_grob("Mean Error", rot = 90, size = 10))

data_bigger %>% 
  group_by(Value, Method, n_groups, sample_size, slope) %>% 
  summarise(mean_error = mean(error)) %>% 
  ungroup() %>%
  mutate(slope = as.factor(slope)) %>% 
  filter(Method == "DESeq2") %>% 
  ggplot(aes(y=mean_error, x=sample_size, group = slope, col = slope)) + 
  geom_point()

top_plot_data = data %>% 
  filter(sample_size == 22, slope == 0.5) 

top_plot = top_plot_data %>% 
  ggplot(aes(y=error,
             x = n_groups,
             group = Value,
             fill = Method)) +
  geom_boxplot() +
  #facet_grid(rows = vars(sample_size), cols = vars(slope)) + 
  #stat_summary(fun=mean, geom="line", aes(group=1))  +
  stat_summary(fun=mean, geom="line", aes(x=0, col = Method), data= top_plot_data %>% filter(Method == "TABI")) + 
  stat_summary(fun=mean, geom="line", aes(group=1, col = Method), data= top_plot_data %>% filter(Method == "edgeR"))  +
  stat_summary(fun=mean, geom="line", aes(group=1, col = Method), data= top_plot_data %>% filter(Method == "limma + voom"))  +
  stat_summary(fun=mean, geom="line", aes(group=1, col = Method), data= top_plot_data %>% filter(Method == "DESeq2"))  +
  #stat_summary(fun=mean, geom="point", aes(group = Value), position = position_dodge(0.75), shape = 18, size = 2) +
  labs(x = "Number of Groups",
       y = expression(atop("Error", "(Absolute Difference between True and Estimated Time of Change"))) +
  scale_fill_discrete(name =  "Method") +
  scale_linetype_discrete(name =  "Method") + 
  #facet_grid(rows = vars(sample_size), cols = vars(slope)) + 
  stat_summary(fun="mean", geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group = Value), linetype = "dashed") + 
  custom_theme

#labs(title = "Error in Estimating Most Rapid Point of Change to Transcript Abudance") + 

# top_plot = annotate_figure(top_plot, 
#                 top =  text_grob("Slope", vjust = 2, size = 10),
#                 right = text_grob("Sample Size", rot = 90, size = 10, vjust = - 0.5))

ggarrange(top_plot,
          bottom_plot,
          ncol = 1,
          heights = c(2,1),
          labels = c("A", "B"))
  
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dittoSeq")
  

  friendly_cols <- dittoSeq::dittoColors()
#friendly_cols <- brewer.pal(n =4, name = "Dark2")
friendly_cols_fil<-brewer.pal(n =4, name = "Dark2")


# Set theme
custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols,   na.value = "black"),
    scale_color_manual(values = friendly_cols,   na.value = "black"),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 9),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  )


custom_theme_2 <-
  list(
    scale_fill_manual(values = friendly_cols_fil,   na.value = "black"),
    scale_color_manual(values = friendly_cols_fil,   na.value = "black"),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 9),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  )

ggsave("Figure_2_v4.pdf",
       width = 183,
       height = 238,
       units = "mm"
       ) 

ggsave("Figure_2_v4.pdf",
       width = 183,
       height = 203,
       units = "mm"
) 


saveRDS(data, file = "boxplot_error_df.rds")


data_smaller %>%
  ggplot(aes(y=error,
             x = n_groups,
             group = Value,
             fill = as.factor(Method))) +
  geom_boxplot() +
  facet_grid(rows = vars(sample_size), cols = vars(slope)) + 
  #stat_summary(fun=mean, geom="line", aes(group=1))  +
  stat_summary(fun=mean, geom="line", aes(x=0, col = Method), data= data_smaller %>% filter(Method == "TABI")) + 
  stat_summary(fun=mean, geom="line", aes(group=1, col = Method), data= data_smaller %>% filter(Method == "edgeR"))  +
  stat_summary(fun=mean, geom="line", aes(group=1, col = Method), data= data_smaller %>% filter(Method == "limma + voom"))  +
  stat_summary(fun=mean, geom="line", aes(group=1, col = Method), data= data_smaller %>% filter(Method == "DESeq2"))  +
  #stat_summary(fun=mean, geom="point", aes(group = Value), position = position_dodge(0.75), shape = 18, size = 2) +
  labs(x = "Number of Groups",
       y = "Error (Absolute Difference between True Time of Change and Estimated Time of Change") +
  scale_fill_discrete(name =  "Method") +
  scale_linetype_discrete(name =  "Method") + 
  facet_grid(rows = vars(sample_size), cols = vars(slope)) + 
  stat_summary(fun="mean", geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group = Value), linetype = "dashed")




# 
# 
# 
# data.table::rbindlist(DESeq2_inflection) %>% 
#   mutate(error = abs(jump_est-alpha)) %>% 
#   ggplot(aes(y=error, x = N_groups, group = N_groups, fill = as.factor(N_groups))) + 
#   geom_boxplot()
# edgeR_inflection = list()



