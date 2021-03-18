
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Intention of Script 

#Create Functions which run DESeq2 / edgeR analysis 
#And compare their utility to TABI 
#Based on FDR and ROC curves


#(Section 1)
#Create a function which takes Posterior_df from TABI and simulated dataset
#Returns PEP Table, FDR statistics and ROC Statistis

#(Section 2)
#Create a function which runs edgeR analysis on simulated datasets 

#(Section 3)
#Create a function which performs DESeq2 analysis on simulated datasets

#(Section 4)
#Create a function which performs Bayseq analysis on simulated datasets

#(Section 5)
#Create function which 
#Returns FDR statistics and ROC statistics 
#For pvalues returned from DESeq2 and edgeR analysis 



##         ##         ##         ##         ##         ##        ##




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


#Section 1 

#Create a function which takes Posterior_df from TABI and simulated dataset
#Returns PEP Table, FDR statistics and ROC Statistics


TABI_fdr_roc_stats<-function(Posterior_df, 
                             simulated_data) {
  
  require(dplyr)
  
  #Calculate probabilities / PEP for each test / gene
  
  PEP_table<-Posterior_df %>% 
    filter(covariate_idx == 1) %>% 
    group_by(gene) %>% #For each gene / test
    summarise(prob_less = sum(.value<0)/n(), #Calculate the probability that the slope is less than zero, as the proportion of calculated values < 0
              prob_more = sum(.value>0)/n(),  #Calculate the probability that the slope is greater than zero, as the proportion of calculated values > 0
              mean = mean(.value),   #Calculate the mean calculated slope as the mean of all calculated slope values
              n = n()) %>% 
    rowwise %>% 
    mutate(prob = #Calculate the probability the slope is not zero as 2 times the min of prob > 0, prob < 0 
             2*min(
               unlist(
                 c(prob_less, prob_more)
                 )
               )
           ) %>% 
    mutate(PEP = prob) %>%
    ungroup %>%
    arrange(PEP) %>% 
    mutate(qvalue = cummean(PEP)) #Calculating FDR as per http://varianceexplained.org/r/bayesian_fdr_baseball/
  
  #Theoretical FDR (Using local FDR / PEP definition)
  
  #Calculate FDR using method in article http://varianceexplained.org/r/bayesian_fdr_baseball/
  
  #Check distribution of qvalues vs PEP
  
  PEP_table<-inner_join(PEP_table, 
             simulated_data %>% 
               select(Gene_number, Null_test) %>% 
               dplyr::rename(gene = Gene_number) %>% 
               dplyr::distinct())
  
  
  # Calculate number of tests declared significant at each level of FDR (expected number of false discoveries)
  
  N_sig_expFDR<-sapply(seq(from =0, to = 1, by = 1/1000), # Each FDR 
                       function(x) PEP_table %>% 
                         filter(qvalue<x) %>% 
                         nrow() 
  )
  
  # True proportion of False Discoveries at each qvalue
  
  Null_tests_PEP_table<-inner_join( #Filtering PEP table so that it only has Null results
    PEP_table,
    simulated_data %>% 
      filter(Null_test == "TRUE") %>% 
      select(Gene_number) %>% 
      dplyr::rename(gene = Gene_number) %>% 
      dplyr::distinct()
  )
  
  Null_tests_sig<-sapply(seq(from =0, to = 1, by = 1/1000), # Each FDR 
                         function(x) Null_tests_PEP_table %>% 
                           filter(qvalue<x) %>% 
                           nrow()
                         
  )
  
  # Actual proportion of false discoveries at each qvalue
  
  True_FD_TABI<-Null_tests_sig/N_sig_expFDR
  True_FD_TABI[1]<-0
  
  # Plot FDR Curves (Predicted, and Actual)
  
  FDR_df<-data.frame(
    FDR_pred = seq(from =0, to = 1, by = 1/1000), 
    FDR_true = True_FD_TABI,
    N_sig = N_sig_expFDR
  )
  
  
  #Create an ROC Curve (for simulated data, previously run on TABI)
  
  # Number of postives at each PEP
  
  True_tests_PEP_table<-inner_join( #Filtering PEP table so that it only has Null results
    PEP_table,
    simulated_data %>% 
      filter(Null_test == "FALSE") %>% 
      select(Gene_number) %>% 
      dplyr::rename(gene = Gene_number) %>% 
      dplyr::distinct()
  )
  
  TP<-sapply(
    seq(from=0, to=1, by=1/1000), 
    function(x) {
      True_tests_PEP_table %>% 
        filter(PEP <x) %>% 
        nrow()
    }
  )
  
  # Total number of Positives / True Tests ( = 100)
  
  P<-simulated_data %>% 
    filter(Null_test == "FALSE") %>% 
    select(Gene_number) %>% 
    dplyr::distinct() %>% 
    nrow()
  
  # TPR True Positive Rate 
  
  TPR<-TP/P
  
  #False Positive Rate 
  
  # True Negatives 
  TN<-sapply(
    seq(from=0, to=1, by=1/1000), 
    function(x) Null_tests_PEP_table %>% 
      filter(PEP>=x) %>% 
      nrow()
  )
  
  #Total numer of Negative / Null Tests ( = 900)
  
  N<-simulated_data %>% 
    filter(Null_test == "TRUE") %>% 
    select(Gene_number) %>% 
    dplyr::distinct() %>% 
    nrow()
  
  
  #False Postive Rate
  
  FPR<-1 - TN/N
  
  ROC_df<- data.frame(
    TPR = TPR,
    FPR = FPR
  )
  
  ROC_FDR<-list(FDR =FDR_df, ROC = ROC_df) %>% 
    bind_cols()
  
  return(list(PEP_table = PEP_table, ROC_FDR_table = ROC_FDR))
}






##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 


#Section 2 

#Function for edgeR analysis (runs edgeR analysis, returns results (Pvalues) from edgeR analysis with 2 groups, 3 groups, 4 groups)

edgeR_group_analysis<- function(simulated_data) {
  
  require(dplyr)
  require(tibble)
  require(tidyr)
  require(edgeR)
  require(qvalue)

  test<- simulated_data
  
  #     #    #     #     #     #     #    #     #     #     #      #      #      #    #  
  
  # 2 Groups 
  
  # P1: 2 Group EdgeR Analysis 
  
  edgeR_df_A_2<-test %>% 
    filter(CAPRA_S <0) %>% 
    select(Gene_number, value, sample) %>% 
    pivot_wider(names_from = sample, values_from = value) %>% 
    column_to_rownames(var= "Gene_number")
  
  edgeR_df_B_2<-test %>% 
    filter(CAPRA_S >=0) %>% 
    select(Gene_number, value, sample) %>% 
    pivot_wider(names_from = sample, values_from = value) %>% 
    column_to_rownames(var= "Gene_number")
  
  edgeR_df_2<-cbind(
    edgeR_df_A_2,
    edgeR_df_B_2
  ) 
  
  n_A_2<-ncol(edgeR_df_A_2)
  n_B_2<-ncol(edgeR_df_B_2)
  
  group_2<-c(rep("A", n_A_2), 
             rep("B", n_B_2))
  
  
  y_2 <- DGEList(counts = edgeR_df_2, 
                 group = group_2)
  
  y_2 <- calcNormFactors(y_2)

  y_2 <- estimateDisp(y_2)

  et_2 <- exactTest(y_2)
 
  edgeR_2group<-et_2$table %>% 
    rownames_to_column(var = "Gene_number") %>% 
    dplyr::rename(Pval_final = PValue) %>% 
    select(Gene_number, Pval_final) %>% 
    as_tibble()

  #     #    #     #     #     #     #    #     #     #     #      #      #      #    #  
  
  # 3 Group EdgeR analysis 
  
  edgeR_df_A_3<-test %>% 
    filter(CAPRA_S <2*-1) %>% 
    select(Gene_number, value, sample) %>% 
    pivot_wider(names_from = sample, values_from = value) %>% 
    column_to_rownames(var= "Gene_number")
  
  edgeR_df_B_3<-test %>% 
    filter(CAPRA_S >=(2*-1)&CAPRA_S<(2)) %>% 
    select(Gene_number, value, sample) %>% 
    pivot_wider(names_from = sample, values_from = value) %>% 
    column_to_rownames(var= "Gene_number")
  
  edgeR_df_C_3<-test %>% 
    filter(CAPRA_S >=2) %>% 
    select(Gene_number, value, sample) %>% 
    pivot_wider(names_from = sample, values_from = value) %>% 
    column_to_rownames(var= "Gene_number")
  
  edgeR_df_3<-cbind(
    edgeR_df_A_3,
    edgeR_df_B_3,
    edgeR_df_C_3
  ) 
  
  
  n_A_3<-ncol(edgeR_df_A_3)
  n_B_3<-ncol(edgeR_df_B_3)
  n_C_3<-ncol(edgeR_df_C_3)
  
  group_3<-c(rep("A", n_A_3), 
             rep("B", n_B_3),
             rep("C", n_C_3))
  
  
  y_3 <- DGEList(edgeR_df_3)
  
  y_3 <- calcNormFactors(y_3)
  
  design_3 <- model.matrix(~group_3)
  
  
  y_3 <- estimateDisp(y_3, 
                      design_3)
  
  fit_3 <- glmQLFit(y_3, 
                design_3)
  
  lrt_early_3 <- glmQLFTest(fit_3, 
                      coef = 2)
  
  lrt_late_3 <- glmQLFTest(fit_3, 
                       coef = 3)
  
  
  Pval_3_group<-(inner_join(lrt_late_3$table %>% 
                              select(PValue) %>% 
                              rownames_to_column(var = "Gene_number") %>% 
                              dplyr::rename(Pval_late = PValue),
                            lrt_early_3$table %>% select(PValue) %>% 
                              rownames_to_column(var = "Gene_number") %>% 
                              dplyr::rename(Pval_early = PValue)))
  

  
  Pval_3_group_total<-Pval_3_group %>% 
    group_by(Gene_number) %>% 
    summarise(Pval_final = 
                ifelse(Pval_early<0.05&Pval_late>=0.05,
                       Pval_early,
                       ifelse(Pval_early>=0.05&Pval_late<0.05,
                              Pval_late,
                              sample(c(Pval_early, Pval_late), 1)
                              )
                       ))
              

    #     #    #     #     #     #     #    #     #     #     #      #      #      #    #  
              
    # 4 group analysis 
    
    edgeR_df_A_4<-test %>% 
      filter(CAPRA_S <(2.5*-1)) %>% 
      select(Gene_number, value, sample) %>% 
      pivot_wider(names_from = sample, values_from = value) %>% 
      column_to_rownames(var= "Gene_number")
    
    edgeR_df_B_4<-test %>% 
      filter(CAPRA_S >=(2.5*-1)&CAPRA_S<0) %>% 
      select(Gene_number, value, sample) %>% 
      pivot_wider(names_from = sample, values_from = value) %>% 
      column_to_rownames(var= "Gene_number")
    
    edgeR_df_C_4<-test %>% 
      filter(CAPRA_S >=0&CAPRA_S<2.5) %>% 
      select(Gene_number, value, sample) %>% 
      pivot_wider(names_from = sample, values_from = value) %>% 
      column_to_rownames(var= "Gene_number")
    
    edgeR_df_D_4<-test %>% 
      filter(CAPRA_S >=2.5) %>% 
      select(Gene_number, value, sample) %>% 
      pivot_wider(names_from = sample, values_from = value) %>% 
      column_to_rownames(var= "Gene_number")
    
    edgeR_df_4<-cbind(
      edgeR_df_A_4,
      edgeR_df_B_4,
      edgeR_df_C_4,
      edgeR_df_D_4
    ) 
    
    n_A_4<-ncol(edgeR_df_A_4)
    n_B_4<-ncol(edgeR_df_B_4)
    n_C_4<-ncol(edgeR_df_C_4)
    n_D_4<-ncol(edgeR_df_D_4)
    
    group_4<-c(rep("A", n_A_4), 
               rep("B", n_B_4),
               rep("C", n_C_4),
               rep("D", n_D_4))
    
    y_4 <- DGEList(edgeR_df_4)
    
    y_4 <- calcNormFactors(y_4)
    
    design_4 <- model.matrix(~group_4)
    
    
    y_4 <- estimateDisp(y_4, 
                        design_4)
    
    fit_4 <- glmQLFit(y_4, 
                    design_4)
    
    lrt_early_4 <- glmQLFTest(fit_4, 
                          coef = 2)
    
    lrt_mid_4 <- glmQLFTest(fit_4, 
                         coef = 3)  
    
    lrt_late_4 <- glmQLFTest(fit_4, 
                         coef = 4)   
    
    Pval_4_group<-inner_join(
      (
        inner_join(
          lrt_late_4$table %>% 
                                select(PValue) %>% 
                                rownames_to_column(var = "Gene_number") %>% 
                                dplyr::rename(Pval_late = PValue),
                              lrt_early_4$table %>% select(PValue) %>% 
                                rownames_to_column(var = "Gene_number") %>% 
                                dplyr::rename(Pval_early = PValue)
          )
        ),
                              lrt_mid_4$table %>% select(PValue) %>% 
                              rownames_to_column(var = "Gene_number") %>% 
                              dplyr::rename(Pval_mid = PValue)
      )
      
  
  Pval_4_group_total<-Pval_4_group %>% 
    group_by(Gene_number) %>% 
    summarise(Pval_final = 
                ifelse(Pval_early<0.05&Pval_late>=0.05&Pval_mid>=0.05,
                       Pval_early,
                       ifelse(Pval_early>=0.05&Pval_late<0.05&Pval_mid>=0.05,
                              Pval_late,
                        ifelse(Pval_early>=0.05&Pval_late>=0.05&Pval_mid<0.05,
                                 Pval_mid,
                        ifelse(Pval_early<0.05&Pval_late<0.05&Pval_mid>=0.05,
                               sample(
                                 c(
                                   Pval_early, 
                                   Pval_late),
                                 1), 
                        ifelse(Pval_early>=0.05&Pval_late<0.05&Pval_mid<0.05,
                               sample(
                                 c(
                                   Pval_mid, 
                                   Pval_late),
                                 1), 
                        ifelse(Pval_early<0.05&Pval_late>=0.05&Pval_mid<0.05,
                               sample(
                                 c(
                                   Pval_mid, 
                                   Pval_early),
                                 1),     
                          sample(c
                                 (
                                   Pval_early,
                                   Pval_mid,
                                   Pval_late),
                                 1)
                          )
                        )
                        )              
                        )         
                        )      
                       )
                )
      
  return(list(
    Group_2_Pvalue = edgeR_2group,
    Group_3_Pvalue = Pval_3_group_total,
              Group_4_Pval =Pval_4_group_total))        
}









##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 









#Section 3

#Create a function which performs DESeq2 analysis on simulated datasets


# Create a function which runs DE analysis with 
# 2, 3, 4 groups (to compare "when" DE occurs)
# From a simulated dataset 


# DESeq2 

{
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
} 


#Function which takes simulate dataset 
# with CAPRA_S from -5 to 5 
# and uses DESeq2 to perform DE analysis 
# With 2, 3 and 4 groups 
# returns the corresponding pvalue (unadj)

DESeq_analysis<-function(simulated_data) {
  
  
  {
    require(dplyr)
    require(tibble)
    require(tidyr)
    require(DESeq2)
    require(foreach)
  }
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  
  
  #Set up Read Count Data
  # Convert simulated data in form for DESeq2 analysis 
  # Data must be in the form – each column is a replicate / condition / CAPRA_S value
  # Each row is a transcript / gene 
  
  countData <-simulated_data%>% 
    select(Gene_ref, Sample, value) %>% 
    mutate(Sample = paste0("ID_", Sample)) %>% 
    tidyr::pivot_wider(names_from = Sample, values_from = value) %>%
    tibble::column_to_rownames(var = "Gene_ref") %>% 
    as.matrix()
  
  #Set up column data
  #Rows in column data must correspond columns in count matrix
  
  coldata<-simulated_data%>% 
    select(Sample, CAPRA_S) %>% 
    dplyr::distinct() %>% 
    mutate(Sample = paste0("ID_", Sample)) %>% 
    mutate(CAPRA_S = as.numeric(CAPRA_S)) %>% 
    mutate(Two_group = ifelse(CAPRA_S<0, "A", "B")) %>% 
    mutate(Three_group = ifelse(CAPRA_S<=2*-1, "A", ifelse(CAPRA_S>=2, "C", "B"))) %>% 
    mutate(Four_group = ifelse(CAPRA_S<2.5*-1, "A", ifelse(CAPRA_S<0, "B", ifelse(CAPRA_S<2.5, "C", "D")))) %>% 
    tibble::column_to_rownames(var = "Sample") %>% 
    as.data.frame()
  
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  
  
  #Two group comparasion
  
  
  { #Set DESeq data type
    dds_2 <- DESeqDataSetFromMatrix(countData=countData, 
                                    colData=coldata, 
                                    design=~Two_group)
    
    
    dds_2 <- DESeq(dds_2)
    
    
    res_2 <- results(dds_2, 
                     tidy = TRUE)
    
    #Extract unadjusted pvalue 
    Pvalue_2_group<- res_2 %>% 
      select(Pval_final = pvalue) %>% 
      as_tibble()
    
    }
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  
  
  # Three group comparasion
  
  
  {
    dds_3 <- DESeqDataSetFromMatrix(countData=countData, 
                                    colData=coldata, 
                                    design=~Three_group)
    
    dds_3 <- DESeq(dds_3)
    
    res_3_early<-results(dds_3, 
                         contrast=c("Three_group","A","B"),
                         tidy = TRUE)
    
    res_3_late<-results(dds_3, 
                        contrast=c("Three_group","C","B"),
                        tidy = TRUE)
    
    Pvalue_3_group<-cbind(
      res_3_early %>% 
        select(Early_3 = pvalue),
      res_3_late  %>% 
        select(Late_3 = pvalue)
    ) %>% 
      mutate(Pval_final = ifelse(
        Early_3<0.05&Late_3>=0.05, Early_3, 
        ifelse(
          Late_3<0.05&Early_3>=0.05, Late_3,
          sample(c(Late_3, Early_3, size =1)))) 
      ) %>% 
      as_tibble()
  }
  
  # Four group comparasion
  {
    dds_4 <- DESeqDataSetFromMatrix(countData=countData, 
                                    colData=coldata, 
                                    design=~Four_group)
    
    dds_4 <- DESeq(dds_4)
    
    
    res_4_AB_diff <- results(dds_4, 
                             contrast=c("Four_group", "A", "B"), 
                             tidy = TRUE) %>% 
      select(AB_4 = pvalue)
    
    
    res_4_BC_diff<- results(dds_4, 
                            contrast=c("Four_group", "B", "C"), 
                            tidy = TRUE) %>% 
      select(BC_4  = pvalue)
    
    
    res_4_CD_diff<- results(dds_4, 
                            contrast=c("Four_group", "C", "D"), 
                            tidy = TRUE) %>% 
      select(CD_4 = pvalue) 
    
    Pvalue_4_group<-cbind(res_4_AB_diff,
                          res_4_BC_diff,
                          res_4_CD_diff) %>% 
      rowwise() %>% 
      mutate(Pval_final = ifelse(AB_4<0.05&BC_4>=0.05&CD_4>=0.05, 
                                 AB_4,
                                 ifelse(AB_4>=0.05&BC_4>=0.05&CD_4<0.05, 
                                        CD_4,
                                        sample(c(AB_4, BC_4, CD_4), size = 1)
                                 )
      )
      ) %>% 
      as_tibble()
  }
  
  
  
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  
  
  return(list(
    Group_2_Pvalue = Pvalue_2_group,
    Group_3_Pvalue = Pvalue_3_group,
    Group_4_Pval = Pvalue_4_group))
  
}





##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 







#  Section 4 
#  Create a function which performs Bayseq analysis on simulated datasets


#BiocManager::install("baySeq")

{library(baySeq)
  library(dplyr)
  library(tidyr)
  library(tibble)}


# Currently needs number of simulated genes to be multiples of 1000 (to work on)

Bayseq_analysis<-function(simulated_data) {
  {require(baySeq)
    require(dplyr)
    require(tidyr)
    require(tibble)}
  
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
    distinct() %>% 
    max()
  
  
  #Calculate total number of simulated genes
  sim_genes_size<-simulated_data %>% 
    select(Gene_ref) %>% 
    dplyr::distinct() %>% 
    nrow()
  
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
  
  
  {library(doParallel)
    library(bigstatsr)
    
    
    cl <- parallel::makeCluster(bigstatsr::nb_cores())
    doParallel::registerDoParallel(cl)
    
    
    library(foreach)
    library(doSNOW)
    
    registerDoSNOW(cl)} 
  
  
  
  # Run analysis comparing Bayseq likelihood of 2 groups to null (of no DE)
  
  Bayseq_2_group=foreach(i=c(seq(from = 1, to = sim_genes_size, by =1000)),
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
                                                bootStraps = 3, #Recommonded bootstrap number of 3
                                                verbose = TRUE)
                           
                           #saveRDS(CD@posteriors, compress = "gzip", file = paste0("Bayseq", i, "2_group_2000_2.rds"))
                           
                           CD@posteriors
                           
                         }
 
  Bayseq_2_group<-lapply(Bayseq_2_group, function(x) {
    x %>% 
      as_tibble() } 
    ) %>% 
    dplyr::bind_rows() 
    
  
  Bayseq_3_group=foreach(i=c(seq(from = 1, to = sim_genes_size, by =1000)),
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
    
  Bayseq_3_group<-lapply(Bayseq_3_group, function(x) {
    x %>% 
      as_tibble() } 
  ) %>% 
    dplyr::bind_rows() 
  
  
  Bayseq_4_group=foreach(i=c(seq(from = 1, to = sim_genes_size, by =1000)),
                         .verbose = T) %do%  {
                           
                           CD <- new("countData", 
                                     data = simdata[i:(999+i),] %>% as.matrix(), 
                                     replicates = replicates %>% c(),
                                     groups = list(groups[[1]], groups[[6]], groups[[7]], groups[[8]], groups[[9]]))
                           
                           
                           
                           
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
  
  
  Bayseq_4_group<-lapply(Bayseq_4_group, function(x) {
    x %>% 
      as_tibble() } 
  ) %>% 
    dplyr::bind_rows() 
  
  return(list(group_2 = Bayseq_2_group, group_3 = Bayseq_3_group, group_4 = Bayseq_4_group))
   
}





#Section 4

# Create a function which takes pvalue results 
# From EdgeR and DESeq2 
# And returns table of for plotting FDR and ROC curves


fdr_roc<-function(pval_list, simulated_data, test_type = "EdgeR") {
  
  require(dplyr)
  require(magrittr)
  require(taRifx)
  require(tibble)
  require(qvalue)
  
  N_P<-simulated_data %>% 
    filter(Null_test == "FALSE") %>% 
    select(Gene_number) %>% 
    dplyr::distinct() %>% 
    nrow()
  
  N_N<-simulated_data %>% 
    filter(Null_test == "TRUE") %>% 
    select(Gene_number) %>% 
    dplyr::distinct() %>% 
    nrow()
  
  
  #   #    #    #    #     #    #    #     #    #    # 
  
  #FDR Set UP 
  
  require(foreach)
  
  qvalue_list=foreach(i=c(1,2,3)) %do% {
    
    qvalue(
      pval_list[[i]]  %$%
        Pval_final) %$%
      qvalues
    
  }
  
  
  qvalue_list<-lapply(
    pval_list,
    function(x){ 
      qvalue(
        x %$%
          Pval_final) %$%
        qvalues
      
    }
  )
  
  
  
  
  qvalue_list<-lapply(
    qvalue_list,
    function(x) {
      data.frame(qvalue = x) %>%
        as_tibble()
    }
  )
  
  
  
  FDR_list<-merge.list(
    pval_list,
    qvalue_list
  )
  
  
  #Number of significant hits at each qvalue (calculated by qvalue)  
  N_sig_pred<-lapply(FDR_list,
                     function(y) {
                       sapply(
                         seq(from=0, to=1, by=1/1000), #Qvalues / FDR values
                         function(x) {
                           sum(
                             (
                               y %$%
                                 qvalue)<x, 
                             na.rm = T)
                         }
                       )
                     }
  )
  
  #Convert above vectors into data frame form, to allow later merging
  
  N_sig_pred<-lapply(
    N_sig_pred,
    function(x) {
      data.frame(N_sig = x) %>% #Call column of Number of sig hits at each qval N_sig
        as_tibble() 
    }
  )
  
  #Calculate how many 
  
  
  if (test_type == "DESeq") {
    FDR_list<-lapply(FDR_list,
                     function(x) cbind(
                       x, simulated_data %>% 
                         select(Gene_number) %>% 
                         dplyr::distinct()
                     )) }
  
  
  
  Null_pval_list<-lapply(
    FDR_list,
    function(x) {
      inner_join(x,
                 simulated_data %>% 
                   filter(Null_test == "TRUE") %>% 
                   select(Gene_number) %>% 
                   dplyr::distinct()
                 
      )
    }
  )
  
  
  
  False_Positive<-lapply(Null_pval_list,
                         function(y) {
                           sapply(
                             seq(from=0, to=1, by=1/1000),
                             function(x) {
                               sum(
                                 (
                                   y %$%
                                     qvalue)<x, 
                                 na.rm = T)
                             }
                           )
                         }
  )
  
  
  False_Positive<-lapply(
    False_Positive,
    function(x) {
      data.frame(N_FP = x) %>% 
        as_tibble()
    }
  )
  
  FDR_table<-merge.list(
    N_sig_pred,
    False_Positive
  )
  
  FDR_table<-lapply(
    FDR_table,
    function(x) {
      x %>% 
        mutate(FDR_pred = seq(from=0, to=1, by=1/1000)) %>% 
        mutate(FDR_true = N_FP/N_sig) %>% 
        mutate_if(is.numeric, replace_na, 0) %>% 
        mutate(TP = N_sig-N_FP) %>% 
        mutate(TPR = TP/N_P) %>% 
        mutate(FPR  = N_FP/N_N)
    }
  )
  
  FDR_ROC_table<-FDR_table %>% 
    bind_rows() %>% 
    mutate(group = c(
      rep("Group_2", n()/3),
      rep("Group_3", n()/3), 
      rep("Group_4", n()/3)))
  
  return(FDR_ROC_table)
} 

#      #      #      #      #      #      #      #      #      # 





##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 







Fig_2_sim_data_sample_105 <- readRDS("~/TABI/dev/Article_Figures/Simulation_Data/Fig_2_sim_data_sample_105.rds")




#ROC Curve Calculations 

#TPR (True Positive Rate)


