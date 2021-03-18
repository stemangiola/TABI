
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












##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



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

#ROC Curve Calculations 

#TPR (True Positive Rate)


heck %>% 
  ggplot(aes(x=FPR, y=TPR, col = group)) 

heck %>% 
  ggplot(aes(x=FDR_pred, y= N_sig, col=group)) +
  geom_line() +
  geom_line(aes(x=FDR_true, y=N_sig), linetype = 2)

#Create a list of data frames where each data frame is the same dataset, but with a different sized sample

#Add sample identifier column
Sim_test<-sigmodial_test %>% 
  mutate(sample = rep(1:150, 1000))


Sim_test_sample<-lapply(
  seq(from = 10, to = 140, by = 10),
  function(x) {
    #Ensuring same result / sample between reproductions
    set.seed(15)
    #Select which samples to keep
    filter_by<-sample.int(150, size = x)
    #Flter out only those samples in the data
    inner_join(Sim_test, 
               data.frame(sample = filter_by))
  }
)

#Create names for the above list elements which correspond to how many samples are in each dataset

Sample_Names<-sapply( 
  seq(from = 10, to = 140, by = 10),
  function(x) paste0("Sim_test_", x, "sample"))

#Set names as above 
names(Sim_test_sample)<-Sample_Names

EdgeR_Analysis<-lapply(
  Sim_test_sample,
  function(x) simuledgeR_fdr_roc(
    edgeR_group_analysis(x), 
    x)
)

n_len<-EdgeR_Analysis %>% 
  bind_rows() %>% 
  nrow()

EdgeR_Analysis %>% 
  bind_rows() %>% mutate(Sample_size = 
c(sapply(seq(from = 10, to = 140, by = 10), 
        function(x) rep(x, n_len/14)))) %>% 
  group_by(group, Sample_size) %>% 
  ggplot(aes(x=FPR, y=TPR, col= as.factor(Sample_size))) +
  facet_grid(. ~ group) +
  geom_line() +
  geom_abline(slope =1, intercept =0, linetype = 2) + 
  labs(title = "ROC Curve EdgeR Analysis", 
       subtitle = "Dashed line is line of Randomness", 
       x = "FPR (False Postive Rate (1-Specificity)", 
       y= "TPR (True Positive Rate, Sensitivity)")

EdgeR_Analysis %>% 
  bind_rows() %>% mutate(Sample_size = 
                           c(sapply(seq(from = 10, to = 140, by = 10), 
                                    function(x) rep(x, n_len/14)))) %>% 
  group_by(group, Sample_size) %>% 
  ggplot(aes(x=FDR_pred, y= N_sig, col = as.factor(Sample_size))) + 
  facet_grid(. ~ group) +
  geom_line() +
  geom_line(aes(x=FDR_true, y=N_sig, group = Sample_size), linetype = 2) + 
  labs(title = "FDR Curve EdgeR Analysis (Varing Sample Size)", 
       subtitle = "Dashed Line Indicates True / Actual Curve, Filled line Indicates Expected Curve",
       y="Number of genes declared significant", x= "FDR")



EdgeR_Analysis %>% 
  bind_rows() %>% mutate(Sample_size = 
                           c(sapply(seq(from = 10, to = 140, by = 10), 
                                    function(x) rep(x, n_len/14)))) %>% 
  filter(Sample_size<100&Sample_size>20) %>% 
  ggplot(aes(x=FPR, y=TPR, col= as.character(Sample_size))) +
  facet_grid(. ~ group) +
  geom_line() 

EdgeR_Analysis %>% 
  bind_rows() %>% mutate(Sample_size = 
                           c(sapply(seq(from = 10, to = 140, by = 10), 
                                    function(x) rep(x, n_len/14)))) %>% 
  group_by(group, Sample_size) %>% 
  filter(Sample_size<100&Sample_size>20) %>% 
  ggplot(aes(x=FDR_pred, y= N_sig, col = as.factor(Sample_size))) + 
  facet_grid(. ~ group) +
  geom_line() +
  geom_line(aes(x=FDR_true, y=N_sig, group = Sample_size), linetype = 2)


sample_posterior_list<-list(Sim_test_30sample_2000_iter,
                            Sim_test_40sample_2000_iter,
                            Sim_test_50sample_2000_iter,
                            Sim_test_60sample_2000_iter,
                            Sim_test_70sample_2000_iter,
                            Sim_test_80sample_2000_iter,
                            Sim_test_90sample_2000_iter,
                            Sim_test_100sample_2000_iter,
                            Sim_test_110sample_2000_iter,
                            Sim_test_120sample_2000_iter)

TABI_FDR_ROC<-lapply(
  sample_posterior_list,
  function(x) TABI_fdr_roc_stats(x, sigmodial_test)
)



TABI_FDR_ROC_Res<-lapply(
  TABI_FDR_ROC,
  function(x) x %$%
    ROC_FDR_table
)

TABI_FDR_ROC_Res %>% 
  bind_rows() 

n_len_TABI<-TABI_FDR_ROC_Res%>% 
  bind_rows() %>% 
  nrow()

list(EdgeR_Analysis %>% 
  bind_rows() %>% mutate(Sample_size = 
                           c(sapply(seq(from = 10, to = 140, by = 10), 
                                    function(x) rep(x, n_len/14)))) %>% 
  filter(Sample_size<130&Sample_size>20) %>% 
  mutate(Test_type = paste0("edgeR_", group)),


TABI_FDR_ROC_Res %>% 
  bind_rows() %>% mutate(Sample_size = 
                           c(sapply(seq(from = 30, to = 120, by = 10), 
                                    function(x) rep(x, n_len_TABI/10)))) %>% 
  mutate(Test_type = "TABI_2000_Iter")) %>% 
  bind_rows() %>% 
  mutate(slope = 1, intercept =0) %>% 
  ggplot(aes(x=FPR, y=TPR, col= Test_type)) +
  facet_wrap(. ~ as.factor(Sample_size)) +
  geom_line(aes(col = Test_type)) + 
  scale_color_manual(values=c("green3", "dodgerblue", "black", "red"))  +
  geom_abline(aes(slope = slope, intercept = intercept), linetype = 2) 

list(EdgeR_Analysis %>% 
       bind_rows() %>% mutate(Sample_size = 
                                c(sapply(seq(from = 10, to = 140, by = 10), 
                                         function(x) rep(x, n_len/14)))) %>% 
       filter(Sample_size<130&Sample_size>20) %>% 
       mutate(Test_type = paste0("edgeR_", group)),
     
     
     TABI_FDR_ROC_Res %>% 
       bind_rows() %>% mutate(Sample_size = 
                                c(sapply(seq(from = 30, to = 120, by = 10), 
                                         function(x) rep(x, n_len_TABI/10)))) %>% 
       mutate(Test_type = "TABI_2000_Iter")) %>% 
  bind_rows() %>% 
  ggplot(aes(x=FDR_true, y=N_sig, col= Test_type)) +
  facet_wrap(. ~ as.factor(Sample_size)) +
  geom_line()  +
  geom_line(aes(x=FDR_pred), linetype = 2) +
  scale_color_manual(values=c("green3", "dodgerblue", "black", "red")) 





TABI<-TABI_FDR_ROC_Res %>% 
  bind_rows() %>% mutate(Sample_size = 
                           c(sapply(seq(from = 30, to = 90, by = 10), 
                                    function(x) rep(x, n_len_TABI/7)))) %>% 
  dplyr::rename(FDR_true = True_FDR_TABI, FDR_pred = FDR) %>% 
  mutate(Test_type = "TABI")

Edge_R<-EdgeR_Analysis %>% 
  bind_rows() %>% mutate(Sample_size = 
                           c(sapply(seq(from = 10, to = 140, by = 10), 
                                    function(x) rep(x, n_len/14)))) %>% 
  filter(Sample_size<100&Sample_size>20) %>% 
  mutate(Test_type = paste0("edgeR_", group)) %>% 
  select(N_sig, FDR_pred, FDR_true, FPR, TPR, Test_type, Sample_size) 

bind_rows(TABI, Edge_R) %>% 
  #filter(Sample_size==30|Sample_size==60|Sample_size==90) %>% 
  ggplot(aes(col = Test_type)) +
  facet_wrap(. ~ Sample_size) + 
  geom_line(aes(x=FDR_pred, y= N_sig)) +
  geom_line(aes(x=FDR_true, y= N_sig), linetype = 2) + 
  scale_color_manual(values=c("green3", "dodgerblue", "black", "red"))



bind_rows(TABI, Edge_R) %>% 
  ggplot(aes(x=FPR, y=TPR, col= Test_type)) +
  facet_wrap(. ~ Sample_size) + 
  geom_line() + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_color_manual(values=c("green3", "dodgerblue", "black", "red"))

