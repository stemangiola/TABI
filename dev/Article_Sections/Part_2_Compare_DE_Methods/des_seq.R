
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Intention of Script 


# (Section 1)
# Create a function which runs DE analysis with 
# 2, 3, 4 groups (to compare "when" DE occurs)
# From a simulated dataset 


# (Section 2)
# Use this function on a simulated dataset 
# And to evalute/ create FDR / and ROC curves 



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 


# Section 1 
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






#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


# Section 2



test
DESeq2_pvalue<-DESeq_analysis(test %>% mutate(Sample = rep(1:(21*5), 15000)))



Bayseq %>% dim()
FPR<-sapply(seq(from=0, to=1, by=1/1000), 
       function(x) {
full_table<-cbind(Bayseq_small %>% 
        exp(), names = simdata %>% rownames()) %>% 
  as_tibble() %>% 
  filter(V2<x)

False<-full_table %>% filter(grepl("V",names)) %>% nrow()

Total<-full_table %>% nrow()

return(False/Total)

})

N_sig<-sapply(seq(from=0, to=1, by=1/1000), 
       function(x) {
         full_table<-cbind(Bayseq_small %>% 
                             exp(), names = simdata %>% rownames()) %>% 
           as_tibble() %>% 
           filter(V2<x)
         
         False<-full_table %>% filter(grepl("V",names)) %>% nrow()
         
         Total<-full_table %>% nrow()
         
         return(Total)
         
       })

rbind(TABI_fdr$ROC_FDR_table %>% select(FDR_true, N_sig) %>% mutate(group ="TABI"), 
data.frame(FDR_true = FPR, N_sig = N_sig) %>% 
  mutate(group = "Bayseq")) %>% 
  ggplot(aes(x=N_sig, y=FDR_true, group = group, col = group)) + 
  geom_path()



Bayseq %>% exp() 
library(ggplot2)


TABI_fdr<-TABI_fdr_roc_stats()


FDR_ROC_DEseq<-fdr_roc(DESeq2_pvalue, 
                       test %>% 
                         mutate(Sample = rep(1:(21*5), 15000)) %>% 
                        mutate(Gene_number = Gene_ref),
                 test_type = "DESeq"
  )


edgeR_pvalues<-edgeR_group_analysis(test %>% 
                                      mutate(sample = rep(1:(21*5), 15000)) %>% 
                                      mutate(Gene_number = Gene_ref))


FDR_ROC_edgeR<-fdr_roc(edgeR_pvalues,
                       test %>% 
                         mutate(Sample = rep(1:(21*5), 15000)) %>% 
                         mutate(Gene_number = Gene_ref)
                       )
FDR_ROC_k3$ROC_FDR_table %>% mutate(Test_type = TABI)



rbind(TABI_fdr$ROC_FDR_table %>% 
        mutate(Test_type = "TABI") %>% 
        mutate(group  = "NA") %>% 
        mutate(N_FP = "NA") %>% 
        mutate(TP = "NA"),
  FDR_ROC_edgeR %>% 
        mutate(Test_type = paste("edgeR", group)),
FDR_ROC_DEseq %>% 
  mutate(Test_type = paste("DESeq2", group))) %>% 
  ggplot(aes(x=N_sig, y= FDR_pred, group = Test_type, col = Test_type)) + 
  geom_line(linetype = 2, alpha = 0.5) + 
  geom_line(aes(y=FDR_true)) + 
  labs(y= "FDR",
       x = "Number of Tests (Simulated Genes) Declared Significant",
       title = "Figure 2. FDR Curve Power Comparision") + 
  ylim(0, 0.1) + 
  xlim(0, 1000) + 
  scale_color_manual(values=c("forestgreen"," springgreen4", "springgreen3", "steelblue1", "steelblue2", "steelblue3", "red")) + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
          panel.grid.minor= element_line(colour='white', linetype = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))

library(ggplot2)

coord_cartesian(xlim =c(0, 0.1), ylim = c(0, 1000))















