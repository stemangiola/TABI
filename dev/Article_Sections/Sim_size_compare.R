



Fig_2_sim_data_sample_105 <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Presentation_30_July/ROC_FDR_plots/Fig_2_sim_data_sample_105.rds")


{

set.seed(40)
Fig_2_sim_data_sample_21<-inner_join(
  Fig_2_sim_data_sample_105,
  Fig_2_sim_data_sample_105 %>% 
    dplyr::select(CAPRA_S, sample) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(CAPRA_S) %>% 
    sample_n(1))  

set.seed(40)
Fig_2_sim_data_sample_42<-inner_join(
  Fig_2_sim_data_sample_105,
  Fig_2_sim_data_sample_105 %>% 
    dplyr::select(CAPRA_S, sample) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(CAPRA_S) %>% 
    sample_n(2))

set.seed(40)
Fig_2_sim_data_sample_63<-inner_join(
  Fig_2_sim_data_sample_105,
  Fig_2_sim_data_sample_105 %>% 
    dplyr::select(CAPRA_S, sample) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(CAPRA_S) %>% 
    sample_n(3))

set.seed(40)
Fig_2_sim_data_sample_84<-inner_join(
  Fig_2_sim_data_sample_105,
  Fig_2_sim_data_sample_105 %>% 
    dplyr::select(CAPRA_S, sample) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(CAPRA_S) %>% 
    sample_n(4))
}



ROC_FDR_non_TABI<-function(sim_data){
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #DESeq2
  
  #Run DESeq2 analysis and obtain raw pvalues
  
  DESeq2_pvalue<-DESeq_analysis(sim_data %>% 
                                  dplyr::rename(Sample = sample))
  
  #Use raw pvalues from DESeq2 obtain FDR / ROC curve values
  
  FDR_ROC_DEseq<-fdr_roc(DESeq2_pvalue, 
                         sim_data %>% 
                           dplyr::rename(Sample = sample)%>% 
                           mutate(Gene_number = Gene_ref),
                         test_type = "DESeq"
  )
  
  
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #edgeR
  
  
  # Run edgeR analysis and use this to obtain raw pvalues
  
  edgeR_pvalues<-edgeR_group_analysis(sim_data %>% 
                                        mutate(Gene_number = Gene_ref))
  
  
  #Use raw pvalues from edgeR analtsis to obtain FDR / ROC curve values
  
  FDR_ROC_edgeR<-fdr_roc(edgeR_pvalues,
                         sim_data %>% 
                           dplyr::rename(Sample = sample) %>% 
                           mutate(Gene_number = Gene_ref)
  )
  
  
  
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #Bayseq
  
  # Run Bayseq analysis and use this obtain log likelihoods / posteriors
  

  


  
  Posterior_Bayseq<-Bayseq_analysis(sim_data %>% 
                                      dplyr::rename(Sample = sample))
  
  #Use log likelihoods to obtain FDR / ROC curve values  
  
  FDR_ROC_Bayseq_df<-lapply(Posterior_Bayseq, function(x) {
    
    FDR_ROC_Bayseq(x, sim_data %>% 
                     dplyr::rename(Sample = sample))
    
  }) 
  
  return(list(FDR_ROC_Bayseq_df  = FDR_ROC_Bayseq_df,
    FDR_ROC_edgeR = FDR_ROC_edgeR,
    FDR_ROC_DEseq = FDR_ROC_DEseq))

}   

{
  slope_1 <- readRDS("~/slope_1.rds")
  stats_non_TABI_slope_1<-non_TABI_stats(slope_1)
  saveRDS(stats_non_TABI_slope_1, "stats_non_TABI_slope_1.rds")
}

{
  slope_2 <- readRDS("~/slope_2.rds")
  stats_non_TABI_slope_2<-non_TABI_stats(slope_2)
  saveRDS(stats_non_TABI_slope_2, "stats_non_TABI_slope_2.rds")
}


{
  FDR_sim_105<-ROC_FDR_non_TABI(Fig_2_sim_data_sample_105)
saveRDS(FDR_sim_105, compress= "gzip", file = "FDR_sim_105.rds")

FDR_sim_84<-ROC_FDR_non_TABI(Fig_2_sim_data_sample_84)
saveRDS(FDR_sim_84, compress= "gzip", file = "FDR_sim_84.rds")


FDR_sim_63<-ROC_FDR_non_TABI(Fig_2_sim_data_sample_63)
saveRDS(FDR_sim_63, compress= "gzip", file = "FDR_sim_63.rds")


FDR_sim_42<-ROC_FDR_non_TABI(Fig_2_sim_data_sample_42)
saveRDS(FDR_sim_42, compress= "gzip", file = "FDR_sim_42.rds")


FDR_sim_21<-ROC_FDR_non_TABI(Fig_2_sim_data_sample_21)
saveRDS(FDR_sim_21, compress= "gzip", file = "FDR_sim_21.rds")


} 



non_TABI_stats<-function(sim_data){
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #DESeq2
  
  #Run DESeq2 analysis and obtain raw pvalues
  
  DESeq2_pvalue<-DESeq_analysis(sim_data)
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #edgeR
  
  
  # Run edgeR analysis and use this to obtain raw pvalues
  
  edgeR_pvalues<-edgeR_group_analysis(sim_data %>% 
                                        mutate(Gene_number = Gene_ref) %>% 
                                        mutate(sample = Sample))
  
  
  
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  #Bayseq
  
  # Run Bayseq analysis and use this obtain log likelihoods / posteriors
  
  
  
  
  Posterior_Bayseq<-Bayseq_analysis(sim_data)
  

  
  return(list(
    DESEq2 = DESeq2_pvalue,
              edgeR = edgeR_pvalues,
              Bayseq = Posterior_Bayseq))
}   

