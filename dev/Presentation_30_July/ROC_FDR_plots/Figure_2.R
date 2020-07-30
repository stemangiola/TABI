
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #





#Intention of Script 

#Use functions defined in script Functions_for_FDR_ROC_curve to run
#edgeR, DESeq2, Bayseq analysis on defined dataset 
#And extract FDR results to produce Figure 2



#(Section 1)
#Run analysis / obtain FDR results 




#(Section 2)
#Produce Figure 2 (FDR Curves)



#(Section 3)
#Produce ROC Curves





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #





##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 




# Section 1 

#Run analysis of simulated data and obtain FDR results



#Load Simulated Data

Fig_2_sim_data_sample_105 <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Presentation_30_July/ROC_FDR_plots/Fig_2_sim_data_sample_105.rds")




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

#DESeq2

#Run DESeq2 analysis and obtain raw pvalues

DESeq2_pvalue<-DESeq_analysis(Fig_2_sim_data_sample_105 %>% 
                                dplyr::rename(Sample = sample))

#Use raw pvalues from DESeq2 obtain FDR / ROC curve values

FDR_ROC_DEseq<-fdr_roc(DESeq2_pvalue, 
                         Fig_2_sim_data_sample_105 %>% 
                         dplyr::rename(Sample = sample)%>% 
                         mutate(Gene_number = Gene_ref),
                       test_type = "DESeq"
)




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

#edgeR


# Run edgeR analysis and use this to obtain raw pvalues

edgeR_pvalues<-edgeR_group_analysis(Fig_2_sim_data_sample_105 %>% 
                                      mutate(Gene_number = Gene_ref))


#Use raw pvalues from edgeR analtsis to obtain FDR / ROC curve values

FDR_ROC_edgeR<-fdr_roc(edgeR_pvalues,
                       Fig_2_sim_data_sample_105 %>% 
                         dplyr::rename(Sample = sample) %>% 
                         mutate(Gene_number = Gene_ref)
)





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

#Bayseq

# Run Bayseq analysis and use this obtain log likelihoods / posteriors

Posterior_Bayseq<-Bayseq_analysis(Fig_2_sim_data_sample_105 %>% 
                                    dplyr::rename(Sample = sample))



#Use log likelihoods to obtain FDR / ROC curve values  

FDR_ROC_Bayseq_df<-lapply(Posterior_Bayseq, function(x) {
  
  FDR_ROC_Bayseq(x, Fig_2_sim_data_sample_105 %>% 
               dplyr::rename(Sample = sample))
  
}) 

  
  #Load bayseq results to save time 


#FDR_ROC_Bayseq_df <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Presentation_30_July/ROC_FDR_plots/FDR_ROC_Bayseq_df.rds")



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


#TABI 

#TABI has been prerun (see  Simulation_Data Folder, script : Sig_sigmulation_trad_parametisation, Section 2 )


#Load TABI posterior / fit results 

#TABI_2000_iter_Fig_2_sim_data <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Presentation_30_July/ROC_FDR_plots/TABI_2000_iter_Fig_2_sim_data.rds")


# Use TABI Posterior Results to calculate FDR curve values for TABI

# FDR_TABI<-TABI_fdr_roc_stats(TABI_2000_iter_Fig_2_sim_data[[2]] %>% 
#                                dplyr::select(-gene) %>% 
#                                dplyr::rename(gene = Gene_ref),                        
#                              Fig_2_sim_data_sample_105 %>% 
#                                dplyr::rename(Sample = sample) %>% 
#                                mutate(Gene_number = Gene_ref))


# Load TABI FDR Results (TABI posterior results too large for github)

FDR_TABL <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Presentation_30_July/ROC_FDR_plots/FDR_TABL.rds")





##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 






# Section 2
# Combine FDR results into single dataframe
# To plot FDR curves for the different methods


{ library(ggplot2)
  library(ggpubr)
}




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Combine FDR curve results into single dataframe for plotting Figure 2



FDR_table<-rbind(FDR_TABI$ROC_FDR_table %>% 
                   mutate(Test_type = "TABI") %>% 
                   select(FDR_true, N_sig, Test_type),
                 FDR_ROC_edgeR %>% 
                   mutate(Test_type = paste("edgeR", group)) %>% 
                   select(FDR_true, N_sig, Test_type),
                 FDR_ROC_DEseq %>% 
                   mutate(Test_type = paste("DESeq2", group)) %>% 
                   select(FDR_true, N_sig, Test_type),
                 FDR_ROC_Bayseq_df[[1]] %>% 
                   mutate(FDR_true = FP/N_sig) %>% 
                   select(N_sig, FDR_true) %>% 
                   mutate(Test_type = "Bayseq_2_group"),
                 FDR_ROC_Bayseq_df[[2]]  %>%
                   mutate(FDR_true = FP/N_sig) %>% 
                   select(N_sig, FDR_true) %>% 
                   mutate(Test_type = "Bayseq_3_group"),
                 FDR_ROC_Bayseq_df[[3]] %>% 
                   mutate(FDR_true = FP/N_sig) %>% 
                   select(N_sig, FDR_true) %>% 
                   mutate(Test_type = "Bayseq_4_group")
)





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



# Zoomed in FDR Plot Comparasion

# Realistic Range of FDR 


FDR_table %>% 
  ggplot(aes(x=N_sig, 
             y= FDR_true, 
             group = Test_type, 
             col = Test_type))  +
  geom_line() + 
  labs(y= "False Discovery Rate",
       x = "Number of Tests (Simulated Genes) Declared Significant",
       title = "FDR Curve Power Comparison"
  ) + 
  coord_cartesian(xlim =c(0, 1000), ylim = c(0, 0.1)) + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))


#ggsave("FDR_comparison_Plot_Realistic_FDR_range.pdf")


zoom_in<-FDR_table %>% 
  ggplot(aes(x=N_sig, 
             y= FDR_true, 
             group = Test_type, 
             col = Test_type))  +
  geom_line() + 
  labs(y= "",
       x = ""
       #title = "Figure 2. FDR Curve Power Comparision"
  ) + 
  coord_cartesian(xlim =c(0, 1000), ylim = c(0, 0.1)) + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))






#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



# Full range of FDR 

FDR_table %>% 
  ggplot(aes(x=N_sig,
             y = FDR_true, 
             group = Test_type, 
             col = Test_type))  +
  geom_line() + 
  labs(y= "False Discovery Rate",
       x = "Number of Tests (Simulated Genes) Declared Significant",
       title = "FDR Curve Power Comparison"
  ) + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))


#ggsave("FDR_comparison_full_range.pdf")


full<- FDR_table %>% 
  ggplot(aes(x=N_sig,
             y = FDR_true, 
             group = Test_type, 
             col = Test_type))  +
  geom_line() + 
  labs(y= "",
       x = ""
       #title = "Figure 2. FDR Curve Power Comparision"
  ) + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))




annotate_figure(
annotate_figure(ggarrange(full, zoom_in,
                          ncol = 2,  
                          label.y = "FDR",
                          label.x = "Number of Tests (Simulated Genes) Declared Significant",
                          common.legend = TRUE, 
                          legend = "bottom"), 
                top = text_grob("Figure 2. FDR Curve Power Comparision", hjust = 1.5),
                bottom = text_grob("Number of Tests (Simulated Genes) Declared Significant"),
                left  = text_grob("FDR", rot = 90)),
                bottom = text_grob("TABI outperforms DESeq2 and edgeR in False Discovery Rate for realistic False Discoveries and simulated data (150,000 Simulated Genes, 15,000 Truly Differentially Expressed). 
                                   Non-TABI methods were performed under recommended settings (DESeq2: sample size = 10000,bootstrap = 3) ", size = 8))



#ggsave("Figure_2_FDR_curve.pdf")




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ## ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



# Section 3 

# ROC Curves 


# Combine ROC results into single dataframe
# To plot ROC curves for the different methods


{ library(ggplot2)
  library(ggpubr)
}



#Combine ROC curve results into single dataframe for plotting Figure 2


ROC_table<-rbind(FDR_TABI$ROC_FDR_table %>% 
                   mutate(Test_type = "TABI") %>% 
                   select(TPR, FPR, Test_type),
                 FDR_ROC_edgeR %>% 
                   mutate(Test_type = paste("edgeR", group)) %>% 
                   select(TPR, FPR, Test_type),
                 FDR_ROC_DEseq %>% 
                   mutate(Test_type = paste("DESeq2", group)) %>% 
                   select(TPR, FPR, Test_type),
                 FDR_ROC_Bayseq_df[[1]] %>% 
                   select(TPR, FPR) %>% 
                   mutate(Test_type = "Bayseq_2_group"),
                 FDR_ROC_Bayseq_df[[2]]  %>%
                   select(TPR, FPR) %>% 
                   mutate(Test_type = "Bayseq_3_group"),
                 FDR_ROC_Bayseq_df[[3]] %>% 
                   select(TPR, FPR) %>% 
                   mutate(Test_type = "Bayseq_4_group")
)


# ROC Curve Plot 


ROC_table %>% 
  ggplot(aes(x=FPR, 
             y=TPR, 
             group = Test_type, 
             col = Test_type)) + 
  geom_line() + 
  labs(x = "1 − Specificity",
       y= "Sensitivity", 
       title = "Receiver Operating Characteristic (ROC) Curve") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))


#ggsave("ROC_curve.pdf")


