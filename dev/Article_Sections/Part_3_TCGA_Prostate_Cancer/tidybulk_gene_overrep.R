




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #





# Intention of Script 

# Show Bargraph / quadrant plot of pathway analysis for TCGA Analysis


non_zero_slope_genes<- fit_table_purity_complete %>% 
  filter(parameters == "beta[2,1]") %>% 
  filter(`2.5%`<0&`97.5%`<0|`2.5%`>0&`97.5%`>0) %>% #Selecting either all > 0 for 95% interval / or all < 0
  select(Gene_name) #Which genes


fit_table_purity_complete_sig <-inner_join(non_zero_slope_genes, fit_table_purity_complete)

inflection_plot<-fit_table_purity_complete %>% #Data Table of sig fit data
  mutate(Slope_significance = ifelse(`2.5%`<0&`97.5%`<0|`2.5%`>0&`97.5%`>0, "Significant (n=3695)", "Non-Significant (n=31167)")) %>% 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% #Select required parameters
  select(parameters, mean, Gene_name, Slope_significance) %>%
  filter(!Gene_name == "TUBBP10") %>% #Two different val for TUBP10 / possible dup?
  pivot_wider(names_from = parameters, values_from = mean) %>% #Reshape table for plotting
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]")%>% 
  ggplot(aes(x=Inflection_mean, y= Slope_mean, color=Slope_significance)) +
  geom_point(alpha=0.7) + 
  scale_color_manual(name = "", values = c("Significant (n=3695)" = "dodgerblue","Non-Significant (n=31167)"= "black")) +
  geom_point(col="black", alpha=0.7) + 
  geom_point(aes(x=Inflection_mean, y= Slope_mean), alpha=0.7, col="dodgerblue", data = fit_table_purity_complete_sig %>% #Data Table of sig fit data
               filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% #Select required parameters
               select(parameters, mean, Gene_name) %>%
               filter(!Gene_name == "TUBBP10") %>% #Two different val for TUBP10 / possible dup?
               pivot_wider(names_from = parameters, values_from = mean) %>% #Reshape table for plotting
               rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]")) +
  stat_smooth(se = FALSE, method = 'lm', col="red") + #Plot fitted linear eq
  labs(title = "Fig 2. TABI Calculated Slope mean vs Inflection mean",
       y="Calculated Slope Mean", 
       x="Calculated Inflection Mean",
       caption ="Each point represents a single gene analysed by TABI (n = 34862)") +
  custom_theme





custom_theme <-
  list(
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





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



{

  
  
  
small_inflection_early<-fit_table_purity_complete %>% 
  filter(parameters == "inflection[1]") %>% 
  filter(sd<0.7375) %>% 
  filter(mean < 0) %>% 
  select(Gene_name)




small_inflection_later<-fit_table_purity_complete %>% 
  filter(parameters == "inflection[1]") %>% 
  filter(sd<0.7375) %>% 
  filter(mean > 0) %>% 
  select(Gene_name)




increase_genes<- fit_table_purity_complete %>% 
  filter(parameters == "beta[2,1]") %>% 
  filter(`2.5%`>0&`97.5%`>0) %>% #Selecting either all > 0 for 95% interval / or all < 0
  select(Gene_name)



decrease_genes<- fit_table_purity_complete %>% 
  filter(parameters == "beta[2,1]") %>% 
  filter(`2.5%`<0&`97.5%`<0) %>% #Selecting either all > 0 for 95% interval / or all < 0
  select(Gene_name)




IE_names<-inner_join(small_inflection_early, 
               increase_genes) %>% 
  rename(transcript = Gene_name ) %>% 
  mutate(Type = "IE")



IL_names<-inner_join(small_inflection_later, 
               increase_genes) %>% 
  rename(transcript = Gene_name ) %>% 
  mutate(Type = "IL")



DE_names<-inner_join(small_inflection_early, 
               decrease_genes) %>% 
  rename(transcript = Gene_name ) %>% 
  mutate(Type = "DE")




DL_names<-inner_join(small_inflection_later, 
               decrease_genes) %>% 
  rename(transcript = Gene_name) %>% 
  mutate(Type = "DL")




}



TCGA_simples<-ed_PC_TCGA %>% 
  select(Sample.ID, ens_iso, transcript, `read count normalised`) %>% 
  separate(ens_iso, c("ens", "iso"), "[.]") 


DL_table<-full_join(TCGA_simples, DL_names) %>% 
  mutate(do_test = ifelse(is.na(Type), FALSE, TRUE))

tt_DL<-tidybulk::tidybulk(
  .data = DL_table,
  .sample = Sample.ID,
  .transcript = transcript,
  .abundance = `read count normalised`
)


tidybulk::test_gene_overrepresentation(
  tt_DL,
  .entrez =  transcript,
  .do_test = do_test,
  species = "Homo sapiens") 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

tidybulk::test_gene_overrepresentation(
  .data,
  .sample = NULL,
  .entrez,
  .do_test,
  species,
  gene_set = NULL
)
tt<-tidybulk::tidybulk(.data = ed_PC_TCGA, .sample = Sample.ID, .transcript = transcript, .abundance = `read count normalised`)
  
ed_PC_TCGA %>% names()

#%>% 
  mutate(do_test = ifelse(transcript == "A1BG", TRUE, FALSE))


tt_TCGA<-tidybulk::tidybulk(TCGA_simples, .sample = Sample.ID, .transcript = transcript, .abundance = `read count normalised`)

tidybulk::test_gene_overrepresentation(
  tt_TCGA,
  .entrez =  transcript,
  .do_test = do_test,
  species = "Homo sapiens") 








non_zero_slope_genes<- fit_table_purity_complete %>% 
  filter(parameters == "beta[2,1]") %>% 
  filter(`2.5%`<0&`97.5%`<0|`2.5%`>0&`97.5%`>0) %>% #Selecting either all > 0 for 95% interval / or all < 0
  select(Gene_name)




fit_table_purity_complete_sig <-inner_join(non_zero_slope_genes, 
                                           fit_table_purity_complete)


Increase_sig<-fit_table_purity_complete_sig %>% 
  filter

  