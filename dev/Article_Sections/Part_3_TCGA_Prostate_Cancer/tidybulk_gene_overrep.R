




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #





# Intention of Script 

# Show Bargraph / quadrant plot of pathway analysis for TCGA Analysis







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

  