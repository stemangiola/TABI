


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



# Intention of script

# Plotting of 6 examples of sigmoidal trends in TCGA dataset





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #





#Create Plotting Function which takes read count table, gene name and fit summary table 

#And plots graph of read count normalised vs scaled CAPRA-S


plot_ed_with_line<-function(gene_name, read_count_data, fit_table) {
  
  #Required packages / functions 
  require(ggplot2)
  require(tibble)
  require(dplyr)
  require(magrittr)
  require(rstan)
  require(tidyr)
  
  #Requires scale_design and parse_formula from TABI to create scaled design matrix 
  
  #Select only paramters for specified gene
  Gene_table <- fit_table %>% 
    filter(Gene_name == gene_name)
  
  plot_line <- function(x) {
    
    eta <- as.numeric(Gene_table %>%
                        filter(parameters == "inflection[1]") %>% 
                        select(mean_value))
    
    
    beta_1 <- as.numeric(Gene_table %>%
                           filter(parameters == "beta1_z[1,1]") %>%
                           select(mean_value))
    
    # beta_2<-as.numeric(Gene_table %>%
    #                      filter(parameters == "beta1_z[2,1]") %>%
    #                      select(mean_value))
    
    #beta_3
    
    y_0 <- as.numeric(Gene_table %>%
                        filter(parameters == "y_cross[1]") %>%
                        select(mean_value))
    
    A<- as.numeric(Gene_table %>%
                     filter(parameters == "A[1]") %>%
                     select(mean_value))
    
    lin_eq<-(-x*beta_1+eta*beta_1)
    
    top<-(y_0-A)*(1+exp(eta*beta_1))
    
    bottom<-(1+exp(lin_eq))
    
    return((exp(A+(top/bottom))+1))}
  
  parse_formula <- function(fm) {
    if(attr(terms(fm), "response") == 1) stop("The formula must be of the kind \"~ covariates\" ")
    else as.character(attr(terms(fm), "variables"))[-1]
  }
  
  scale_design = function(df, formula){
    df %>%
      setNames(c("sample_idx", "(Intercept)", 
                 parse_formula(formula))) %>%
      gather(cov, value, -sample_idx) %>%
      group_by(cov) %>%
      mutate( value = ifelse( !grepl("Intercept", cov) & length(union(c(0,1), value)) != 2, scale(value), value )) %>%
      ungroup() %>%
      spread(cov, value) %>%
      arrange(as.integer(sample_idx)) %>%
      select(`(Intercept)`, one_of(parse_formula(formula)))
  }
  
  
  #Create table with scaled CAPRA_S and normalised gene count for graphing
  graphing_table <- data.frame( #Scale CAPRA_S using TABI method 
    CAPRA_S = as.vector(model.matrix(object = ~CAPRA_S, data = read_count_data %>% 
                                       filter(is.na(CAPRA_S)==F, is.na(transcript)==F) %>% 
                                       filter(transcript ==gene_name)) %>%
                          as_tibble(rownames="sample_idx") %>%
                          scale_design(~CAPRA_S) %>% select(CAPRA_S)),
    #Combine with read counts 
    read_count_normalised = as.vector(read_count_data %>% filter(is.na(CAPRA_S)==F) %>% 
                                        filter(transcript ==gene_name) %>% 
                                        select(`read count normalised`)) )
  
  line_data<-sapply(seq(from = -2, to= 2, by = 0.1), 
                    function(x) plot_line(x))
  
  line_df<-data.frame(y= line_data, 
                      CAPRA_S = seq(from = -2, to= 2, by = 0.1))
  
  #Plot curve of specified parameters and of count points
  graphing_table %>% 
    ggplot(aes(x=CAPRA_S, y=read.count.normalised)) + 
    geom_point(col="dodgerblue") +
    stat_function(fun=plot_line, geom="line") +
    geom_path(aes(y=y), data = line_df) + 
    scale_y_log10() +
    labs(title=paste(gene_name),
         x = "",
         y = "") #+ 
  # theme(axis.text.y = element_blank(),
  #      axis.text.x = element_blank())
  
  

}





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



# Load TCGA read count table with reduced number of transcripts  (only those plotted below)

df <- readRDS("~/TABI/dev/Presentation_30_July/Regression_Eg/example_TCGA.rds")

# Load fit table 

fit_table_purity_complete <- readRDS("~/TABI/dev/Presentation_30_July/Regression_Eg/fit_table_purity_complete.rds")





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



# Plotting individual genes with their fitted equations





{ # Increasing (Top Row)
  
  #Increasing Early 
  
  #Alternative Increasing Early - SPRY4 more dramatic trend, but not as early as C12orf49
IE<-plot_ed_with_line("C12orf49", 
                      df %>% 
                        dplyr::rename(CAPRA_S = `CAPRA-S`), 
                      fit_table_purity_complete %>% 
                        rename(mean_value = mean))  + 
    labs(subtitle="Early Increase") + 
    theme(panel.grid.major = element_line(colour='white', linetype = 1),
          panel.grid.minor= element_line(colour='white', linetype = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))

IE<-plot_ed_with_line("SPRY4",
                      df%>%
                        dplyr::rename(CAPRA_S = `CAPRA-S`),
                      fit_table_purity_complete %>%
                        rename(mean_value = mean)) +
  labs(subtitle="Early Decrease") +
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

 #Middle Increase
IM<-plot_ed_with_line("IQGAP3", 
                      df%>% 
                        dplyr::rename(CAPRA_S = `CAPRA-S`), 
                      fit_table_purity_complete %>% 
                        rename(mean_value = mean))  + 
  labs(subtitle="Middle Increase") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

# Late Increase

IL<-plot_ed_with_line("CXCR5", 
                      df%>% 
                        dplyr::rename(CAPRA_S = `CAPRA-S`), 
                      fit_table_purity_complete %>% 
                        rename(mean_value = mean))  + 
  labs(subtitle="Late Increase") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

} 

{# Decreasing
  
  #Decreasing Early
DE<-plot_ed_with_line("SLC15A1", 
                        df%>% 
                          dplyr::rename(CAPRA_S = `CAPRA-S`), 
                        fit_table_purity_complete %>% 
                          rename(mean_value = mean)) +   
    labs(subtitle="Early Decrease") + 
    theme(panel.grid.major = element_line(colour='white', linetype = 1),
          panel.grid.minor= element_line(colour='white', linetype = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  
  #Decreasing Middle 
 DM<-plot_ed_with_line("PROK1", 
                    df%>% 
                      dplyr::rename(CAPRA_S = `CAPRA-S`), 
                    fit_table_purity_complete %>% 
                      rename(mean_value = mean))  + 
   labs(subtitle="Middle Decrease") + 
   theme(panel.grid.major = element_line(colour='white', linetype = 1),
         panel.grid.minor= element_line(colour='white', linetype = 1),
         panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  
  #Deacreasing Late
DL<-plot_ed_with_line("LGI3", 
                      df %>% 
                        dplyr::rename(CAPRA_S = `CAPRA-S`), 
                      fit_table_purity_complete %>% 
                        rename(mean_value = mean)) + 
  labs(subtitle="Late Decrease") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
  panel.grid.minor= element_line(colour='white', linetype = 1),
  panel.border = element_rect(colour = "black", fill=NA, size=0.5))


} 

{
library(ggplot2)
library(ggpubr)
}
  
norm_plot<- ggarrange(IE, IM, IL, DE, DM, DL, nrow =2, 
                      ncol =3, 
                      labels = c("I", "II", "III", "IV", "V", "VI")
) 



plot_B<-annotate_figure(norm_plot,
                         left = text_grob("Normalised Read Count Value", rot = 90, hjust = 0.5, x = 1, size = 16),
                         #fig.lab= "CAPRA-S Score",
                         bottom = text_grob("Scaled CAPRA-S score", size = 16),
                         right = text_grob("", size = 35)
                         # top = text_grob("B. Normalised Read Count vs Prostate Cancer Stage", hjust =1.1, face = "bold", size = 15)
)

# ggsave("TABI_examples_plot_2.pdf", width = 15, height = 10)





# Removing DE as it is problematic looking fit
blank<-ggplot(x=-2:2, y= 0:1000) + 
  theme_bw() +
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_blank())

unprob_norm_plot<-ggarrange(IE, IM, IL, blank, DM, DL, nrow =2, 
          ncol =3, 
          labels = c("I", "II", "III", "", "IV", "V")
) 

plot_B1<-annotate_figure(unprob_norm_plot,
                        left = text_grob("Normalised Read Count Value", rot = 90, hjust = 0.5, x = 1, size = 16),
                        #fig.lab= "CAPRA-S Score",
                        bottom = text_grob("Scaled CAPRA-S score", size = 16),
                        right = text_grob("", size = 35)
                        # top = text_grob("B. Normalised Read Count vs Prostate Cancer Stage", hjust =1.1, face = "bold", size = 15)
)


plot_B1

#ggsave("TABI_examples_DE_rm.pdf", width = 15, height =10)





