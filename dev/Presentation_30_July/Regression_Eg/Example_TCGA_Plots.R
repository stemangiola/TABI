


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
                        dplyr::select(mean_value))
    
    
    beta_1 <- as.numeric(Gene_table %>%
                           filter(parameters == "beta1_z[1,1]") %>%
                           dplyr::select(mean_value))
    
    # beta_2<-as.numeric(Gene_table %>%
    #                      filter(parameters == "beta1_z[2,1]") %>%
    #                      dplyr::select(mean_value))
    
    #beta_3
    
    y_0 <- as.numeric(Gene_table %>%
                        filter(parameters == "y_cross[1]") %>%
                        dplyr::select(mean_value))
    
    A<- as.numeric(Gene_table %>%
                     filter(parameters == "A[1]") %>%
                     dplyr::select(mean_value))
    
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
      dplyr::select(`(Intercept)`, one_of(parse_formula(formula)))
  }
  
  
  #Create table with scaled CAPRA_S and normalised gene count for graphing
  graphing_table <- data.frame( #Scale CAPRA_S using TABI method 
    CAPRA_S = as.vector(model.matrix(object = ~CAPRA_S, data = read_count_data %>% 
                                       filter(is.na(CAPRA_S)==F, is.na(transcript)==F) %>% 
                                       filter(transcript ==gene_name)) %>%
                          as_tibble(rownames="sample_idx") %>%
                          scale_design(~CAPRA_S) %>% dplyr::select(CAPRA_S)),
    #Combine with read counts 
    read_count_normalised = as.vector(read_count_data %>% filter(is.na(CAPRA_S)==F) %>% 
                                        filter(transcript ==gene_name) %>% 
                                        dplyr::select(`read count normalised`)) )
  
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
    #ylim(c(200, 3500)) + 
    labs(
      #title=paste(gene_name),
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


p4<-df %>% filter(transcript == "C12orf49") %>% 
  dplyr::rename(CAPRA_S = `CAPRA-S`) %>% 
  filter(is.na(CAPRA_S)==F) %>% 
  rowwise() %>% 
  dplyr::rename(value = 'read count normalised') %>% 
  mutate(Group = ifelse(CAPRA_S<=1, 1, ifelse(CAPRA_S>1&CAPRA_S<=3, 2, ifelse(CAPRA_S>=6, 4, 3)) )) %>% 
  ggplot(aes(x=Group, y=value)) + 
  geom_point(col = "dodgerblue") + 
  scale_y_log10(breaks=c(1000,2000,3000, 4000)) +
  scale_x_continuous(labels=scaleFUN) + 
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black") +  
  labs(x = "", y= "") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  labs(title = "Four Groups Approach") + 
  custom_theme



p2<-df %>% filter(transcript == "C12orf49") %>% 
  dplyr::rename(CAPRA_S = `CAPRA-S`) %>% 
  rowwise() %>% 
  filter(is.na(CAPRA_S)==F) %>% 
  dplyr::rename(value = 'read count normalised') %>% 
  mutate(Group = ifelse(CAPRA_S<3.5, 1, 2)) %>% 
  mutate(Group = as.double(Group)) %>% 
  ggplot(aes(x=Group, y=value)) + 
  geom_point(col = "dodgerblue") + 
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black") +  
  scale_y_log10(breaks=c(1000,2000,3000, 4000)) +
  scale_x_continuous(labels=scaleFUN) + 
  labs(x = "", y= "") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  labs(title = "Two Groups Approach") + 
  custom_theme 


scaleFUN <- function(x) sprintf("%.f", x)



p3<-df %>% filter(transcript == "C12orf49") %>% 
  dplyr::rename(CAPRA_S = `CAPRA-S`) %>% 
  rowwise() %>% 
  filter(is.na(CAPRA_S)==F) %>% 
  dplyr::rename(value = 'read count normalised') %>% 
  mutate(Group = ifelse(CAPRA_S<=2, 1, ifelse(CAPRA_S<4, 2, 3))) %>% 
  ggplot(aes(x=Group, y=value)) + 
  geom_point(col = "dodgerblue") + 
  scale_y_log10(breaks=c(1000,2000,3000, 4000)) + 
  scale_x_continuous(labels=scaleFUN) + 
  #scale_y_log10(limits = c(1, 4000)) + 
  geom_smooth(aes(Group), formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(Group<1.5)) +  
  geom_smooth(aes(Group), formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(Group>=1)) + 
  labs(x = "", y= "") + 
  labs(title = "Three Groups Approach") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  custom_theme 



p7<-df %>% filter(transcript == "SPRY4") %>% 
  dplyr::rename(CAPRA_S = `CAPRA-S`) %>% 
  rowwise() %>% 
  filter(is.na(CAPRA_S)==F) %>% 
  rename(value = 'read count normalised') %>% 
  ggplot(aes(x=CAPRA_S, y=value)) + 
  geom_point(col = "dodgerblue") + 
  scale_y_log10() + 
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(CAPRA_S>6)) +  
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(CAPRA_S<=7) %>% filter(CAPRA_S>=6)) +  
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(CAPRA_S<=6) %>% filter(CAPRA_S>=5)) + 
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(CAPRA_S<=5) %>% filter(CAPRA_S>=4)) +  
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(CAPRA_S<=4) %>% filter(CAPRA_S>=3)) + 
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(CAPRA_S<=3) %>% filter(CAPRA_S>=2)) +
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(CAPRA_S<=2) %>% filter(CAPRA_S>=1)) + 
  geom_smooth(formula = y ~ x, method = "lm", se =FALSE, col = "black", data = . %>% filter(CAPRA_S<=1) %>% filter(CAPRA_S>=0)) + 
  labs(x = "", y= "") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))


IE
friendly_cols <- dittoSeq::dittoColors()

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


{ # Increasing (Top Row)
  
  #Increasing Early 
  
  #Alternative Increasing Early - SPRY4 more dramatic trend, but not as early as C12orf49
IE<-plot_ed_with_line("C12orf49", 
                      df %>% 
                        dplyr::rename(CAPRA_S = `CAPRA-S`), 
                      fit_table_purity_complete %>% 
                        dplyr::rename(mean_value = mean))  + 
    labs(title="TABI Approach", subtitle = "C12orf49") + 
    theme(panel.grid.major = element_line(colour='white', linetype = 1),
          panel.grid.minor= element_line(colour='white', linetype = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
    custom_theme + 
    scale_y_log10(breaks=c(1000,2000,3000, 4000)) 
  
  ggarrange(
  annotate_figure(
    ggarrange(IE, 
            p2, 
            p3, 
            p4,
            nrow = 2,
            ncol = 2)), 
    bottom = text_grob("Scaled CAPRA-S Score", size = 10),
    left = text_grob("Normalised Transcript Abundance", size = 10, rot = 90),
  ggplot() + theme_void(),
  nrow = 2,
  heights = c(1,2)) 
  
  

IE<-plot_ed_with_line("SPRY4",
                      df%>%
                        dplyr::rename(CAPRA_S = `CAPRA-S`),
                      fit_table_purity_complete %>%
                        dplyr::rename(mean_value = mean)) +
  labs(subtitle="Early Decrease", title = "SPRY4") +
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  custom_theme

 #Middle Increase
IM<-plot_ed_with_line("IQGAP3", 
                      df%>% 
                        dplyr::rename(CAPRA_S = `CAPRA-S`), 
                      fit_table_purity_complete %>% 
                        dplyr::rename(mean_value = mean))  + 
  labs(subtitle="Middle Increase", title = "IQGAP3") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  custom_theme

# Late Increase

IL<-plot_ed_with_line("CXCR5", 
                      df%>% 
                        dplyr::rename(CAPRA_S = `CAPRA-S`), 
                      fit_table_purity_complete %>% 
                        dplyr::rename(mean_value = mean))  + 
  labs(subtitle="Late Increase", title = "CXCR5") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  custom_theme

} 

{# Decreasing
  
  #Decreasing Early
DE<-plot_ed_with_line("SLC15A1", 
                        df%>% 
                          dplyr::rename(CAPRA_S = `CAPRA-S`), 
                        fit_table_purity_complete %>% 
                          dplyr::rename(mean_value = mean)) +   
    labs(subtitle="Early Decrease", title = "SLC15A1") + 
    theme(panel.grid.major = element_line(colour='white', linetype = 1),
          panel.grid.minor= element_line(colour='white', linetype = 1),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
    custom_theme
  
  #Decreasing Middle 
 DM<-plot_ed_with_line("PROK1", 
                    df%>% 
                      dplyr::rename(CAPRA_S = `CAPRA-S`), 
                    fit_table_purity_complete %>% 
                      dplyr::rename(mean_value = mean))  + 
   labs(subtitle="Middle Decrease", title = "PROK1") + 
   theme(panel.grid.major = element_line(colour='white', linetype = 1),
         panel.grid.minor= element_line(colour='white', linetype = 1),
         panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
   custom_theme
  
  #Deacreasing Late
DL<-plot_ed_with_line("LGI3", 
                      df %>% 
                        dplyr::rename(CAPRA_S = `CAPRA-S`), 
                      fit_table_purity_complete %>% 
                        dplyr::rename(mean_value = mean)) + 
  labs(subtitle="Late Decrease", title = "LG13") + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
  panel.grid.minor= element_line(colour='white', linetype = 1),
  panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  custom_theme


} 

{
library(ggplot2)
library(ggpubr)
}

plot_examples= ggarrange(IE, IM, IL, DE, DM, DL, nrow =2, 
          ncol =3, 
          labels = c("I", "II", "III", "IV", "V", "VI"))
  
  norm_plot<- annotate_figure(
  plot_examples,
            #bottom = text_grob("Scaled CAPRA-S Risk Score", hjust = -0.75, vjust = 20, size = 10),
            left = text_grob("Normalised Transcript Abundance", 
                             rot = 90, 
                             vjust = 0, #-35, 
                             hjust = 0.85, 
                             size = 10)
) 

  ggarrange(ggplot() + theme_void(),
            norm_plot,
            complex_heat_map,
            nrow = 3,
            labels = c("A", "B", "C"),
            heights  = c(1, 1.5,1.25))


p3<-p3 %>% labs(title = "Three Groups")
p2<-p2 %>% labs(title = "Two Groups")
p4<-p4 %>% labs(title = "Four Groups")

test<-ggarrange(p2 , p3, p4,
          ncol =3, 
          labels = c("III.", "IV.", "V"),
          vjust = 1
) 

I<-ggarrange(IE, eg, nrow =1, vjust = 1, labels = c("I.", "II. "))

IE<-IE %>% 
  labs(title = "TABI Approach", 
       subtitle  = "SPRY4") 
  
  

t1<-ggarrange(I, test, 
             # label.x = 0.9,
              nrow = 2) 

annotate_figure(t1,
  top = text_grob("Figure 2. Underlying Methodolgy", size = 16)
)


ggsave("TABI_test.pdf", width = 22, height =10)

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





