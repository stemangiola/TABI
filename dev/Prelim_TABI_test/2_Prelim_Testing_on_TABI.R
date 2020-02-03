library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

plot_1 <- fit_p1 %>% #Data Table 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  select(parameters, "mean", Gene_name) %>% 
  pivot_wider(names_from = parameters, values_from = "mean")%>%
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]") %>% 
  ggplot(aes(x=Inflection_mean, y= Slope_mean)) +
  geom_point() + 
  stat_smooth(se = FALSE, method = 'lm', col="red") +
  labs(title = "CAPRA_S 8 as separate group") +
  ylim(-2.5, 4.5) +
  xlim(-3, 3)

  


plot_2 <-fit_ed_p1 %>% #Data Table 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  select(parameters, "mean", Gene_name) %>% 
  pivot_wider(names_from = parameters, values_from = "mean")%>%
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]") %>% 
  ggplot(aes(x=Inflection_mean, y= Slope_mean)) +
  geom_point() + 
  stat_smooth(se = FALSE, method = 'lm', col="red") +
  labs(title = "CAPRA_S 8 placed in group CAPRA_S 7") +
  ylim(-2.5, 4.5) +
  xlimm(-3, 3)

grid.arrange(plot_1, plot_2, nrow=2)

#Plotting Function with paramters




fit_p1 <- fit_p1 %>% rename(mean_value = mean)

plot_ed_with_line <- function(gene_name, data, fit_table) {
  
  #Required packages / functions 
  require(ggplot2)
  require(dplyr)
  require(plotly)
  #Requires scale_design and parse_formula from TABI to create scaled design matrix 
  
  #Select only paramters for specified gene
  Gene_table <- fit_table %>% 
    filter(Gene_name == gene_name)
  
  #Stanfit (for converiting to fit table)
  #fit_table<-as.data.frame(summary( #Summarise stanfit object
    # Call Stanfit file
    #stan_fit)$summary) %>% #Extract summary of those 
    #rownames_to_column("parameters") %>% #Make row names (which are names of paramters) explicit column 
    #as_tibble() %>% 
    #filter(!grepl("y_hat", parameters), !grepl("y_gen", parameters), !grepl("phi", parameters)) %>% #Remove rows with y_hat and y_gen 
    #mutate(Gene_name = gene_name)   

plot_line <- function(x) {
  
  
  eta <- as.numeric(Gene_table %>%
                      filter(parameters == "inflection[1]") %>% 
                      select(mean_value))
                      
  
  beta <- as.numeric(Gene_table %>%
                       filter(parameters == "beta[2,1]") %>%
                       select(mean_value))
  
  y_0 <- as.numeric(Gene_table %>%
                      filter(parameters == "y_cross[1]") %>%
                      select(mean_value))
  
  A<- as.numeric(Gene_table %>%
                   filter(parameters == "A[1]") %>%
                   select(mean_value))
  
  top<-(y_0-A)*(1+exp(eta*beta))
  
  bottom<-(1+exp(-x*beta+eta*beta))
  
  return(exp(A+(top/bottom)))}

#Create table with scaled CAPRA_S and normalised gene count for graphing
graphing_table <- data.frame( #Scale CAPRA_S using TABI method 
  CAPRA_S = as.vector(model.matrix(object = ~CAPRA_S, data = data %>% filter(is.na(CAPRA_S)==F, is.na(transcript)==F) %>% filter(transcript ==gene_name)) %>%
                        as_tibble(rownames="sample_idx") %>%
                        scale_design(~CAPRA_S) %>% select(CAPRA_S)),
  #Combine with read counts 
  read_count_normalised = as.vector(data %>% filter(is.na(CAPRA_S)==F) %>% filter(transcript ==gene_name) %>% select(`read count normalised`)) )

#Plot curve of specified parameters and of count points
graphing_table %>% 
  ggplot(aes(x=CAPRA_S, y=read.count.normalised)) + 
  geom_point(col="dodgerblue") +
  stat_function(fun=plot_line, geom="line")

#ggplotly(plot)
}

to_save <-normalised_PC_TCGA %>% select(`read_count normalised`, transcript, purity.score, CAPRA_S)
saveRDS(to_save, compress = "gzip", file="download_PC")


plot_ed_with_line("CDC20B", normalised_PC_TCGA, fit_ed_p1) + geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

compressed_TCGA %>%
  ggplot(aes(`read count normalised` +1, group=Sample.ID)) + geom_density() + scale_x_log10()


library(tidyr)

library(dplyr)
               



