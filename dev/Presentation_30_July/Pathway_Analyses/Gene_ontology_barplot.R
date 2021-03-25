

#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #





# Intention of Script 

# Show Bargraph / quadrant plot of pathway analysis for Gene Ontology TCGA Analysis





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #





##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 




#Gene Ontology result plot

{library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
} 



#Number of results to show per quadrant
n_res<-10



#Load Table of gene ontology results 

GO_res<-read.csv("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Presentation_30_July/Pathway_Analyses/Linear_group_GO_results.csv") %>% 
  as_tibble() %>% 
  dplyr::rename(group = "Gene_Function_Set") %>% 
  # separate(group, c("group", "GO"), sep = "[(]") %>% 
  select(group, Direction_of_Expression_Change, 
         Region_of_Change, 
         Raw_P_value, 
         Method, 
         Fold_enrichment) %>% 
  filter(Method == "TABI")  %>% 
  group_by(Direction_of_Expression_Change,
           Region_of_Change) %>% 
  top_n(n=n_res, wt = Raw_P_value)

{
#Set up bargraphs for each quandrant
#Based on log_Raw_Pvalue

#Quadrant 2 Results (Increasing Late)
IL<-GO_res %>% 
    filter(Direction_of_Expression_Change == "Increase", 
           Region_of_Change == "Late") %>% 
    arrange(Raw_P_value) %>% 
    mutate(log_Raw_Pvalue = ifelse(Region_of_Change == "Early", 
                                   log(Raw_P_value),
                                   -log(Raw_P_value))) 


#Quadrant 1 Results (Increasing Early)
IE<-GO_res %>% 
  filter(Direction_of_Expression_Change == "Increase", 
         Region_of_Change == "Early") %>% 
  arrange(Raw_P_value) %>% 
  mutate(log_Raw_Pvalue = ifelse(Region_of_Change == "Early", 
                                 log(Raw_P_value), 
                                 -log(Raw_P_value))) 

#Quadrant 3 Results (Decreasing Early)
DE<-GO_res %>% 
  filter(Direction_of_Expression_Change == "Decrease", 
         Region_of_Change == "Early") %>% 
  arrange(Raw_P_value) %>% 
  mutate(log_Raw_Pvalue = ifelse(Region_of_Change == "Early", 
                                 log(Raw_P_value), 
                                 -log(Raw_P_value))) 


#Quadrant 4 Results (Decreasing Late)
DL<-GO_res %>% 
  filter(Direction_of_Expression_Change == "Decrease", 
         Region_of_Change == "Late") %>% 
  arrange(Raw_P_value) %>% 
  mutate(log_Raw_Pvalue = ifelse(Region_of_Change == "Early", 
                                 log(Raw_P_value), 
                                 -log(Raw_P_value))) 

  } 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #





#Qudrant plot which plots raw p-values pvalue 


ggplot() + 
  sapply(1:15, function(x)
    annotate("rect", xmin=0, xmax=IL$log_Raw_Pvalue[x], ymin=x-0.9 , ymax=x+0.1, alpha=0.2, color="red", fill="red")
  ) + 
  sapply(1:15, function(x)
    annotate("rect", xmin=0, xmax=IE$log_Raw_Pvalue[x], ymin=x-0.9 , ymax=x+0.1, alpha=0.2, color="blue", fill="blue")  
  ) + 
  sapply(1:15, function(x)
    annotate("rect", xmin=0, xmax=DE$log_Raw_Pvalue[x], ymin=-x+0.9 , ymax=-x-0.1, alpha=0.2, color="darkorange1", fill="orange")  
  ) + 
  sapply(1:15, function(x)
    annotate("rect", xmin=0, xmax=DL$log_Raw_Pvalue[x], ymin=-x+0.9 , ymax=-x-0.1, alpha=0.2, color="Dark Green", fill="green")  
  ) +
  sapply(1:15, function(x)
    annotate("text",x = 60, y= x-0.5, label = IL$group[x], size = 3, hjust = 1, colour = "firebrick")) + 
  sapply(1:15, function(x)
    annotate("text",x = -60, y= x-0.5, label = IE$group[x], size = 3, hjust = 0, colour = "dodgerblue4")) +
  sapply(1:15, function(x)
    annotate("text", x = -60, y= -x+0.5, label = DE$group[x], size = 3, hjust = 0, colour = "darkorange3"))  +
  sapply(1:15, function(x)
    annotate("text", x = 60, y= -x+0.5, label = DL$group[x], size = 3, hjust = 1, colour = "forestgreen")) + 
  annotate("text", x = -10, y= 1, label = "Direction of Gene Expression Change", size = 5.5, fontface = "italic", angle = 90, hjust = 0) + 
  annotate("text", x = 10, y= -16, label = "Stage of Prostate Cancer", size = 5.5, fontface = "italic", hjust = 0) + 
  annotate("text",x = 60, y= 17, label = "Late Increasing", size = 6, fontface = 4, hjust = 1)  + 
  annotate("text", x = 60, y= -17, label = "Late Decreasing", size = 6, fontface = 4, hjust = 1)  + 
  annotate("text", x = -60, y= 17, label = "Early Increasing", size = 6, fontface = 4, hjust = 0)  + 
  annotate("text",x = -60, y= -17, label = "Early Decreasing", size = 6, fontface = 4, hjust = 0)  + 
  annotate("segment", x = -19, xend = 19, y = 0, yend = 0, colour = "black", size=1.5, alpha=1, arrow=arrow()) + 
  annotate("segment", x = -60, xend = 60, y = 0, yend = 0, colour = "black", size=0.5, alpha=1, linetype = 2) +  
  annotate("segment", x = 0, xend = 0, y = -20, yend = 20, colour = "black", size=0.5, alpha=1, linetype = 2) +
  annotate("segment", x = 19, xend = -19, y = 0, yend = 0, colour = "black", size=1.5, alpha=1, arrow=arrow()) +
  annotate("segment", x = 0, xend = 0, y = -16, yend = 16, colour = "black", size=1.5, alpha=1, arrow=arrow()) + 
  annotate("segment", x = 0, xend = 0, y = 16, yend = -16, colour = "black", size=1.5, alpha=1, arrow=arrow()) +
  annotate("text", x = 5, y= -1.8, label = "-log(unadjusted pvalue)", size = 3, hjust = 0) +
  labs(y= "", x = "") + 
  sapply(seq(from = -15, to = 15, by =2.5), function(x) 
    annotate("text", x = x, y= -0.75, label = paste0(abs(x)), size = 3, hjust = 0)) + 
  sapply(seq(from = -15, to = 15, by =2.5), function(x)
    annotate("segment", x = x, xend = x, y = 0, yend = -0.25, colour = "black", size=0.5, alpha=1)) + 
  xlim(-60, 60) + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



# Quadrant Plot with gene enrichment measured 


# Set up rank by gene enrichment for plotting

{ DL<-DL %>% 
    arrange(desc(Fold_enrichment)) %>% 
    mutate(Fold_enrichment = as.numeric(Fold_enrichment)) 
  
  DE<- DE %>% 
    arrange(desc(Fold_enrichment)) %>% 
    mutate(Fold_enrichment = -1*as.numeric(Fold_enrichment))
  
  IL<-IL %>% 
    arrange((Fold_enrichment)) %>% 
    mutate(Fold_enrichment = as.numeric(Fold_enrichment))
    
  
  IE<-IE %>% 
    arrange((Fold_enrichment)) %>% 
    mutate(Fold_enrichment = -1*as.numeric(Fold_enrichment)) 
  }


# Plot gene enrichement pattern plot

ggplot() + 
  sapply(1:15, function(x)
    annotate("rect", xmin=0, xmax= IL$Fold_enrichment[x], ymin=x-0.9 , ymax=x+0.1, alpha=0.2, color="red", fill="red")
  ) + 
  sapply(1:15, function(x)
    annotate("rect", xmin=0, xmax=IE$Fold_enrichment[x], ymin=x-0.9 , ymax=x+0.1, alpha=0.2, color="blue", fill="blue")  
  ) + 
  sapply(1:15, function(x)
    annotate("rect", xmin=0, xmax=DE$Fold_enrichment[x], ymin=-x+0.9 , ymax=-x-0.1, alpha=0.2, color="darkorange1", fill="orange")  
  ) + 
  sapply(1:15, function(x)
    annotate("rect", xmin=0, xmax=DL$Fold_enrichment[x], ymin=-x+0.9 , ymax=-x-0.1, alpha=0.2, color="Dark Green", fill="green")  
  ) +
  sapply(1:15, function(x)
    annotate("text",x = 1300, y= x-0.5, label = IL$group[x], size = 5, hjust = 1, colour = "firebrick")) + 
  sapply(1:15, function(x)
    annotate("text",x = -1300, y= x-0.5, label = IE$group[x], size = 5, hjust = 0, colour = "dodgerblue4")) +
  sapply(1:15, function(x)
    annotate("text", x = -1300, y= -x+0.5, label = DE$group[x], size = 5, hjust = 0, colour = "darkorange3"))  +
  sapply(1:15, function(x)
    annotate("text", x = 1300, y= -x+0.5, label = DL$group[x], size =5, hjust = 1, colour = "forestgreen")) + 
  annotate("text", x = -50, y= 1, label = "Direction of Gene Expression Change", size = 7, fontface = "italic", angle = 90, hjust = 0) + 
  annotate("text", x = 50, y= -11, label = "Stage of Prostate Cancer", size = 7, fontface = "italic", hjust = 0) + 
  annotate("text",x = 950, y= 13, label = "Late Increasing", size = 8, fontface = 4, hjust = 1)  + 
  annotate("text", x =950, y= -13, label = "Late Decreasing", size = 8, fontface = 4, hjust = 1)  + 
  annotate("text", x = -950, y= 13, label = "Early Increasing", size = 8, fontface = 4, hjust = 0)  + 
  annotate("text",x = -950, y= -13, label = "Early Decreasing", size = 8, fontface = 4, hjust = 0)  + 
  annotate("segment", x = -550, xend = 550, y = 0, yend = 0, colour = "black", size=1.5, alpha=1, arrow=arrow()) + 
  annotate("segment", x = -1300, xend = 1300, y = 0, yend = 0, colour = "black", size=0.5, alpha=1, linetype = 2) +  
  annotate("segment", x = 0, xend = 0, y = -13, yend = 13, colour = "black", size=0.5, alpha=1, linetype = 2) +
  annotate("segment", x = 550, xend = -550, y = 0, yend = 0, colour = "black", size=1.5, alpha=1, arrow=arrow()) +
  annotate("segment", x = 0, xend = 0, y = -11, yend = 12, colour = "black", size=1.5, alpha=1, arrow=arrow()) + 
  annotate("segment", x = 0, xend = 0, y = 12, yend = -11, colour = "black", size=1.5, alpha=1, arrow=arrow()) +
  annotate("text", x = 25, y= 0.5, label = "Fold Enrichment", size = 5, hjust = 0) +
  labs(y= "", x = "") + 
  sapply(seq(from = -550, to = 550, by =100), function(x) 
    annotate("text", x = x, y= -0.75, label = paste0(abs(x)/100), size = 4, hjust = 1)) + 
  sapply(seq(from = -550, to = 550, by =50), function(x)
    annotate("segment", x = x, xend = x, y = 0, yend = -0.25, colour = "black", size=0.5, alpha=1)) + 
  xlim(-1300, 1300) + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  annotate("text",x = -950, y= 15.5, label = "Pattern of Enriched Gene Ontology Classes in Differential Gene Expression", size = 10, fontface = 2, hjust = 0.25)




# Save plot
# Plot saved is that of top 10 signficant / plotted with gene enrichement
# I.e. of above code
ggsave("GO_quad_pat_consistent_sig.pdf", width = 25, height = 13)







