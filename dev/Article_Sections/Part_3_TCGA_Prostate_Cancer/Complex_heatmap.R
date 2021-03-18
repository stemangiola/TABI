
box::use(dplyr[...],
         tidyr['seperate'],
         ComplexHeatmap[...],
         simplifyEnrichment[...])

# Load panther GO enrichment Results (DE genes, )


# Filter Panther GO Enrichment Results by Stage (Early) and Direction of Change (Increase)
GO_ID_ear_up= Linear_group_GO_results %>% 
  filter(
    Direction_of_Expression_Change == "Increase"&Region_of_Change == "Early") %>% 
    select(Gene_Function_Set)  %>% 
    tidyr::separate(col = Gene_Function_Set, 
                             into = c(NA, "GO_ID"), 
                                             sep = "[()]")

# 
Q1 = GO_similarity(GO_ID_ear_up$GO_ID)

# Set up complex heatmap options
# Adding black borders / setting font size
ht_opt(
         heatmap_column_names_gp = gpar(fontface = "italic"), 
                 heatmap_column_title_gp = gpar(fontsize = 24),
                 legend_border = "black",
                 heatmap_border = "black",
                 annotation_border = "black")


# Setting range of colours for word cloud in heatmap
col_fun = colorRamp2(c(8, 12, 16), 
                        viridis(n = 4)[1:3])

plot_q1 = grid.grabExpr(ht_clusters(Q1, 
                               cluster_terms(Q1), 
                               column_title = "Early Increase",
                               order_by_size = TRUE, 
                               col = magma(n = 50),  
                               dend = TRUE,
                               depth = 1,
                               word_cloud_grob_param = list(
                                 col = function(fs) col_fun(fs)), 
                               fontsize_range = c(8, 16)))



GO_ID_ear_dec= Linear_group_GO_results %>% 
  filter(
    Direction_of_Expression_Change == "Decrease"&Region_of_Change == "Early") %>% 
  select(Gene_Function_Set)  %>% 
  tidyr::separate(col = Gene_Function_Set, 
                                   into = c(NA, "GO_ID"), 
                                                  sep = "[()]")


Q2 = GO_similarity(GO_ID_ear_dec$GO_ID)


plot_q2 = grid.grabExpr(ht_clusters(Q2, 
                                    cluster_terms(Q2), 
                                    column_title = "Early Decrease",
                                    order_by_size = TRUE, 
                                    col = magma(n = 50),  
                                    word_cloud_grob_param = list(
                                      col = function(fs) col_fun(fs)), 
                                    fontsize_range = c(8, 16)))


GO_ID_late_up= Linear_group_GO_results %>% 
  filter(
    Direction_of_Expression_Change == "Increase"&Region_of_Change == "Late") %>% 
    select(Gene_Function_Set)  %>% 
    tidyr::separate(col = Gene_Function_Set, 
                                       into = c(NA, "GO_ID"), 
                                                        sep = "[()]")



Q3 = GO_similarity(GO_ID_late_up$GO_ID)


plot_q3 = grid.grabExpr(ht_clusters(Q3, 
                                    cluster_terms(Q3), 
                                    column_title = "Late Increase",
                                    order_by_size = TRUE, 
                                    col = magma(n = 50),  
                                    word_cloud_grob_param = list(
                                      col = function(fs) col_fun(fs)), 
                                    fontsize_range = c(8, 16)))



GO_ID_late_up= Linear_group_GO_results %>% 
  filter(
    Direction_of_Expression_Change == "Decrease"&Region_of_Change == "Late") %>% 
  select(Gene_Function_Set)  %>% 
  tidyr::separate(col = Gene_Function_Set, 
                  into = c(NA, "GO_ID"), 
                  sep = "[()]")

Q4 = GO_similarity(GO_ID_late_dec$GO_ID)

plot_q4 = grid.grabExpr(ht_clusters(Q4, 
                                    cluster_terms(Q4), 
                                    column_title = "Late Decrease",
                                    order_by_size = TRUE, 
                                    col = magma(n = 50),  
                                    word_cloud_grob_param = list(
                                      col = function(fs) col_fun(fs)), 
                                    fontsize_range = c(8, 16)))


ggarrange(plot_q1,
          plot_q3,
          plot_q2,
          plot_q4,
          ncol = 2,
          nrow = 2)





