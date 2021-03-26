


Linear_group_GO_results <-read.csv("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Linear_group_GO_results.csv")


box::use(dplyr[...],
         ComplexHeatmap[...],
         simplifyEnrichment[...])


library(dplyr)
library(viridis)
library(magick)
library(circlize)
library(simplifyEnrichment)
library(ComplexHeatmap)
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
ht_opt(heatmap_column_names_gp = gpar(size = 7,
  fontface = "italic"), 
                 legend_border = "black",
                 heatmap_border = "black",
                 annotation_border = "black")


# Setting range of colours for word cloud in heatmap
col_fun = colorRamp2(c(4, 6, 8), 
                        viridis(n = 4)[1:3])

# grid.grabExpr(

plot_q1 = grid.grabExpr(simplifyEnrichment::ht_clusters(Q1, 
                               cluster_terms(Q1), 
                               column_title = "Early Increase",
                               order_by_size = TRUE, 
                               col = magma(n = 50),  
                               min_term = 15,
                               max_words = 5,
                               show_heatmap_legend = FALSE,
                               bg_gp = gpar(fill = "#FFFFFF", col = "#AAAAAA"),
                               column_title_gp = gpar(fontsize = 6),
                               word_cloud_grob_param = list(
                                 line_space = 0.1,
                                 word_space = 0.5,
                                 max_width = 25,
                                 col = function(fs) col_fun(fs)), 
                               fontsize_range = c(4, 8))) 





GO_ID_ear_dec= Linear_group_GO_results %>% 
  filter(
    Direction_of_Expression_Change == "Decrease"& Region_of_Change == "Early") %>% 
  dplyr::select(Gene_Function_Set)  %>% 
  tidyr::separate(col = Gene_Function_Set, 
                                   into = c(NA, "GO_ID"), 
                                                  sep = "[()]")


Q2 = GO_similarity(GO_ID_ear_dec$GO_ID)


plot_q2 = grid.grabExpr(simplifyEnrichment::ht_clusters(Q2, 
                                    cluster_terms(Q2), 
                                    column_title = "Early Decrease",
                                    order_by_size = TRUE, 
                                    col = magma(n = 50),  
                                    min_term = 15,
                                    max_words = 5,
                                    show_heatmap_legend = FALSE,
                                    heatmap_legend_side = "bottom",
                                   # heatmap_legend_param = list(legend_direction = "horizontal"),
                                    column_title_gp = gpar(fontsize = 8),
                                    bg_gp = gpar(fill = "#FFFFFF", col = "#AAAAAA"),
                                    word_cloud_grob_param = list(
                                      line_space = 0.1,
                                      word_space = 0.5,
                                      max_width = 25,
                                      col = function(fs) col_fun(fs)), 
                                    fontsize_range = c(4, 8)))


GO_ID_late_up= Linear_group_GO_results %>% 
  filter(
    Direction_of_Expression_Change == "Increase"&Region_of_Change == "Late") %>% 
    dplyr::select(Gene_Function_Set)  %>% 
    tidyr::separate(col = Gene_Function_Set, 
                                       into = c(NA, "GO_ID"), 
                                                        sep = "[()]")



Q3 = GO_similarity(GO_ID_late_up$GO_ID)


plot_q3 = grid.grabExpr(simplifyEnrichment::ht_clusters(Q3, 
                                    cluster_terms(Q3), 
                                    column_title = "Late Increase",
                                    order_by_size = TRUE, 
                                    col = magma(n = 50),  
                                    min_term = 15,
                                    max_words = 5,
                                    show_heatmap_legend = FALSE,
                                    column_title_gp = gpar(fontsize = 8),
                                    bg_gp = gpar(fill = "#FFFFFF", col = "#AAAAAA"),
                                    word_cloud_grob_param = list(
                                      line_space = 0.1,
                                      word_space = 0.5,
                                      max_width = 25,
                                      col = function(fs) col_fun(fs)), 
                                    fontsize_range = c(4, 8)))



GO_ID_late_up= Linear_group_GO_results %>% 
  filter(
    Direction_of_Expression_Change == "Decrease"&Region_of_Change == "Late") %>% 
  dplyr::select(Gene_Function_Set)  %>% 
  tidyr::separate(col = Gene_Function_Set, 
                  into = c(NA, "GO_ID"), 
                  sep = "[()]")

Q4 = GO_similarity(GO_ID_late_up$GO_ID)

  plot_q4 = grid.grabExpr(simplifyEnrichment::ht_clusters(Q4, 
                      cluster_terms(Q4), 
                                    column_title = "Late Decrease",
                                    column_title_gp = gpar(fontsize = 6),
                                    order_by_size = TRUE, 
                                    col = magma(n = 50),
                                    show_heatmap_legend = FALSE,
                                    min_term = 15,
                                    max_words = 5,
                                    fontsize_range = c(4,8),
                                    bg_gp = gpar(fill = "#FFFFFF",
                                                 col = "#AAAAAA",
                                                 lineheight = 200),
                                    word_cloud_grob_param = list(
                                      line_space = 0.1,
                                      word_space = 0.5,
                                      max_width = 25,
                                      col = function(fs) col_fun(fs)))) 
ht_clusters(Q3, 
               cluster_terms(Q3), 
               column_title = "Late Increase",
              column_title_gp = gpar(fontsize = 12),
               order_by_size = TRUE, 
               col = magma(n = 50),
               min_term = 10,
               max_words = 5,
             # dend = TRUE,
              show_heatmap_legend = FALSE,
              show_annotation_legend = FALSE, 
              #heatmap_legend_side = "bottom",
              #align_heatmap_legend = "heatmap_bottom",
               bg_gp = gpar(fill = "#FFFFFF", col = "#AAAAAA", lineheight = 200),
               word_cloud_grob_param = list(
                 line_space = 0.1,
                 word_space = 
                 max_width = 25,
                 col = function(fs) col_fun(fs)), 
               fontsize_range = c(4, 10))


ht_list = plot_q1 + plot_q2 

# + plot_q3 + plot_q4

complex_heat_map = ggpubr::ggarrange(plot_q1,
          plot_q3,
          plot_q2,
          plot_q4,
          #labels = c("D", "", "", ""),
          font.label = list(size = 12),
          ncol = 2,
          nrow = 2)


ggplot2::ggsave(filename = "complex_heat_map.pdf", width=183, height=238/2.5, units = "mm")




