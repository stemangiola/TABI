
#Comparing TABI with alternative method
# Grouping and then performing linear regression analysis 

library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)
library(gridExtra)

# TCGA set which has had gene counts normalised by tidyBulk/ na rm
compress_normalised_narm_TCGA <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/compress_normalised_narm_TCGA.rds")


# With Grouping Into 4 groups
# Define groups by CAPRA_S score (0-1 = Group_A, 2-3 = Group B, 4-5 = Group C, 6-8 = Group D)
Group_test_table<-compress_normalised_narm_TCGA %>%
  mutate(CAPRA_S_group = ifelse(CAPRA_S>=0&CAPRA_S<=1, "A"
                                , ifelse(CAPRA_S>=2&CAPRA_S<=3,
                                         "B",
                                         ifelse(CAPRA_S>=4&CAPRA_S<=5, "C",
                                                "D"))))  %>%
  select(CAPRA_S_group, 
         `read_count normalised`, 
         purity.score, 
         transcript)%>%
  rename(`read_count_normalised` = `read_count normalised`) %>%
  dplyr::mutate_if(is.character, 
                   as.factor) 

# Linear Equation testing differential expression Group A / B
test_AB<-Group_test_table %>%
  filter(CAPRA_S_group=="A"|CAPRA_S_group=="B") %>%
  group_by(transcript) %>%
  summarise(p_val = (tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,5],
            CAPRA_S_slope = (tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,4])

# Tranform result into a usable form
test_AB_trans<-data.frame(pval_AB = t(test_AB$p_val), 
                          slope_AB =t(test_AB$CAPRA_S_slope), 
                          transcript = test_AB$transcript)

# Linear Equation testing differential expression Group B / C
test_BC<-Group_test_table %>%
  filter(CAPRA_S_group=="B"|CAPRA_S_group=="C") %>%
  group_by(transcript) %>%
  summarise(CAPRA_S_p_val = (broom::tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,5],
            CAPRA_S_slope = (broom::tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,4])


# Tranform result into a usable form
test_BC_trans<-data.frame(pval_BC = t(test_BC$CAPRA_S_p_val), 
                          slope_BC =t(test_BC$CAPRA_S_slope), 
                          transcript = test_BC$transcript)

# Linear Equation testing differential expression Group C / D
test_CD<-Group_test_table %>%
  filter(CAPRA_S_group=="C"|CAPRA_S_group=="D") %>%
  group_by(transcript) %>%
  summarise(CAPRA_S_p_val = (broom::tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,5],
            CAPRA_S_slope = (broom::tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,4])


# Tranform result into a usable form
test_CD_trans<-data.frame(pval_CD = t(test_CD$CAPRA_S_p_val), 
                          slope_CD =t(test_CD$CAPRA_S_slope), 
                          transcript = test_CD$transcript)



#Calculate / select a single pvalue for each gene - by the below method 
#significance threshold is taken as 0.05
P_value_group_4<-all_test_group_4 %>% 
  group_by(transcript) %>%  #For each gene
  summarise(pval_total = ifelse(
    (pval_AB>=0.05)&(pval_BC>=0.05)&(pval_CD>=0.05), # If all tests are non-sig
    sample(c(pval_AB, pval_BC, pval_CD), size =1),  #Take a random sample of associated P-values
    ifelse(
      (pval_AB<0.05)&(pval_BC>=0.05)&(pval_CD>=0.05), #If Only A/B is sig
      pval_AB,                                         #Then select the pvalue for A/B
      ifelse(
        (pval_AB>=0.05)&(pval_BC<0.05)&(pval_CD>=0.05), #If Only B/C is sig 
        pval_BC,                                         #Then select the pvalue for B/C 
        ifelse(
          (pval_AB>=0.05)&(pval_BC>=0.05)&(pval_CD<0.05), #If Only C/D is sig
          pval_CD,                                           #Then select the pvalue for C/D
          ifelse(
            (pval_AB<0.05)&(pval_BC<0.05)&(pval_CD>=0.05)&(sign(pval_AB)==sign(pval_BC)),  #If Two slopes are sig (AB BC) and same direction
            sample(c(pval_AB, pval_BC)),                                                  #Take a random sample of the two sig pvalues 
            ifelse(
              (pval_AB>=0.05)&(pval_BC<0.05)&(pval_CD<0.05)&(sign(pval_BC)==sign(pval_CD)),  #Two slopes are sig (BC CD) and same direction
              sample(c(pval_BC, pval_CD)),                                                   #Take a random sample of the two sig pvalues 
              ifelse(
                (pval_AB<0.05)&(pval_BC>=0.05)&(pval_CD<0.05)&(sign(pval_AB)==sign(pval_CD)),  #Two slopes are sig (AB CD) and same direction
                sample(c(pval_AB, pval_CD)),                                                    #Take a random sample of the two sig pvalues 
                ifelse(
                  (pval_AB<0.05)&(pval_BC<0.05)&(pval_CD>=0.05)&(sign(pval_AB)==sign(pval_BC)),  #Two slopes are sig (AB BC) and same direction
                  sample(c(pval_AB, pval_BC)),                                                    #Take a random sample of the two sig pvalues 
                  sample(c(pval_AB, pval_BC, pval_CD))                                              #Otherwise take a random sample of associated Pvalues
                ) 
              ) 
            ) 
          ) 
        )
      )
    )
  )
  )







test_AB_trans %>% 
  ggplot(aes(x=pval_AB)) +
  geom_histogram() +
  labs(title = "Four Group Linear Regression Analysis Early Changes P-values (Same Sample)",
       x = "P-value (Between Group A B, early changes")


test_BC_trans %>% 
  ggplot(aes(x=pval_BC)) +
  geom_histogram() 

test_CD_trans %>% 
  ggplot(aes(x=pval_CD)) +
  geom_histogram() 

inner_join(
  inner_join(test_AB_trans,
      test_BC_trans),
  test_CD_trans) %>% 
  rowwise %>% 
  mutate(min_pval = min(unlist(c(pval_AB, pval_BC, pval_CD)))) %>% 
  ggplot(aes(x= min_pval)) +
  geom_histogram() +
  labs(x = "Smallest p-value", title = "4 Group Analysis (Same 400 Genes) ")


# Number of significant changes for each comparasion
sig_AB<-test_AB_trans %>%
  drop_na() %>%
  filter(pval_AB<0.05) %>%
  select(transcript)


  
sig_BC<-test_BC_trans %>%
  drop_na() %>%
filter(pval_BC<0.05) %>%
  select(transcript)

sig_CD<-test_CD_trans %>%
  drop_na() %>%
filter(pval_CD<0.05) %>%
  select(transcript)

# Collection of all statistically significant differentially expressed genes
sig_genes<-rbind(sig_AB, 
                 sig_BC, 
                 sig_CD) %>% 
  distinct()

intersect(sig_genes %>% select(transcript) %>% rename(Gene_name = transcript), TABI_result)


nrow(sig_genes)




# Total number of statically signifcant genes
n_sig_genes<-nrow(sig_genes)

# Genes that change late / for those in high risk group
# Defined by changing between groups C / D but not between A / B or B / C
late_sig_genes<-anti_join(anti_join(sig_CD, 
                                    sig_BC), 
                          sig_AB) %>% 
  distinct()

n_late_sig_gene<-nrow(late_sig_genes)

# Gene with late Significant change and increasing (postive slope)
late_increase_sig_genes<-inner_join(test_CD_trans, 
                                    late_sig_genes) %>% 
  filter(slope_CD>0) %>% 
  select(transcript)

nrow(late_decrease_sig_genes)

cat(as.vector(late_increase_sig_genes$transcript)[501:920], sep="\n")

# Genes with late sig change and decrease (negative slope)
late_decrease_sig_genes<-inner_join(test_CD_trans, 
                                    late_sig_genes) %>% 
  filter(slope_CD<0) %>%   
  select(transcript)


#Same as above but for early gene changes 
# Defined by differential gene expression between A /B  but not C/D or B /C 
early_sig_genes<-anti_join(anti_join(sig_AB, 
                                     sig_BC), 
                           sig_CD) %>% 
  distinct()

early_increase_sig_genes<-inner_join(test_AB_trans, early_sig_genes) %>% 
  filter(slope_AB>0) %>% 
  select(transcript)

early_decrease_sig_genes<-inner_join(test_AB_trans, early_sig_genes) %>% 
  filter(slope_AB<0) %>% 
  select(transcript)


#Plot Calculated Slope vs Time (Region i.e. Early / Vs Late)

rbind(inner_join(test_CD_trans, late_sig_genes) %>% 
        mutate(Change = 1) %>%
        select(Change, slope_CD) %>%
        rename(slope = slope_CD), 
      inner_join(test_AB_trans, early_sig_genes) %>% 
        mutate(Change = 0) %>%
        select(Change, slope_AB) %>%
        rename(slope = slope_AB)) %>%
  ggplot(aes(x=Change, 
             y= slope)) + 
  geom_point(alpha = 0.5, 
             colour = "dodgerblue") +
  stat_smooth(method="lm", 
              se=F,
              col="red") +
  labs(title = "Calculated Slope vs Time of differential Gene Change",
       subtitle = "Signifcant Genes Only (alpha < 0.05)",
       x= "Region of Change (0 = Early, 1 = Late)",
       y= "Calculated value of Slope")

# Adding in middle changes 
mid_sig_genes<-anti_join(anti_join(sig_BC, sig_AB), 
                         sig_CD) %>% 
  distinct()

# Plot as above but with middle changes included
rbind(inner_join(test_CD_trans, late_sig_genes) %>% 
        mutate(Change = 2) %>%
        select(Change, slope_CD) %>%
        rename(slope = slope_CD), 
      inner_join(test_AB_trans, early_sig_genes) %>% 
        mutate(Change = 0) %>%
        select(Change, slope_AB) %>%
        rename(slope = slope_AB),
      inner_join(test_BC_trans, mid_sig_genes) %>%
       mutate(Change = 1) %>%
        select(Change, slope_BC) %>%
        rename(slope = slope_BC)) %>%
  ggplot(aes(x=Change, y= slope)) + 
  geom_point(alpha = 0.5, colour = "dodgerblue") +
  stat_smooth(method="lm", se=F, col="red") +
  labs(title = "Calculated Slope vs Time of differential Gene Change",
       subtitle = "Signifcant Genes Only (alpha < 0.05)",
       x= "Region of Change (0 = Early, 1 = Middle, 2 = Late)",
       y= "Calculated value of Slope")



# With grouping in three groups
# Grouping and then performing linear regression analysis 

# Define Group A as  CAPRA_S 0-2, B 3-5, C > 5
Group_test_table_3<-compress_normalised_narm_TCGA %>%
  mutate(CAPRA_S_group = ifelse(CAPRA_S>=0&CAPRA_S<=2, "A"
                                , ifelse(CAPRA_S>=3&CAPRA_S<=5,
                                         "B","C")))  %>%
  select(CAPRA_S_group, `read_count normalised`, purity.score, transcript)%>%
  rename(`read_count_normalised` = `read_count normalised`) %>%
  dplyr::mutate_if(is.character, as.factor) 

#Test differential gene expression between group A / B
test_AB_3_group<-Group_test_table_3 %>%
  filter(CAPRA_S_group=="A"|CAPRA_S_group=="B") %>%
  group_by(transcript) %>%
  summarise(p_val = (tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,5],
            CAPRA_S_slope = (tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,4])

test_AB_trans_3_group<-data.frame(pval_AB = t(test_AB_3_group$p_val), 
                                  slope_AB =t(test_AB_3_group$CAPRA_S_slope), 
                                  transcript = test_AB_3_group$transcript)

#Test differential gene expression between group B / C

test_BC_3_group<-Group_test_table_3 %>%
  filter(CAPRA_S_group=="B"|CAPRA_S_group=="C") %>%
  group_by(transcript) %>%
  summarise(CAPRA_S_p_val = (broom::tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,5],
            CAPRA_S_slope = (broom::tidy(summary(lm(read_count_normalised ~ CAPRA_S_group + purity.score))))[2,4])

test_BC_trans_3_group<-data.frame(pval_BC = t(test_BC_3_group$CAPRA_S_p_val), 
                            slope_BC =t(test_BC_3_group$CAPRA_S_slope), 
                            transcript = test_BC_3_group$transcript)

library(ggplot2)

#? Log

test_BC_3_group_log<-Group_test_table_3 %>%
  filter(CAPRA_S_group=="B"|CAPRA_S_group=="C") %>%
  group_by(transcript) %>%
  summarise(CAPRA_S_p_val = (broom::tidy(summary(lm(log(read_count_normalised+1) ~ CAPRA_S_group + purity.score))))[2,5],
            CAPRA_S_slope = (broom::tidy(summary(lm(log(read_count_normalised+1) ~ CAPRA_S_group + purity.score))))[2,4])

test_BC_trans_3_group_log<-data.frame(
  pval_BC = t(test_BC_3_group_log$CAPRA_S_p_val), 
                                  slope_BC =t(test_BC_3_group_log$CAPRA_S_slope), 
                                  transcript = test_BC_3_group_log$transcript)

test_BC_trans_3_group_log %>%
  ggplot() +
  geom_histogram(aes(x=pval_BC))


test_AB_trans_3_group %>%
  ggplot() +
  geom_histogram(aes(x=pval_AB))

test_AB_3_group_log<-Group_test_table_3 %>%
  filter(CAPRA_S_group=="A"|CAPRA_S_group=="B") %>%
  group_by(transcript) %>%
  summarise(p_val = (tidy(summary(lm(log(read_count_normalised+1) ~ CAPRA_S_group + purity.score))))[2,5],
            CAPRA_S_slope = (tidy(summary(lm(log(read_count_normalised+1) ~ CAPRA_S_group + purity.score))))[2,4])

test_AB_trans_3_group<-data.frame(
  pval_AB = t(test_AB_3_group_log$p_val), 
                                  slope_AB =t(test_AB_3_group_log$CAPRA_S_slope), 
                                  transcript = test_AB_3_group_log$transcript)



test_AB_trans_3_group %>%
  ggplot() +
  geom_histogram(aes(x=pval_AB))


test_BC_trans_3_group %>%
  ggplot() +
  geom_histogram(aes(x=pval_BC))

cbind(test_AB_trans_3_group, test_BC_trans_3_group %>% select(-transcript)) %>%
  
  ggplot() +
  geom_point(aes(x=pval_AB, y=pval_BC))


library(plotly)
ggplotly()

(test_AB_trans_3_group$pval_AB)



test_AB_trans_3_group

#Significant genes for each comparasion 
sig_AB_3_group<- test_AB_trans_3_group %>%
  drop_na() %>%
  filter(pval_AB <0.05) %>%
  select(transcript)

sig_BC_3_group<- test_BC_trans_3_group %>%
  drop_na() %>%
  filter(pval_BC<0.05) %>%
  select(transcript)


# Collection of all statistically significant differentially expressed genes
sig_3_group<-rbind(sig_AB_3_group, sig_BC_3_group) %>%
  distinct()

# Define early changes (only A/B, not in B/C)
early_sig_3_group<-anti_join(sig_AB_3_group, sig_BC_3_group) %>% 
  distinct()

#Increasing + Early
early_increase_sig_genes_group_3 <-inner_join(early_sig_3_group %>% 
             select(transcript),
             test_AB_trans_3_group %>% 
               filter(slope_AB>0))

# Decreasing + Early 
early_decrease_sig_genes_group_3<-inner_join(early_sig_3_group %>% 
                                               select(transcript),
                                             test_AB_trans_3_group %>% 
                                               filter(slope_AB<0))

# Define late changes  (only B/C not in A/B)

late_sig_3_group<-anti_join(sig_BC_3_group, 
                            sig_AB_3_group)

# Increasing + Late

late_increase_sig_genes_group_3 <-inner_join(late_sig_3_group %>% 
                                                select(transcript),
                                              test_BC_trans_3_group %>% 
                                               filter(slope_BC>0))

# Decreasing + Late

late_decrease_sig_genes_group_3 <-inner_join(late_sig_3_group %>% 
                                               select(transcript),
                                             test_BC_trans_3_group %>% 
                                               filter(slope_BC<0))


# Plot Direction of Change (Slope) againts Time of Change (Early vs Late)
nrow(inner_join(late_sig_3_group, 
                test_BC_trans_3_group) %>% 
       filter(slope_BC>0))

nrow(inner_join(test_BC_trans_3_group, 
                late_sig_3_group %>% 
                  select(transcript)) %>%
  mutate(Group =1) %>%
  select(Group, slope_BC) %>%
  rename(slope = slope_BC) )

rbind(inner_join(test_AB_trans_3_group, 
                 early_sig_3_group) %>%
        mutate(Group = 0) %>%
        select(Group, slope_AB) %>%
        rename(slope = slope_AB),
    inner_join(test_BC_trans_3_group, 
               late_sig_3_group) %>%
      mutate(Group =1) %>%
      select(Group, slope_BC) %>%
      rename(slope = slope_BC)) %>%
  ggplot(aes(x=Group, y=slope)) + 
geom_point(alpha = 0.5, colour = "dodgerblue") +
  stat_smooth(method="lm", se=F, col="red") +
  labs(title = "Calculated Slope vs Time of differential Gene Change",
       subtitle = "Signifcant Genes Only (alpha < 0.05)",
       x= "Region of Change (0 = Early, 1 = Late)",
       y= "Calculated value of Slope",
       caption = "Three groups (A (CAPRA-S 0-2), B (CAPRA-S 3-5), and C (CAPRA-S >5)")

#Comparing signifcant genes (method with 3 groups, method with 4 groups)
#Overlapping genes marked as significant
n_same<-nrow(dplyr::intersect(sig_3_group, 
                              sig_genes))

#Comparing genes early / late 
#Number of genes detected as being late changes in both analysis 
n_same_late<-nrow(dplyr::intersect(late_sig_3_group, 
                                   late_sig_genes))

#Number different
n_diff_late<-nrow(dplyr::setdiff(late_sig_3_group, 
                                 late_sig_genes))

#Percentage same between analysis (late genes)
n_same_late/(n_diff_late + n_same_late)

#Repeat as above but for early
n_same_early<-nrow(dplyr::intersect(early_sig_3_group, 
                                    early_sig_genes))

n_diff_early<-nrow(dplyr::setdiff(early_sig_3_group, 
                                  early_sig_genes))
#Percentage same
n_same_early/(n_diff_early + n_same_early)

# Comparing Genes discovered in above methods with those found by TABI
# Complete table of fit results (no outlier removal, purity included as covariate)

fit_table_purity_complete <- readRDS("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/TABI_results_and_analysis_purity/500_iter/fit_table_purity_complete.rds")

TABI_sig_results<-fit_table_purity_complete %>%
  filter(parameters =="beta[2,1]") %>%
  filter(`2.5%`<0&`97.5%`<0|`2.5%`>0&`97.5%`>0) %>% #Selecting either all > 0 for 95% interval / or all < 0
  select(Gene_name)


#All significant genes
#Overlap between TABI sig genes and 4 groups

sig_genes_TABI_4group_overlap<-dplyr::intersect(sig_genes, 
                                                TABI_sig_results %>% 
                                                  rename(transcript = Gene_name))

nrow(sig_genes_TABI_4group_overlap)

#Difference between TABI sig genes and 4 groups
sig_genes_TABI_4group_diff<-dplyr::setdiff(sig_genes, 
                                           TABI_sig_results %>% 
                                             rename(transcript = Gene_name))

nrow(sig_genes_TABI_4group_diff)

#Overlap between TABI sig genes and 3 groups 

sig_genes_TABI_3group_overlap<-dplyr::intersect(sig_3_group, 
                                                TABI_sig_results %>% 
                                                  rename(transcript = Gene_name))

nrow(sig_genes_TABI_3group_overlap)

#Difference between TABI sig genes and 3 groups 
sig_genes_TABI_3group_diff<-dplyr::setdiff(sig_3_group, 
                                           TABI_sig_results %>% 
                                             rename(transcript = Gene_name))

nrow(sig_genes_TABI_3group_diff)

# Overlap between all three methods 
sig_genes_all_methods<-dplyr::intersect(sig_genes_TABI_3group_overlap, 
                                        sig_genes_TABI_4group_overlap)

nrow(sig_genes_all_methods)

library(VennDiagram)
grid.newpage()
draw.triple.venn(
  area1 = nrow(TABI_sig_results), #Number of TABI marked sig genes
  area2 = nrow(sig_3_group), # Number of method with 3 groups marked sig genes
  area3 = nrow(sig_genes), # Number of method with 4 groups marked sig genes
  n12 = nrow(sig_genes_TABI_3group_overlap), # number of genes marked sig by 3 group method, and TABI
  n23 = n_same, #Number of overlaping sig marked genes (Method 3 group, Method 4 group)
  n13 = nrow(sig_genes_TABI_4group_overlap), #Number of genes marked sig by 4 group method and TABI
  n123 = nrow(sig_genes_all_methods), #Number of genes marked sig in all methods
  category = c("TABI", "3 Groups", "4 Groups"),
  scaled = T,
  fill = c("red", "dodgerblue", "yellow")
)

#Genes marked as early changes only
#Overlap with TABI and 4 group method
early_localised_TABI<-inner_join(fit_table_purity_complete %>%
             filter(parameters=="inflection[1]", mean<0),
           TABI_sig_results) %>% 
  select(Gene_name)

nrow(early_localised_TABI)

early_decrease_local_TABI<-inner_join(fit_table_purity_complete %>%
             filter(parameters=="beta[2,1]", mean<0), 
           early_localised_TABI) %>%
  select(Gene_name)

early_increase_local_TABI<-inner_join(fit_table_purity_complete %>%
                                        filter(parameters=="beta[2,1]", mean>0), 
                                      early_localised_TABI) %>%
  select(Gene_name)


n_early_overlap_TABI_4group<-nrow(dplyr::intersect(early_localised_TABI %>% 
                                                     rename(transcript = Gene_name), 
                                                   early_sig_genes %>% 
                                                     select(transcript)))

#Overlap with TABI and 3 group method
n_early_overlap_TABI_3group<-nrow(dplyr::intersect(early_localised_TABI %>% 
                                                    rename(transcript = Gene_name), 
                                                  early_sig_3_group %>% 
                                                    select(transcript)))

#Overlap between all three methods 
n_early_total_overlap<-nrow(dplyr::intersect(n_early_overlap_TABI_3group, 
                                             n_early_overlap_TABI_4group))


#Venn Diagram Options

#rea1 The size of the first set
#area2 The size of the second set
#area3 The size of the third set
#n12 The size of the intersection between the first and the second set
#n23 The size of the intersection between the second and the third set
#n13 The size of the intersection between the first and the third set
#n123 The size of the intersection between all three sets

#Overlap in early changes between methods

grid.newpage()
draw.triple.venn(
  area1 = nrow(early_localised_TABI), #Number of TABI marked sig genes + localised + early
  area2 = nrow(early_sig_3_group), # Number of method with 3 groups marked sig genes +early 
  area3 = nrow(early_sig_genes), # Number of method with 4 groups marked sig genes + early
  n12 = n_early_overlap_TABI_3group, # number of genes marked sig by 3 group method, and TABI + early
  n23 = n_same_early, #Number of overlaping sig marked genes (Method 3 group, Method 4 group) + early
  n13 = n_early_overlap_TABI_4group, #Number of genes marked sig by 4 group method and TABI
  n123 = 0, #Number of genes marked sig in all methods
  category = c("TABI", "3 Groups", "4 Groups"),
  scaled = T,
  fill = c("red", "dodgerblue", "yellow")
)

#Genes marked as late changes only
#Overlap with TABI and 4 group method
late_localised_TABI<-inner_join(fit_table_purity_complete %>%
                                   filter(parameters=="inflection[1]", mean>0),
                                 TABI_sig_results) %>% 
  select(Gene_name)


late_decrease_local_TABI<-inner_join(fit_table_purity_complete %>%
                                        filter(parameters=="beta[2,1]", mean<0), 
                                      late_localised_TABI) %>%
  select(Gene_name)

nrow(late_decrease_local_TABI)


late_increase_local_TABI<-inner_join(fit_table_purity_complete %>%
                                       filter(parameters=="beta[2,1]", mean>0), 
                                     late_localised_TABI) %>%
  select(Gene_name)

nrow(late_increase_local_TABI)

n_late_overlap_TABI_4group<-nrow(dplyr::intersect(late_localised_TABI %>% 
                                                     rename(transcript = Gene_name), 
                                                   late_sig_genes %>% 
                                                     select(transcript)))

#Overlap with TABI and 3 group method
n_late_overlap_TABI_3group<-nrow(dplyr::intersect(late_localised_TABI %>% 
                                                     rename(transcript = Gene_name), 
                                                   late_sig_3_group %>% 
                                                     select(transcript)))

#Overlap between all three methods 
n_late_total_overlap<-nrow(dplyr::intersect(n_late_overlap_TABI_3group, 
                                            n_late_overlap_TABI_4group))

#Overlap in late changes between methods

grid.newpage()
draw.triple.venn(
  area1 = nrow(late_localised_TABI), #Number of TABI marked sig genes + localised + early
  area2 = nrow(late_sig_3_group), # Number of method with 3 groups marked sig genes +early 
  area3 = nrow(late_sig_genes), # Number of method with 4 groups marked sig genes + early
  n12 = n_late_overlap_TABI_3group, # number of genes marked sig by 3 group method, and TABI + early
  n23 = n_same_late, #Number of overlaping sig marked genes (Method 3 group, Method 4 group) + early
  n13 = n_late_overlap_TABI_4group, #Number of genes marked sig by 4 group method and TABI
  n123 = 0, #Number of genes marked sig in all methods
  category = c("TABI", "3 Groups", "4 Groups"),
  scaled = T,
  fill = c("red", "dodgerblue", "yellow")
)

#




