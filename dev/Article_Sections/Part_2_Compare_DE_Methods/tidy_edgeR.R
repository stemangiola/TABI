

ROC_table<-rbind(FDR_sim_21$FDR_ROC_edgeR %>% 
                   mutate(Test_type = paste("edgeR", group)) %>% 
                   select(TPR, FPR, Test_type),
                 FDR_sim_21$FDR_ROC_DEseq %>% 
                   mutate(Test_type = paste("DESeq2", group)) %>% 
                   select(TPR, FPR, Test_type),
                 FDR_sim_21$FDR_ROC_Bayseq_df[[1]] %>% 
                   select(TPR, FPR) %>% 
                   mutate(Test_type = "Bayseq_2_group"),
                 FDR_sim_21$FDR_ROC_Bayseq_df[[2]]  %>%
                   select(TPR, FPR) %>% 
                   mutate(Test_type = "Bayseq_3_group"),
                 FDR_sim_21$FDR_ROC_Bayseq_df[[3]] %>% 
                   select(TPR, FPR) %>% 
                   mutate(Test_type = "Bayseq_4_group")
)  

# In general edgeR analysis


# EdgeR analysis for factors of 21 

edgeR_analysis<-function(data, n_group){
  
  require(tidyverse)
  require(edgeR)
  require(foreach)
  require(rlang)
  
  assign_df<- function(name , df) {
    assign(rlang::quo_name(dplyr::enquo(name)), 
           df, 
           envir=.GlobalEnv)
  }
  
  
  groups<-c(paste("edgeR_df_", letters, sep=""))[1:n_group]
  
  CAPRA_S_val<-seq(from = -5, to =5, by = 0.5)
  
  if((length(CAPRA_S_val) %% n_group) ==0){

  n_CAPRA_per_group<-length(CAPRA_S_val)/n_group
  
  
  edgeR_dfs<-foreach(i = c(1:n_group))%do%{ 
    if(i==1){
      assign_df(cat(groups[i]), data %>% 
                  filter(CAPRA_S<CAPRA_S_val[n_CAPRA_per_group]) %>% 
                  select(Gene_ref, value, sample) %>% 
                  pivot_wider(names_from = sample, values_from = value) %>% 
                  tibble::column_to_rownames(var= "Gene_ref")
      )
    }
    
  
    
    else{
      assign_df(cat(groups[i]), data %>%
                  filter(CAPRA_S>=CAPRA_S_val[n_CAPRA_per_group*(i-1)]& CAPRA_S<CAPRA_S_val[n_CAPRA_per_group*(i)]) %>%
                  select(Gene_ref, value, sample) %>% 
                  pivot_wider(names_from = sample, values_from = value)%>%
                  tibble::column_to_rownames(var= "Gene_ref")
      ) 
    }
  }
  
  } 
  

  else{
  
    (length(CAPRA_S_val) %% n_group)==j
    
    n_CAPRA_per_group<-((length(CAPRA_S_val) - j) %% n_group)
    

  edgeR_dfs<-foreach(i = c(1:n_group))%do%{ 
    if(i==1){
      assign_df(cat(groups[i]), data %>% 
        filter(CAPRA_S<=CAPRA_S_val[n_CAPRA_per_group + 1]) %>% 
        select(Gene_ref, value, sample) %>% 
        pivot_wider(names_from = sample, values_from = value) %>% 
        tibble::column_to_rownames(var= "Gene_ref")
        )
    }
    
  
    
    else if(i<=j){
      assign_df(cat(groups[i]), data %>% 
                  filter(CAPRA_S_val[(n_CAPRA_per_group + 1)*(i-1)]<CAPRA_S&CAPRA_S<=CAPRA_S_val[(n_CAPRA_per_group + 1)*i]) %>% 
                  select(Gene_ref, 
                         value, 
                         sample) %>% 
                  pivot_wider(names_from = sample, 
                              values_from = value) %>% 
                  tibble::column_to_rownames(var= "Gene_ref")
      )
    }
    
  

    else{
   assign_df(cat(groups[i]), data %>%
      filter(CAPRA_S>CAPRA_S_val[(n_CAPRA_per_group + 1)*(i-1) + n_CAPRA_per_group*(j)]&CAPRA_S<CAPRA_S_val[(n_CAPRA_per_group + 1)*(i+1)+ n_CAPRA_per_group*(j)]) %>%
      select(Gene_ref, value, sample) %>% 
      pivot_wider(names_from = sample, values_from = value)%>%
      tibble::column_to_rownames(var= "Gene_ref")
      ) 
    }
  }
  
  }

  
  edgeR_groups<-bind_cols(edgeR_dfs) %>% as.matrix()

  groups<- list()
  for(i in c(1:n_group)){
  
  assign(
    paste0("n_",letters[i]),
    ncol(edgeR_dfs[i])
  )
    
  groups[[i]]<-rep(LETTERS[i], ncol(edgeR_dfs[[i]]))  
  
  } 
  
  groups<-Reduce(c,groups)
 
  
  y <- DGEList(counts = edgeR_groups, group = groups)
  
  #y<- dgList(edgeR_groups)
  
  #y <- calcNormFactors(y)
  
  design <- model.matrix(~groups)
  
  #my.contrasts<-makeContrasts()
  
  y <- estimateDisp(y, 
                      design)
  
  fit<- glmQLFit(y, 
                    design)
  
  #glmQLFTest
  
  res_tables<-foreach(i=c(1:(n_group-1)))%do% {
  
    contrast_val<-rep(0, n_group)
    
    contrast_val[c(i, i+1)]<-c(-1,1)
    
    
    res<-glmQLFTest(fit, 
                  contrast = contrast_val) 
    
    pval<- res %$%
      table %>%
      rownames_to_column(var = "Gene_ref") %>% 
      select(Gene_ref, PValue) %>% 
     dplyr::rename(!!paste0("Pval_", LETTERS[i], LETTERS[i+1]) := PValue)
    
    fold_change<-res %$%
      table %>%
      rownames_to_column(var = "Gene_ref") %>% 
      select(Gene_ref, logFC) %>% 
      dplyr::rename(!!paste0("logFC_", LETTERS[i], LETTERS[i+1]) := logFC)
               
  return(inner_join(pval, fold_change)) 
  }       

  if((length(CAPRA_S_val) %% n_group) ==0) {
    mult= length(CAPRA_S_val) %/% n_group

    jump_below<-CAPRA_S_val[1:(n_group-1)*mult]

    jump_up<-CAPRA_S_val[1:(n_group-1)*mult+1]
    
    let_vec<-c()

    for(i in 1:(n_group-1)){
      let_vec[i]<-paste0(LETTERS[i], LETTERS[i+1])}

    when_data<-data.frame(Area = let_vec,
                          jump_up = jump_up,
                          jump_below = jump_below) %>%
      mutate(point_est = rowMeans(select(., jump_up, jump_below)))

  }
  
  else{
    (length(CAPRA_S_val) %% n_group)==j
    
    mult = (length(CAPRA_S_val) - j) %/% n_group
    
    jump_up<-CAPRA_S_val[c(1:j*(mult+1), length(CAPRA_S_val) - 0:(n_group - j) *mult)]
    
    jump_below<-CAPRA_S_val[c(1:j*(mult+1)-1, length(CAPRA_S_val) - 0:(n_group - j) *mult-1)]

    let_vec<-c()
    for(i in 1:(n_group - 1)){
      let_vec[i]<-paste0(LETTERS[i], LETTERS[i+1])}
    
    when_data<-data.frame(Area = let_vec,
                          jump_up = jump_up,
                          jump_below = jump_below) %>%
      mutate(jump_up  = as.numeric(jump_up)) %>% 
      mutate(jump_down = as.numeric(jump_down)) %>% 
      mutate(point_est = rowMeans(select(., jump_up, jump_below)))
    
  }
  

  test = Reduce(inner_join, res_tables)

      test_data<-left_join(
    Reduce(inner_join, res_tables) %>%
             select(Gene_ref, contains("logFC")) %>%
               mutate_if(is.numeric, abs) %>%
               gather(which, value, -Gene_ref) %>%
               group_by(Gene_ref) %>%
               filter(rank(-value) == 1) %>%
      separate(which, c("logFC", "Area")),
    when_data)
  
  return(test)  
  
}

test<-Sig_multiple_value_data

edgeR_df_A_3<-test %>% 
  filter(CAPRA_S <2*-1) %>% 
  select(Gene_ref, value, Sample) %>% 
  pivot_wider(names_from = Sample, values_from = value) %>% 
  column_to_rownames(var= "Gene_number")

edgeR_df_B_3<-test %>% 
  filter(CAPRA_S >=(2*-1)&CAPRA_S<(2)) %>% 
  select(Gene_number, value, Sample) %>% 
  pivot_wider(names_from = Sample, values_from = value) %>% 
  column_to_rownames(var= "Gene_number")

edgeR_df_C_3<-test %>% 
  filter(CAPRA_S >=2) %>% 
  select(Gene_number, value, Sample) %>% 
  pivot_wider(names_from = Sample, values_from = value) %>% 
  column_to_rownames(var= "Gene_number")

edgeR_df_3<-cbind(
  edgeR_df_A_3,
  edgeR_df_B_3,
  edgeR_df_C_3
) 


n_A_3<-ncol(edgeR_df_A_3)
n_B_3<-ncol(edgeR_df_B_3)
n_C_3<-ncol(edgeR_df_C_3)

group_3<-c(rep("A", n_A_3), 
           rep("B", n_B_3),
           rep("C", n_C_3))


y_3 <- DGEList(edgeR_df_3)

y_3 <- calcNormFactors(y_3)

design_3 <- model.matrix(~group_3)


y_3 <- estimateDisp(y_3, 
                    design_3)

fit_3 <- glmQLFit(y_3, 
                  design_3)

lrt_early_3 <- glmQLFTest(fit_3, 
                          coef = 2)

lrt_late_3 <- glmQLFTest(fit_3, 
                         coef = 3)


Pval_3_group<-(inner_join(lrt_late_3$table %>% 
                            select(PValue) %>% 
                            rownames_to_column(var = "Gene_number") %>% 
                            dplyr::rename(Pval_late = PValue),
                          lrt_early_3$table %>% select(PValue) %>% 
                            rownames_to_column(var = "Gene_number") %>% 
                            dplyr::rename(Pval_early = PValue)))

edgeR_analysis(Sig_multiple_value_data %>% rename(sample = Sample), 4) %>% 
      select(contains("logFC")) %>% 
  abs() %>% 
  rownames_to_column('id') %>%
  gather(which, value, -id) %>% 
  group_by(id) %>% 
  filter(rank(-value) == 1)  %>% 
  separate(which, c("logFC", "Area")) %>% 
  
  if((length(CAPRA_S_val) %% n_group) ==0) {
    mult= length(CAPRA_S_val) %/% n_group
    
    jump_up<-CAPRA_S_val[1:mult*n_group]
    
    jump_below<-CAPRA_S_val[1:mult*n_group-1]
    let_vec<-c()
    
    for(i=1:mult){
    let_vec[i]<-paste0(LETTERS[i], LETTERS[i+1])} 
    
    when_data<-data.frame(Area = let_vec, 
               jump_up = jump_up, 
               jump_below = jump_below)
    
  }

inner_join(when_data,
select(contains("logFC")) %>% 
  abs() %>% 
  rownames_to_column('id') %>%
  gather(which, value, -id) %>% 
  group_by(id) %>% 
  filter(rank(-value) == 1)) 





make_groups<-function(table, 
                      n_group){
  
  CAPRA_S_val<-seq(from = -5, to =5, by = 0.5)
  
  if((length(CAPRA_S_val) %% n_group) ==0){
    
    n_per_group<-length(CAPRA_S_val)/n_group
    
    groups<-rep(n_per_group, n_group)
    
    CAPRA_numbers<-c(1:(n_group-1)*n_per_group)
    
    contrast_values<-(CAPRA_S_val[CAPRA_numbers]+CAPRA_S_val[CAPRA_numbers-1])/2
    
    
    n_per_CAPRA_S<-(table %>% select(Sample) %>% max())/ length(CAPRA_S_val)
    
    group_let<-list()
    
    for(i in 1:length(groups)){
      
      group_let[[i]]<-rep(LETTERS[i],n_per_CAPRA_S*groups[i]) 
      
    }
    
    } 
  
  else{
    n_per_group<-length(CAPRA_S_val) %/% n_group 
    
    n_left_over<-length(CAPRA_S_val) %% n_group
    
    groups<-c(rep(n_per_group+1,n_left_over), rep(n_per_group, n_group))
    
    CAPRA_numbers<-0
    
    n_per_CAPRA_S<-(table %>% select(Sample) %>% max())/ length(CAPRA_S_val)
    group_let<-list()
    
    for(i in 1:n_group){
    
      if(i<=n_left_over){
        CAPRA_numbers[i]<-(n_per_group+1)*i
        
        group_let[[i]]<-rep(LETTERS[i],(n_per_group+1)*n_per_CAPRA_S) 
        }
      
    else{
      CAPRA_numbers[i]<-CAPRA_numbers[i-1]+n_per_group
      
      group_let[[i]]<-rep(LETTERS[i],(n_per_group)*n_per_CAPRA_S) 
      }
    
    }
    CAPRA_numbers<-CAPRA_numbers[-n_group]
    
    contrast_values<-(CAPRA_S_val[CAPRA_numbers]+CAPRA_S_val[CAPRA_numbers-1])/2
} 
  
 
vector<-Reduce(c,group_let) 

return(list(group = vector,
            contrast_values = contrast_values)) 
 
}



edgeR_test<-function(table, n_group) {

tidy<-data.frame(table, 
                        group = make_groups(table, 
                                            n_group)$group
                 ) %>% 
  as_tibble()  %>% 
  mutate(Sample = as.character(Sample)) %>% 
  mutate(alpha = as.numeric(alpha)) %>% 
  mutate(alpha =  alpha*(-1)) %>% 
  tidybulk::tidybulk(.sample = Sample, 
                     .transcript = Gene_number, 
                     .abundance = value)


contrasts<-sapply(2:n_group, 
                  function(x) paste0("group", LETTERS[x], "-", "group", LETTERS[x-1]))

contrasts_names<-sapply(2:n_group, 
                        function(x) paste0(LETTERS[x-1], LETTERS[x]))

jump<-make_groups(table, 
            n_group)$contrast_values

con_df<-data.frame(contrast = contrasts,
           group = contrasts_names,
           jump_est = jump)



edgeR_analysis<-tidybulk::test_differential_abundance(
  tidy,
  ~0 +group,
  .contrasts = contrasts,
  method = "edgeR_quasi_likelihood",
  scaling_method = "none",
  omit_contrast_in_colnames = FALSE,
  #prefix = "",
  action = "get",
  significance_threshold = 1,
  fill_missing_values = FALSE
)

# 
# coeff<-attr(edgeR_analysis, "internals")$edgeR
# 
# coeff<-coeff$coefficients %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column("Gene_number") %>% 
#   as_tibble() %>% 
#   mutate("AB" = abs(groupB- groupA)) %>% 
#   mutate("BC" =  abs(groupB - groupC)) %>% 
#   mutate("CD" = abs(groupC - groupD)) %>% 
#   select(Gene_number, AB, BC, CD) %>% 
#   pivot_longer(!Gene_number, 
#                names_to= "group", 
#                values_to = "FC") %>%  
#   group_by(Gene_number) %>%
#   filter(rank(-FC) == 1)



logFC<-edgeR_analysis %>% 
  select(Gene_number, #alpha, 
         contains("FC"))  %>% 
  mutate_if(is.numeric, abs) %>%
  pivot_longer(!Gene_number, 
               names_to= "which", 
               values_to = "FC") %>%
  group_by(Gene_number) %>%
  filter(rank(-FC) == 1) %>% 
  separate(which, 
           c("name", "contrast"), 
           sep = "_")


alpha<- table %>% 
  select(Gene_number, True_Slope, alpha) %>% 
  distinct() %>% 
  mutate(alpha = ifelse(True_Slope>0,
                        -1*as.numeric(alpha), 
                        as.numeric(alpha))) 

  df<-left_join(alpha,
  left_join(logFC,
            con_df)) 
  
  #df1<-left_join(left_join(coeff, con_df), alpha)
  
return(df) 
}
Sig_multiple_value_data<-Sig_multiple_value_data %>% 
  filter(grepl("X", Gene_number))

test<-rbind(Slope_neg1 %>% mutate(Gene_number = paste0(Gene_number, "neg")), Slope_1)

grid.arrange(
  edgeR_test(test, 3) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "3 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") + 
    geom_point(), 
  edgeR_test(test, 4) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "4 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_jitter(size = 1),
  edgeR_test(test, 5) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "5 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_jitter(size = 1),
  edgeR_test(test, 6) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "6 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_jitter(size = 1),
  edgeR_test(test, 7) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "7 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_jitter(size = 1),
  edgeR_test(test, 8) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "8 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_jitter(size = 1),
  edgeR_test(test, 9) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "9 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_jitter(size = 1), 
  edgeR_test(test, 10) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "10 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_jitter(size = 1), 
  edgeR_test(test, 15) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "15 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_jitter(size = 1) 
)

grid.arrange(
  edgeR_test(Slope_1, 8) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "Slope 1, 15 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_point(), 
  edgeR_test(Slope_2, 8) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "Slope 2,8 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_point(),
  edgeR_test(Slope_05, 8) %>%
    ggplot(aes(x=alpha, y=jump_est)) + 
    labs(title = "Slope 0.5, 8 Groups", 
         x = "True Inflection",
         y= "Estimated Inflection") +
    geom_point() 
)


  
  