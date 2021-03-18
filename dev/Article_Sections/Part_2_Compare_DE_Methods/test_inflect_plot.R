

TABI_2000_iter_Fig_2_sim_data[[1]] %>% 
  filter(grepl("-1", Gene_ref) & grepl("X", Gene_ref)) %>% 
  select(Gene_ref) %>% 
  distinct()


TABI_2000_iter_Fig_2_sim_data[[1]] %>% 
  filter(grepl("_1_", Gene_ref) & grepl("X", Gene_ref)) %>% 
  select(Gene_ref) %>% 
  distinct()


test<-bind_rows(Fig_2_sim_data_sample_105 %>% 
  filter(grepl("_1_", Gene_ref) & grepl("X", Gene_ref)),
  Fig_2_sim_data_sample_105 %>% 
    filter(True_Slope>=0.9 & grepl("X", Gene_ref)))

ggplot(aes(x=N_groups, y= abs(alpha - jump_est), group = N_groups, fill = as.factor(N_groups))) + 
  +   geom_boxplot() + labs(x= "Number of Groups", y= "Error (Absolute Difference between True and Estimated Inflection)") + geom_boxplot(aes(x=0, y= dif/3.35, group = 0),fill = "dodgerblue",  data = TABI) + scale_fill_discrete(name = "Number of edgeR Groups")
FDR_ROC_Inflect_table(test %>% rename(Sample = sample), 3)



test = sigmodial_sim(n_true_tests = 100, 
              n_false_tests = 0, 
              beta = 1, 
               k = 3, 
               A = 5, 
               sample_size = 84, 
              disp_size = 4)

make_groups(test %>% rename(Sample = sample_id, True_Slope = slope), 3)

FDR_ROC_Inflect_table(test %>% rename(Sample = sample_id, True_Slope = slope), 3)$Inflect_est

