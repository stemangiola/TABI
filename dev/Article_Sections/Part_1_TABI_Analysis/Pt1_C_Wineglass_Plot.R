
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Intention of Script

# Create Wineglass plot for TABI based on realistic simulated data 

# Figure 1 C 

# Requires sigmoidal_sim_df function (see Functions_Simulate_RNAseq script)





#####################################################################



# Section 1: Simulate realistic datasets for TABI
# Run TABI on these dataset 



# Section 2: Wineglass plot of the dataset (Figure_C)



# Section 3: Supplementary Figure of FDR estimate via Posterior Probability of Inclusion 
# (i.e. posterior prob. of differential abundance)



# Section 4: Wine glass plot for Bayseq for comparasion



#####################################################################


# Set up 

# If package box is not installed, then install package 
if (
  !("box" %in% installed.packages()[ , "Package"])) {
  install.packages("box")
}

# Load functions from Functions_Simulate_RNA_seq_df using box 
options(box.path = setwd("~/TABI/dev/Article_Sections"))   


#' @export
box::use(./Functions_Simulate_RNAseq_df[...]) 

# Load required packages 

box::use(tidyr[...], 
         furrr[...],
         ggplot2[...],
         data.table[...]) 


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



# Section 1: Create dataset of realistic data 
# Analyse using TABI 


# 1.1 Simulate data  (using sigmodial function from Functions_Simulated_Sigmoid_df_RNAseq.R
# Parameter + function simulate set up


# Range of x-coordinates to simulate over 


x_cord = seq(from = -5, 
             to = 5, 
             by = 0.5) 


# For differentially expressed genes
# Sigmoidal parameter values to simulate over


# K values to simulate over 

# Use TCGA dataset to create a distribution of realistic ks to simulate over

load("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Article_Sections/Part_1_TABI_Analysis/TCGA_Prostate_Simple.rda")


# Use TCGA data - and take the mean at low CAPRA-S score (0) and high CAPRA-S score (8) for each genes
# Then log of the difference is an empirical estimate of k

empirical_ks= TCGA %>% 
  filter(CAPRA_S == 0 | CAPRA_S>7) %>% 
  group_by(transcript, CAPRA_S) %>% 
  summarise(mean = mean(read_count_normalised, 
                        na.rm = TRUE)) %>% 
  mutate(CAPRA_S = paste0("CAPRA_S_", CAPRA_S)) %>% 
  tidyr::pivot_wider(values_from = mean, 
                     names_from = CAPRA_S) %>% 
  rowwise() %>% 
  mutate(dif_80 = CAPRA_S_8-CAPRA_S_0) %$%
  dif_80 %>% 
  log() %>% 
  abs()

# Remove inappropriate values (NAs and InF)
empirical_ks = empirical_ks[!empirical_ks %in% c(-Inf, NA, NaN, Inf)]

# Plot distribution of k values 
data.frame(x = empirical_ks) %>% 
  ggplot(aes(x=x)) + 
  geom_histogram(bins = 1000) + 
  labs(x = "Empirical k values (From Prostate TCGA dataset)",
       title = "Histogram of Empirical k values (From Prostate TCGA dataset)")

# Quantiles of k values
empirical_ks %>% 
  quantile(c(0.5, 0.75, 0.9))

# Distribution of k values to simulate over 
# Assume top only ~10% of k values are truly DE 
# Hence, the distribution of k values to simulate over 
# is the filtered top ~10% ks (i.e. 90% percentile)
sim_ks = 
  data.frame(sim_ks = empirical_ks) %>% 
  filter(sim_ks > empirical_ks %>% quantile(c(0.9))) 


# Plot distribution of ks in 90th percentile
  sim_ks %>% 
  ggplot(aes(x=sim_ks)) + 
  geom_histogram(bins = 100) + 
  labs(x = "Empirical k values (From Prostate TCGA dataset, 90% percentile)",
      title = "Histogram of Empirical k values (From Prostate TCGA dataset, 90% percentile)")
  
# Sample this resulting distribution of ks as 
#Set of k values to simulate over 
    
ks = sample(sim_ks$sim_ks, 
            size = n_true_genes, 
            replace = T)  


# A values to simulate over

# Filtering out genes which have ks in the bottom 90% 
# I.e. genes likely to be null 
true_genes = TCGA %>% 
  filter(CAPRA_S == 0 | CAPRA_S>7) %>% 
  group_by(transcript, CAPRA_S) %>% 
  summarise(mean = mean(read_count_normalised)) %>% 
  mutate(CAPRA_S = paste0("CAPRA_S_", CAPRA_S)) %>% 
  tidyr::pivot_wider(values_from = mean, 
                     names_from = CAPRA_S) %>% 
  rowwise() %>% 
  mutate(dif_80 = CAPRA_S_8-CAPRA_S_0) %>% 
  mutate(dif_80 = log(dif_80)) %>%
  mutate(dif_80 = abs(dif_80)) %>% 
  filter(dif_80>=empirical_ks %>% quantile(c(0.9))) %>% 
  select(transcript)

# For these true genes estimate k as 
#  log(mean at CAPRA-S 0 ) - log(mean at CAPRA-S 8)
empirical_ks_true=  inner_join(TCGA, true_genes) %>% 
  filter(CAPRA_S == 0 | CAPRA_S>7) %>% 
  group_by(transcript, CAPRA_S) %>% 
  summarise(mean = mean(read_count_normalised, na.rm = TRUE)) %>% 
  mutate(CAPRA_S = paste0("CAPRA_S_", CAPRA_S)) %>% 
  tidyr::pivot_wider(values_from = mean, 
                     names_from = CAPRA_S) %>% 
  rowwise() %>% 
  mutate(dif_80 = CAPRA_S_8-CAPRA_S_0) %$%
  dif_80 %>% 
  log() %>% 
  abs()

# Then k + A for these genes is the log of the mean of the upper plateau 
# estimate this as the maximum of the log mean at CAPRA-S 0 or at CAPRA-S 8
empirical_kA=inner_join(TCGA, true_genes) %>% 
  filter(CAPRA_S == 0|CAPRA_S == 8) %>% 
  group_by(transcript, 
           CAPRA_S) %>% 
  summarise(mean = mean(read_count_normalised, 
                        na.rm = TRUE)) %>% 
  mutate(CAPRA_S = paste0("CAPRA_S_", CAPRA_S)) %>% 
  tidyr::pivot_wider(values_from = mean, 
                     names_from = CAPRA_S) %>% 
  rowwise() %>% 
  mutate(max = max(CAPRA_S_8, 
                   CAPRA_S_0)) %>% 
  filter_if(~is.numeric(.), 
            all_vars(!is.infinite(.))) %$%
  max %>% 
  log() 

# A = k + A - k
empirical_A= abs(empirical_kA - empirical_ks_true)


empirical_A = empirical_A[!empirical_A %in% c(-Inf, NA, NaN, Inf)]


data.frame(sim_A = empirical_A) %>% 
  ggplot(aes(x=sim_A)) + 
  geom_histogram(bins = 100) + 
  labs(x = "Empirical A values (From Prostate TCGA dataset, 90% percentile of ks)",
       title = "Histogram of Empirical A values (From Prostate TCGA dataset, 90% percentile of ks)")


As = sample(empirical_A, 
            size = n_true_genes, 
            replace = T) #Set of A values to simulate


# sample sizes to simulate over

 = 1:135

# slopes to simulate over 

slopes = 

set.seed(106)

# Vector of parameter combinations to simulate over
par_vector= data.frame(slopes,
                        ks,
                        As,
                        sample_size_DE, 
                        disp_DE, 
                        inflection = "all", 
                        set_seed = FALSE)


library(future)


# n workers for parallel 
n_workers = 15

plan(multisession, 
     workers = n_workers)


# Apply function of to all possible parameter combinations (for DE genes)
wine_glass_data_DE = future_pmap_dfr(.l = par_vector, 
                .f =  ~sigmoidal_sim_df(
                  n_true_tests = length(x_cord)*3,
                  n_false_tests = 0,
                  beta = ..1,
                  k = ..2,
                  A = ..3,
                  sample_size = ..4,
                  disp_size = ..5,
                  x_cord =  seq(from = -5, to =5, by =0.5)),
                .progress = TRUE) 



# wine_glass_data_DE= mapply(FUN = function(b, k, a, s, d, i) # slopes, ks, As, dispersion, sample sizes, inflection 
#   sigmoidal_sim_df(
#     1,
#     0,
#     b, #beta 
#     k, #k values 
#     a, # A values 
#     s, # sample size
#     d, # dispersion value
#     x_cord, # x_cord range 
#     i*s), # alpha = inflection*slope
#   par_vector[,1], #slopes
#   par_vector[,2], #ks
#   par_vector[,3], #As
#   par_vector[,4], #sample sizes 
#   par_vector[,5], #dispersions
#   par_vector[,6],  # inflection values
#   SIMPLIFY = FALSE)

wine_glass_data_non_DE = 
  sigmoidal_sim_df(
    n_true_tests = 0,
    n_false_tests = n_false_genes,
    beta = 0,
    k = 0,
    A = 0,
    sample_size = sample_size_non_DE,
    disp_size = disp_non_DE
  )


# Combine simulated data into a single dataframe 
wine_glass_data = rbind(wine_glass_data_non_DE, 
                        wine_glass_data_DE) 



# Check dataset (plot random example genes)



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




# Section 2: Wineglass plot of the dataset 





#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



# Section 3: Create wineglass plot for Bayseq based on newly simulated data 

# Using above simulated data 






