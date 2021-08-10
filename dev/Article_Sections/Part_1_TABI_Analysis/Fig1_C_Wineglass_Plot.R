
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Intention of Script

# Create Wineglass plot for TABI based on realistic simulated data 

# Figure 1 C 

# Requires sigmoidal_sim_df function (see TABI_Article_Functions.R script)



#####################################################################



# Section 1: Simulate realistic datasets for TABI using parameters from TCGA data
# Combine datasets and then
# Run TABI on these dataset 



# Section 2: Figure C Wineglass plot of the dataset (Figure_C)
# which is the posterior probability of exclusion (i.e. the posterior prob of differential abundance), 
# vs. log fold gene change


# Section 3: 
# Supplementary Figure of FDR estimate via Posterior Probability of Inclusion 
# (i.e. posterior prob. of differential abundance)



###############################################################################


# Set up 

# If package box is not installed, then install package 
if (
  !("box" %in% installed.packages()[ , "Package"])) {
  install.packages("box")
}

# Load functions from Functions_Simulate_RNA_seq_df using box 
options(box.path = setwd("~/TABI/dev/Article_Sections"))   


#' @export
box::use(./TABI_Article_Functions[...]) 

# Load required packages 

box::use(tidyr[...], 
         furrr[...],
         ggplot2[...],
         data.table[...],
         ggpubr[...]) 


################################################################################


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


###############################################################################

# Section 1 : Simulate realistic datasets for TABI using parameters from TCGA data



# Obtain set of predominantly DE genes using edgeR 
# to create a set of reasonable parameters for sigmoidal function for DE genes 

# For simplicity sake, split samples in 'high' vs 'low' CAPRA-S scores and test for DE between these two groups

# Load TCGA data
load("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Copy_Of_Article_Sections/TCGA_Prostate_Simple.rda")


# Split samples into two groups + put sample column in format which is recognised by tidybulk
TCGA_two_group = TCGA %>% 
                rowwise() %>% 
                 mutate(group = ifelse(CAPRA_S < 4, "A", "B")) %>% 
                 mutate(sample= as.factor(sample)) 

# Use tidybulk to test for DE using edgeR
TCGA_two_group_DE = tidybulk::test_differential_abundance(
  .data = TCGA_two_group, 
  .formula = ~group + purity.score,  # include purity.score as a covariate
  .sample = sample, #sample column
  .transcript = transcript, #transcript id column
  .abundance = read_count_normalised, # read count column
  method = "edgeR_quasi_likelihood", 
  scaling_method = "none" #data is already normalised, no need to repeat that
)

# Filter top genes by FDR < 0.05, 
top_genes = TCGA_two_group_DE %>% 
  filter(FDR<0.05) %>% 
  group_by(transcript, CAPRA_S) %>% 
  filter(CAPRA_S < 1 | CAPRA_S >= 7) 

# group CAPRA_S 7/8 togheter, by replacing every instance of CAPRA_S = 8 with CAPRA_S = 7
top_genes[top_genes == 8] = 7
  

#calculate mean normalised read count at each 'plateau' 
# (for rough indicator, assume plateaus are at CAPRA_S = 0, or CAPRA_S = 7/8)
top_genes_mean = top_genes %>% 
  summarise(mean = mean(read_count_normalised)) %>% 
  ungroup() %>% 
  group_by(transcript) %>% 
  summarise(min = (min(mean)), 
            max = (max(mean))) %>% 
  rowwise() %>%
  mutate(dif = max - min)

ggarrange(
top_genes_mean %>% 
  ggplot(aes(x=min, y = dif)) + 
  geom_point(), 
top_genes_mean %>% 
  ggplot(aes(x=min, y = dif)) + 
  geom_point() + 
  ylim(0, 20000) + 
  xlim(0, 10000),
top_genes_mean %>% 
  ggplot(aes(x=min, y = dif)) + 
  geom_point() + 
  ylim(0, 2500) + 
  xlim(0, 2500))

  top_genes_mean %>% 
    ggplot(aes(x=log(min), y = log(dif))) + 
    geom_point()

sample_ak<-function(){
    
    repeat {
      
      a = runif(n = 1,
                min = -10,
                max = 10)
      
      k = runif(n = 1,
                min = 0,
                max = 10)
      
      if(a<k) break}
    
    return(data.frame(a = a,k = k))
    
  }  
  

ka =
  Reduce(
        "rbind",
        replicate(n = 1000, 
                  sample_ak(), 
                  simplify = FALSE)
                                    )




ka %>% 
  ggplot(aes(x=a, y=k)) + 
  geom_point()


  

As = c(
      rep(-50, 150),
      runif(n = 1350, 
           min = -10, 
           max  = 10))

As = runif(n = 1500, 
           min = -10, 
           max  = 10)

K_A = runif(n = 1500,
            min = -10,
            max = 10)

TABI_Article_Functions$plot_rsim_df(
TABI_Article_Functions$sigmoidal_sim_df(
  n_true_tests = 1, #Number of different alphas per element in slopes
  n_false_tests = 0, #No null tests
  beta = 1, 
  k = 7,
  A = 0,
  sample_size = 63,
  disp_size = 1, 
  covar =  seq(from = -5, to =5, by =0.5),
  alpha = c(0,0), 
  seed_n  = FALSE))

sim_A_k = data.frame(a = As, 
           k = K_A + As, 
           ka = K_A) 


ggarrange(
  sim_A_k %>% 
  ggplot(aes(x = As, y= ka)) + 
  geom_point(),
  sim_A_k %>% 
    ggplot(aes(x = As, y= ka)) + 
    geom_point() +
    xlim(0, 400) + 
    ylim(0, 400))


data.frame(a = As, 
           k = K_A + As, 
           ka = K_A) %>% 
  ggplot(aes(x = As, y= k)) + 
  geom_point()

summary(top_genes_mean$dif)

ggarrange(
  data.frame(x=log(top_genes_mean$max)) %>% 
  ggplot(aes(x=x)) + 
  geom_histogram() + 
  ylim(0, 200), 
data.frame(x= (rgamma(n = 1049, shape = 3, scale = 1.5))) %>% 
  ggplot(aes(x=x)) + 
  ylim(0, 200) + 
  geom_histogram()) 


ggarrange(
  data.frame(x=log(top_genes_mean$dif)) %>% 
    ggplot(aes(x=x)) + 
    geom_histogram(bins = 100) + 
    ylim(0, 200), 
  data.frame(x= (rgamma(n = 1049, shape = 1, scale = 1.5))) %>% 
    ggplot(aes(x=x)) + 
    ylim(0, 200) + 
    geom_histogram(bins = 100)) 


ggarrange(
  data.frame(x=log(top_genes_mean$min)) %>% 
    ggplot(aes(x=x)) + 
    geom_histogram(bins = 100) + 
    ylim(0, 100), 
  data.frame(x= (rgamma(n = 1049, shape = 1.25, scale = 1.5))) %>% 
    ggplot(aes(x=x)) + 
    ylim(0, 100) + 
    geom_histogram(bins = 100)) 


testing %>%
  ggplot(aes(x=min, y=dif)) + 
  geom_point()


################################################################################


# 1.2 Simulate data  (using sigmodial function from Functions_Simulated_Sigmoid_df_RNAseq.R
# Parameter + function simulate set up

# For each 'true' gene simulate 10 null genes

# Range of x-coordinates to simulate over 


x_cord = seq(from = -5, 
                        to = 5, 
                                 by = 0.5) 


n_true_genes = 1500


As = sample(empirical_A, 
            size = n_true_genes, 
            replace = T) #Set of A values to simulate


n_false_genes = 1350


# For differentially expressed genes
# Sigmoidal parameter values to simulate over


# K values to simulate over 

# Use TCGA dataset to create a distribution of realistic ks to simulate over

load("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Copy_Of_Article_Sections/TCGA_Prostate_Simple.rda")









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
  quantile(c(0.5, 
                  0.75, 
                      0.9))

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
    
ks = sample(empirical_ks,
            size = n_true_genes, 
            replace = T)  


# A values to simulate over



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



# Remove inappropriate values (NAs and InF)
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

sample_size = 1:135

# slopes to simulate over 

slopes = rt(n = n_true_genes, 
                          df = 3)



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






