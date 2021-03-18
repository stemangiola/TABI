

#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Intention of Script



## Test / plot the impact of changing parameters in the Sigmoidal Curve equation
## on TABI Inferences vs. other differential abundance methods 

# Requires sigmoidal_sim_df function (see Functions_Simulated_Sigmoid_df script)



# Section 1 
# Simulate dataset over range of sample sizes, 
# slopes, inflections, dispersion values, A (lower plateau),  k (difference between lower and upper plateau
# and run TABI on dataset

# Section 2
# Impact of slope - to true inference of slope
# 1.1 Parameter Values to simulate curves for testing impact of slope 
# 1.2 Is TABI able to accurately retrieve slope information from curves?


# Section 2 
# Impact of inflection - to true inference of inflection 
# 2.1 Parameter Values to simulate curves for testing impact of slope 
# 2.2 Is TABI able to accurately retrieve inflection information from curves?


# Section 3 
# Impact of dispersion value on retrieving slope / inflection 


# Section 4
# Impact of A (lower plateau) on retrieving slope / inflection 


# Section 5 
# Impact of k (difference between lower and upper plateau) on retrieving slope / inflection 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

# 
# par= expand.grid(setdiff(,
#                         c(11, 22))




hi = future_pmap_dfr(.l = par, 
                              .f =  ~sigmoidal_sim_df(
                                n_true_tests = 2,
                                n_false_tests = 0,
                                beta = .x, 
                                k = median(ks),
                                A = median(As),
                                sample_size = .y,
                                disp_size = 0.85, 
                                alpha = c(0,0),
                                x_cord =  seq(from = -5, to =5, by =0.5)),
                              .progress = TRUE) 








par_vector= expand.grid(c(0.25, -0.25, 0.5, -0.5, 1, -1, 3, -3),
                        c(11,22,55, 77))


setdiff(seq(from=-4, to =4, by =0.5), c(0))

As = As

ks = ks

sample_sizes = seq(from = 10, to = 80, by = 10)
library(furrr)

plan(multisession, workers = 15)

sim_df_slopes_sample_size = future_pmap_dfr(.l = par_vector, 
                              .f =  ~sigmoidal_sim_df(
                                n_true_tests = 5,
                                n_false_tests = 0,
                                beta = .x, 
                                k = 8,
                                A = 3,
                                sample_size = .y,
                                disp_size = 0.85, 
                                alpha = "all",
                                x_cord =  seq(from = -2.5, to = 2.5, by = 0.5)),
                              .progress = TRUE) 


by_Gene_ref = sim_df_slopes_sample_size %>% 
  rename(sample = sample_id) %>% 
  split(.$Gene_ref)

n<- 15

# Set up parallelisation clusters
plan(multisession, 
     workers = n)


TABI_sample_size_slope = by_Gene_ref %>% 
  future_map(
    ~TABI::TABI_glm(
      .data = .,
      ~CAPRA_S,
      .abundance = value, # simulated RNA seq value
      .sample = sample, #sample column is sample_id
      .transcript = Gene_ref, # column of transcript / gene ids 
      chains = 4,
      cores = 1,
      iter = 2000, # total number of iterations
      warmup = 1000
    ) %$% 
      fit, .progress = TRUE) 




par_vector= expand.grid(slopes,
                        sample_sizes)




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 




# Section 1: Test Impact of Slope
# By creating simulated sigmoidal data with all other parameters held fixed
# Running TABI on data 
# Comparing inferences to true values

library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(tibble)




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




### Simulation dataset for testing slopes 

## Slopes from -5 to 5, holding all other factors constant


# Parameter + function simulate set up

#Set of slopes (beta values) to simulate
slopes = setdiff(
                  seq(from=-4, to =4, by =0.5), 
                                                c(0)) 


# Set all other values to a single default 

# To choose default consider TCGA data

load("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Article_Sections/Part_1_TABI_Analysis/TCGA_Prostate_Simple.rda")


#k value to simulate (take genes with >90% percentile of estimated ks as truly differentially expressed)
# then find the median distribution of k for these truly differentially expressed genes
ks_true_boundary =TCGA %>% 
  filter(CAPRA_S == 0 | CAPRA_S>7) %>% 
  group_by(transcript, CAPRA_S) %>% 
  summarise(mean = mean(read_count_normalised)) %>% 
  mutate(CAPRA_S = paste0("CAPRA_S_", CAPRA_S)) %>% 
  tidyr::pivot_wider(values_from = mean, names_from = CAPRA_S) %>% 
  rowwise() %>% 
  mutate(dif_80 = abs(CAPRA_S_8-CAPRA_S_0)) %$%
  dif_80 %>% 
  log() %>% 
  quantile(c(0.9)) %>% 
  round(digits = 1)



# Filtering out genes which have ks in the bottom 90% 
# I.e. genes likely to be null 
true_genes = TCGA %>% 
  filter(CAPRA_S == 0 | CAPRA_S>7) %>% 
  group_by(transcript, CAPRA_S) %>% 
  summarise(mean = mean(read_count_normalised)) %>% 
  mutate(CAPRA_S = paste0("CAPRA_S_", CAPRA_S)) %>% 
  tidyr::pivot_wider(values_from = mean, names_from = CAPRA_S) %>% 
  rowwise() %>% 
  mutate(dif_80 = abs(CAPRA_S_8-CAPRA_S_0)) %>% 
  filter(dif_80>=ks_true_boundary) %>% 
  select(transcript)

# For these true genes estimate k as 
#  log(mean at CAPRA-S 0 ) - log(mean at CAPRA-S 8)
empirical_ks_true= inner_join(TCGA, true_genes) %>% 
  filter(CAPRA_S == 0 | CAPRA_S>7) %>% 
  group_by(transcript, CAPRA_S) %>% 
  summarise(mean = mean(read_count_normalised, na.rm = TRUE)) %>% 
  mutate(CAPRA_S = paste0("CAPRA_S_", CAPRA_S)) %>% 
  tidyr::pivot_wider(values_from = mean, names_from = CAPRA_S) %>% 
  rowwise() %>% 
  mutate(dif_80 = abs(CAPRA_S_8-CAPRA_S_0)) %$%
  dif_80 %>% 
  log()

# Take the value of k to be the median of this distribution of ks 
ks = empirical_ks_true %>% 
  median()


# Then to find the value of A to fix for simulation 

# Estimate k+A for the proxy distribution of selected DE genes
# Then k + A for these genes is the log of the mean of the upper plateau 
# estimate this as the maximum of the log mean at CAPRA-S 0 or at CAPRA-S 8

empirical_kA = inner_join(TCGA, true_genes) %>% 
  filter(CAPRA_S == 0|CAPRA_S == 8) %>% 
  group_by(transcript, CAPRA_S) %>% 
  summarise(mean = mean(read_count_normalised, na.rm = TRUE)) %>% 
  mutate(CAPRA_S = paste0("CAPRA_S_", CAPRA_S)) %>% 
  tidyr::pivot_wider(values_from = mean, names_from = CAPRA_S) %>% 
  rowwise() %>% 
  mutate(max = max(CAPRA_S_8, CAPRA_S_0)) %>% 
  filter_if(~is.numeric(.), all_vars(!is.infinite(.))) %$%
  max %>% 
  log() 


empirical_A= (empirical_kA - empirical_ks_true)

# Value of A to fix in the equation
As = 

# Find the dispersion 

# Dispersion value to simulate
# as per Hardcastle & Kelly 2010 
# simulate from gamma distribution (shape = 0.85, scale = 0.5 )
disp = rgamma(shape = 0.85, scale = 0.5, n = 20) #Dispersion values 



x_cord = seq(from = -5, to = 5, by = 0.5) # Range of x-coordinates to simulate over


sample_sizes = length(x_cord)*4 #Sample size per gene to simulate - 4 samples per x_coordinate

alpha= 0 # Set of alpha to simulate - only simulate alpha at zero (which implies inflection is always zero)
n_true_tests = 5 # Number of sigmoidal curves to simulate - 5 for each slope
n_false_tests  = 0 # Number of null tests to simulate - none 



# Simulation dataset 

slope_sim<-
  lapply(slopes, function(s)
  sigmoidal_sim_df(n_true_tests = n_true_tests, 
         n_false_tests = n_false_tests, 
         beta = s, #Value of beta for equation
         k =ks , #k value ( k + a is upper plateau)
         A = As, #A value (lower plateau level)
         sample_size = sample_sizes, #Number of samples per gene - sample size needs to be a multiple of x_cord length
         disp_size = disp, #value of dispersion - used to simulate for negbinomial distribution 
         x_cord = x_cord, #vector of x-coordinates (simulated range pf CAPRA_S values)
         alpha = alpha, # vector of inflection values to simulate
         null_mean_distribution = TCGA_mean_distribution
)) 

sample_sizes

sample_size_sim<-
  lapply(c(21*1:4), function(s)
    sigmoidal_sim_df(n_true_tests = n_true_tests, 
                     n_false_tests = n_false_tests, 
                     beta = 2, #Value of beta for equation
                     k =ks , #k value ( k + a is upper plateau)
                     A = As, #A value (lower plateau level)
                     sample_size = s, #Number of samples per gene - sample size needs to be a multiple of x_cord length
                     disp_size = disp, #value of dispersion - used to simulate for negbinomial distribution 
                     x_cord = x_cord, #vector of x-coordinates (simulated range pf CAPRA_S values)
                     alpha = alpha, # vector of inflection values to simulate
                     null_mean_distribution = TCGA_mean_distribution
    )) 

slope_sim<- Reduce(rbind, slope_sim) 


sample_size_sim<- Reduce(rbind, sample_size_sim) 




# Structure of dataset 
slope_sim %>% 
  head()

slope_sim %>% 
  dim()

slope_sim %>% 
  select(Gene_ref) %>% 
  distinct() %>% 
  nrow()


# Plot random simulated gene as a check 
plot_rsim_df(slope_sim)

plot_rsim_df(sim_df_slopes_sample_size)
  

# Run TABI on the data 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
  

#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


saveRDS(slope_sim, "slope_sim.rds")


saveRDS(sample_size_sim, "sample_size_sim.rds")
  
# 1.2   
# Is TABI able to accurately retrieve slope information from curves?
  
  
## Run TABI on this data frame

library(furrr)

## number of clusters to run in parallel

n<- 15

# Set up parallelisation clusters
plan(multisession, 
         workers = n)


# Run TABI on each gene, in parallel, 
# Return the summary statistics on the fit of the parameters 
by_Gene_ref= sim_df_slope_sample_size %>% 
  rename(sample = sample_id) %>% 
  split(.$Gene_ref)

by_Gene_ref_sample_size = sample_size_sim %>% 
  rename(sample = sample_id) %>% 
  split(.$Gene_ref)


by_Gene_ref = sim_data_DE %>% 
  rename(sample = sample_id) %>% 
  split(.$Gene_ref)

n<- 15

# Set up parallelisation clusters
plan(multisession, 
     workers = n)


TABI_sample_size_slope = by_Gene_ref %>% 
  future_map(
    ~TABI::TABI_glm(
      .data = .,
      ~CAPRA_S,
      .abundance = value, # simulated RNA seq value
      .sample = sample, #sample column is sample_id
      .transcript = Gene_ref, # column of transcript / gene ids 
      chains = 4,
      cores = 1,
      iter = 2000, # total number of iterations
      warmup = 1000
    ) %$% 
      fit) 

saveRDS(TABI_sample_size_slope, file = "TABI_sample_size_slope.rds")

by_Gene_ref_slope = slope_sim %>% 
    rename(sample = sample_id) %>% 
    split(.$Gene_ref)


Slope_TABI = by_Gene_ref_slope %>% 
  future_map(
    ~TABI::TABI_glm(
      .data = .,
      ~CAPRA_S,
      .abundance = value, # simulated RNA seq value
      .sample = sample, #sample column is sample_id
      .transcript = Gene_ref, # column of transcript / gene ids 
      chains = 4,
      cores = 1,
      iter = 2000, # total number of iterations
      warmup = 1000
    ) %$% 
      fit) 
saveRDS(Slope_TABI, "Slope_TABI.rds")





slope_TABI_fit= (slope_sim %>% 
  rename(sample = sample_id) %>% 
  split(.$Gene_ref)) %>% 
future_map(~TABI::TABI_glm(
  .data = .x,
  ~CAPRA_S,
  .abundance = value, # simulated RNA seq value
  .sample = sample, #sample column is sample_id
  .transcript = Gene_ref, # column of transcript / gene ids 
  chains = 4,
  cores = 1,
  iter = 2300, # total number of iterations
  warmup = 1000
)  %$% 
  fit) %>% 
  map_dfr(~ rstan::summary(.x) %$% 
            summary %>% 
            as.data.frame())
toc()

saveRDS(slope_TABI_fit, "slope_TABI_fit.rds")

# Extract data on the fit of the slope only

slope_fit = slope_TABI_fit %>% 
  tibble::rownames_to_column(var = "Pars") %>% 
  filter(grepl("beta", Pars))

# Get the order of Gene IDs that was parallelised over 

Gene_ref= (slope_sim %>% 
             rename(sample = sample_id) %>% 
             split(.$Gene_ref)) %>% 
  names()

# Combine TABIs fit with their respective Gene IDs
slope_inference = data.frame(slope_fit, Gene_ref) 

# Extract data on the underlying slope of each simulated gene / transcript

slope_true = slope_sim %>% 
  select(Gene_ref, slope) %>% 
  distinct()


# Combine data on TABI's inference on slope
# and true data on the slope 

slope_inference = inner_join(slope_true, 
                             slope_inference)

# Plot true slope data vs inferred slope by TABI
  slope_inference %>% 
  rename(upper = 'X97.5.', lower = 'X2.5.') %>% 
  ggplot(aes(x=slope, y= mean)) + 
  geom_point(col = "dodgerblue") + 
  geom_smooth(method = lm, se = FALSE, col = "red", lty = 2) + 
    geom_errorbar(aes(x= slope,
                      ymin = lower,
                      ymax = upper),
                  col = 'black',
                  alpha = 0.2) 
  

  
  
  
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  
  
  
  
  
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



  
  
  
  
  #   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
  
  
  
  
# Section 2: Test Impact of Inflection 

# Inflection across all possible x-coordinate values, holding all other factors constant  

slopes = 1 #Set of slopes (beta values) to simulate
ks = ks  #Set of k values to simulate
As = As #Set of A values to simulate
disp =  #Dispersion values 
x_cord = seq(from = -5, to = 5, by = 0.5) # Range of x-coordinates to simulate over
sample_sizes = length(x_cord)*1:5 #Sample size per gene to simulate - 4 samples per x_coordinate
alpha = "all" # #Range of alphas to simulate
n_true_tests = length(x_cord)*5 # Number of sigmoidal curves to simulate - 5 for each inflection
n_false_tests  = 0 # Number of null tests to simulate - none 


inflection_sim<-sigmoidal_sim_df(n_true_tests = n_true_tests, 
                         n_false_tests = n_false_tests, 
                         beta = slopes, #Value of beta for equation
                         k =ks , #k value ( k + a is upper plateau)
                         A = As, #A value (lower plateau level)
                         sample_size = sample_sizes, #Number of samples per gene - sample size needs to be a multiple of x_cord length
                         disp_size = disp, #value of dispersion - used to simulate for negbinomial distribution 
                         x_cord = x_cord, #vector of x-coordinates (simulated range pf CAPRA_S values)
                         alpha = alpha # vector of inflection values to simulate
)



# Structure of dataset 
inflection_sim %>% 
  head()

inflection_sim %>% 
  dim()

inflection_sim %>% 
  select(Gene_ref) %>% 
  distinct() %>% 
  nrow()


# Plot random simulated gene as a check 
plot_rsim_df(inflection_sim)



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




# 2.2   
# Is TABI able to accurately retrieve inflection information from curves?


## number of clusters to run in parallel

n<-12

# Set up parallelisation clusters
plan(multisession, workers = n)



library(tictoc)
tic()
# Run TABI on each gene, in parallel, 
# Return the summary statistics on the fit of the parameters 
inflection_TABI_fit= (inflection_sim %>% 
                   rename(sample = sample_id) %>% 
                   split(.$Gene_ref)) %>% 
  future_map(~TABI::TABI_glm(
    .data = .x,
    ~CAPRA_S,
    .abundance = value, # simulated RNA seq value
    .sample = sample, #sample column is sample_id
    .transcript = Gene_ref, # column of transcript / gene ids 
    chains = 4,
    cores = 1,
    iter = 1300, # total number of iterations
    warmup = 300
  )  %$% 
    fit) %>% 
  map_dfr(~ rstan::summary(.x) %$% 
            summary %>% 
            as.data.frame())
toc()


saveRDS(inflection_TABI_fit, "inflection_TABI_fit.rds")

# Extract data on the fit of the slope only

inflection_fit = inflection_TABI_fit %>% 
  tibble::rownames_to_column(var = "Pars") %>% 
  filter(grepl("beta", Pars))

# Get the order of Gene IDs that was parallelised over 

Gene_ref= (inflection_sim %>% 
             rename(sample = sample_id) %>% 
             split(.$Gene_ref)) %>% 
  names()

# Combine TABIs fit with their respective Gene IDs
inflection_inference = data.frame(inflection_fit, Gene_ref) 

# Extract data on the underlying slope of each simulated gene / transcript

inflection_true = inflection_sim %>% 
  select(Gene_ref, inflec) %>% 
  distinct()


# Combine data on TABI's inference on slope
# and true data on the slope 

inflection_inference = inner_join(inflection_true, 
                             inflection_inference)

# Plot true slope data vs inferred slope by TABI
inflection_inference %>% 
  rename(upper = 'X97.5.', lower = 'X2.5.') %>% 
  mutate(upper = upper/10) %>% 
  mutate(lower = lower/10) %>% 
  ggplot(aes(x=inflec, y= mean)) + 
  geom_point(col = "dodgerblue") + 
  geom_smooth(method = lm, se = FALSE, col = "red", lty = 2) + 
  geom_errorbar(aes(x= inflec,
                    ymin = lower,
                    ymax = upper),
                col = 'black',
                alpha = 0.2) 







##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



# Section 3 Test: Impact of dispersion value on inferring slope / inflection 


# Parameter + function simulate set up
slopes = 2 #Set of slopes (beta values) to simulate
ks = ks  #Set of k values to simulate
As = As #Set of A values to simulate
set.seed(102)
disp = rgamma(shape = 0.85, scale = 0.5, n = 20) #Dispersion values 
x_cord = seq(from = -5, to = 5, by = 0.5) # Range of x-coordinates to simulate over
sample_sizes = length(x_cord)*4 #Sample size per gene to simulate - 4 samples per x_coordinate
alpha = 1 # Set of alphas to simulate
n_true_tests = length(disp)*5 # Number of sigmoidal curves to simulate - 5 for each dispersion value
n_false_tests  = 0 # Number of null tests to simulate - none 


disp_slopes_sim<-
  lapply(disp, function(d)
  sigmoidal_sim_df(n_true_tests = n_true_tests, 
                        n_false_tests = n_false_tests, 
                        beta = slopes, #Value of beta for equation
                        k =ks , #k value ( k + a is upper plateau)
                        A = As, #A value (lower plateau level)
                        sample_size = sample_sizes, #Number of samples per gene - sample size needs to be a multiple of x_cord length
                        disp_size = d, #value of dispersion - used to simulate for negbinomial distribution 
                        x_cord = x_cord, #vector of x-coordinates (simulated range pf CAPRA_S values)
                        alpha = alpha # vector of inflection values to simulate
)) 




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 






#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


