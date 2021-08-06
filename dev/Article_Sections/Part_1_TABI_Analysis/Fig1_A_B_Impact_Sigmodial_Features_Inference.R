

# Set up

# n cores / for running in parallel 
n =  availableCores()

#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


#Intention of Script



## Show TABI is able to retrieve slope and inflection for realistic simulated data

# Requires sigmoidal_sim_df function (see Functions_Simulate_RNAseq script)



################################################################################



# Section 1 : Slope / Data for Figure 1A
# Simulate data set over range of slopes
# run TABI on data set



# Section 2: Inflection / Data for Figure 1B
# Simulate data set over range of inflection values
# run TABI on data set


# Section 1: Plot Data for Figure 1A and Figure 1B 
# show TABI is able to retrieve realistic inflection values


# Section 3: Supplementary 
# show the relationship between inflection and correctly identifying DE transcripts

################################################################################

# Load required packages 
# box::use(dplyr[...], 
#          roxygen2[...],
#          furrr[...],
#          future[...],
#          tictoc[...],
#          purrr[...],
#          TABI[...]) 

library(dplyr)
library(roxygen2)
library(furrr)
library(TABI)
library(future)
library(purrr)
library(future)
library(rstan)
library(tictoc)
library(magrittr)
library(broom)
library(tidyr)

# If package box is not installed, then install package 
# if (
#      !("box" %in% installed.packages()[ , "Package"])) {
#   install.packages("box")
# }



# Load functions from Functions_Simulate_RNA_seq_df using box 
options(box.path = "/data/gpfs/projects/punim0586/ibeasley/TABI/dev/Article_Sections")   


#' @export
box::use(./TABI_Article_Functions) 

# If you want details and example usage 
# of the functions used with module TABI_Article_Functions
# use box::help(TABI_Article_Functions$function_name)

# e.g. 

#box::help(TABI_Article_Functions$sigmoidal_sim_df)

################################################################################


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


################################################################################

# Figure 1A

# Section 1: TABI retrieval of slope and slope impact on DE call 


# Change slope, but keep all other factors the same. 


#  #  #  #  #  #  #  


# Set up Parameter Values

# Slope / beta values to simulate over

# from running TABI on TCGA data, estimated |slope| was < 4, so simulate mostly within,
# but also slightly beyond this range i.e. |slope| <= 5

tic("Set up")


slopes = runif(n = 500, min = -5, max = 5) #Exclude slope of = 0 (this is not a DE gene!)



# Number of different slopes simulated

slopes %>%
  length() # 100

# slopes = rep(slopes, 10)


# All other parameters remain constant with:
# (reasonable parameters given TCGA data, see Figure 1C)

A = 6 #lower plateau

k = 6 + 7 #upper plateau


# n workers (for furrr parallelisation)

# cl <- parallel::makeCluster(n/4)
# plan(cluster, workers = cl)

plan(multisession, workers = n)

# Three samples / biological replicates per 'CAPRA-S' / covariate value

sample_size = seq(from = -5,
                  to =5,
                  by =0.5) %>% length()*3


# Simulate data sets of DE transcripts over above slope values
slopes_sim_df = furrr::future_map_dfr(.x = slopes,
                                      .f =  ~TABI_Article_Functions$sigmoidal_sim_df(
                                        n_true_tests = 1, #Number of different alphas per element in slopes
                                        n_false_tests = 0, #No null tests
                                        beta = .x,
                                        k = k,
                                        A = A,
                                        sample_size = sample_size,
                                        disp_size = 0.425,
                                        covar =  seq(from = -5, to =5, by =0.5),
                                        alpha = c(0,0),
                                        seed_n  = 151), # set seed
                                      .progress = TRUE)

# Check: Structure of data frame

slopes_sim_df %>%
  head()

slopes_sim_df %>%
  dim()

# Total number of simulated genes - should be 10*number of elements in slopes vector
slopes_sim_df %>%
  dplyr::select(Gene_ref) %>%
  distinct() %>%
  nrow()

# Check: plot some genes for realism

TABI_Article_Functions$plot_rsim_df(slopes_sim_df)


saveRDS(slopes_sim_df,
        file = "./slopes_sim_df.rds")

toc()

# Run this data set on TABI and extract fit data

plan(multisession, workers = n)

tic("Slope dataset TABI")

slopes_df_by_gene_ref = slopes_sim_df %>%
  dplyr::rename(sample = sample_id) %>%
  split(.$Gene_ref)


TABI_fit_slopes_list = slopes_df_by_gene_ref %>%
  furrr::future_map(
    ~TABI::TABI_glm(
      .data = .x,
      ~CAPRA_S,
      .abundance = value, # simulated RNA seq value
      .sample = sample, #sample column is sample_id
      .transcript = Gene_ref, # column of transcript / gene ids
      chains = 4,
      cores = 1,
      iter = 5000, # total number of iterations
      warmup = 2500
    ) %$%
      fit %>%
      rstan::summary() %$%
      summary %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "par") %>%
      cbind(., .x %>% dplyr::select(Gene_ref, slope, inflect, seed_number) %>% distinct()))


TABI_fit_slopes_df = Reduce("rbind",
                            TABI_fit_slopes_list)

data.table::fwrite(TABI_fit_slopes_df,
                   "./TABI_fit_slopes_df.csv")

toc()


##################################################################################


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


####################################################################################

# Figure 1B
# Section 2: TABI retrieval of inflection, and impact of inflection on DE calls


# # # # # # # Set up Parameters for simulation
tic("Inflection set up")

# Inflection  = - alpha / beta, so for simplicity in this example let beta = 1

slopes = 1 

# alpha values to simulate over / -inflection values to simulate over
# 1 value either side of the range of covariates

alpha = runif(n = 600, min = -6, max = 6)

# alpha = rep(alpha, 10)

# Number of different alpha simulated 

alpha %>% 
  length()


# All other parameters remain constant with 

k  = 6 + 7 # upper plateau 

A = 6 # lower plateau 


sample_size = seq(from = -5,
                  to =5,
                  by =0.5) %>% length()*3

# Set up furrr parrallelisation


# Simulate data sets of DE transcripts over alpha values

inflection_sim_df = furrr::future_map_dfr(.x = alpha,  
                                          .f =  ~TABI_Article_Functions$sigmoidal_sim_df(
                                            n_true_tests = 1, #Number of tests per alpha value
                                            n_false_tests = 0, #No null tests
                                            beta = slopes, 
                                            k = k,
                                            A = A,
                                            sample_size =  sample_size,
                                            disp_size = 0.425, 
                                            alpha = c(.x,.x), #Keep infection constant - and centered in the middle
                                            covar =  seq(from = -5, to =5, by =0.5),
                                            seed_n = 151),
                                          .progress = TRUE,
                                          .options = furrr_options(seed = TRUE)) 


# Check: Structure of data frame 
inflection_sim_df %>% 
  head()

inflection_sim_df %>% 
  dim()

# Total number of simulated genes - should be 5*number of different slopes
inflection_sim_df %>% 
  select(Gene_ref) %>% 
  distinct() %>% 
  nrow()

toc()

# Check: plot some genes for realism

#plot_rsim_df(slopes_sim_df)


saveRDS(inflection_sim_df, 
        file = "./inflection_sim_df.rds")


# Calculate the formula to converted unscaled inflection values to scaled inflection values

scale_xcord<- function(formula, data){
  
  require(tidyr)
  require(dplyr)
  
  scale_design = function(df, formula){
    df %>%
      setNames(c("sample_idx", "(Intercept)", parse_formula(formula))) %>%
      gather(cov, value, -sample_idx) %>%
      group_by(cov) %>%
      mutate( value = ifelse( !grepl("Intercept", cov) & length(union(c(0,1), value)) != 2, scale(value), value )) %>%
      ungroup() %>%
      pivot_wider(names_from = cov, values_from = value) %>%
      arrange(as.integer(sample_idx)) %>%
      select(`(Intercept)`, one_of(parse_formula(formula)))
  }
  
  parse_formula <- function(fm) {
    if(attr(terms(fm), "response") == 1) stop("The formula must be of the kind \"~ covariates\" ")
    else as.character(attr(terms(fm), "variables"))[-1]
  }
  
  X = model.matrix(object = formula, data = data) %>%
    as_tibble(rownames="sample_idx") %>%
    scale_design(formula)
  
  scaled_df = cbind(
    X %>% 
      select(-'(Intercept)') %>% 
      rename(scaled_CAPRA_S = CAPRA_S),
    data %>% 
      select(CAPRA_S)) 
  
  return(scaled_df)
  
}


inflection_df_by_gene_ref = inflection_sim_df %>% 
  rename(sample = sample_id) %>% 
  split(.$Gene_ref)

lm_scaled = scale_xcord(~CAPRA_S, 
                        inflection_df_by_gene_ref[[1]]) %>% 
  lm(scaled_CAPRA_S~CAPRA_S, data = .) 


scaled_inflect =tidy(lm_scaled)$estimate[1]+tidy(lm_scaled)$estimate[2]*(inflection_sim_df %>% 
                                                                           select(inflect) %>% 
                                                                           distinct() %$% inflect)

inflection_scaled = data.frame(inflect = inflection_sim_df %>% 
                                 select(inflect) %>% 
                                 distinct(), 
                               scaled_inflect = scaled_inflect)

inflection_sim_df = inner_join(inflection_sim_df, inflection_scaled)

print("All set up complete")

# Run this data set on TABI and extract fit data 

inflection_df_by_gene_ref = inflection_sim_df %>% 
  rename(sample = sample_id) %>% 
  split(.$Gene_ref)

tic("TABI Inflection")



TABI_fit_inflection_list = inflection_df_by_gene_ref %>% 
  furrr::future_map(
    ~TABI::TABI_glm(
      .data = .x,
      ~CAPRA_S,
      .abundance = value, # simulated RNA seq value
      .sample = sample, #sample column is sample_id
      .transcript = Gene_ref, # column of transcript / gene ids 
      chains = 4,
      cores = 1,
      iter = 5000, # total number of iterations
      warmup = 2000
    ) %$% 
      fit  %>% 
      rstan::summary() %$% 
      summary %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "par") %>%
      cbind(., .x %>% dplyr::select(Gene_ref, slope, inflect, scaled_inflect, seed_number) %>% distinct()),
    .options = furrr_options(seed = TRUE)) 

toc()


TABI_fit_inflection_df = Reduce("rbind",
                                TABI_fit_inflection_list)


data.table::fwrite(TABI_fit_inflection_df, 
                   "./TABI_fit_inflection_df.csv")



################################################################################


# Plotting Figure 1A + B

# Stop running in parallel

future:::ClusterRegistry("stop")

custom_theme <-
  list(
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 9),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  )


# Figure 1 A (Inferred vs True Inflection PLot)

Fig_1A = TABI_fit_inflection_df %>%
  dplyr::filter(par == "inflection[1]") %>%
  dplyr::rename(lower := "2.5%", upper := "97.5%") %>%
  dplyr::select(Gene_ref, mean, lower, upper,inflect, scaled_inflect) %>%
  ggplot(aes(x=scaled_inflect,
             y= mean)) +
  geom_point(col = "dodgerblue") +
  geom_abline(slope = 1, intercept = 0, col = "red", lty = 2) +
  geom_errorbar(aes(x= scaled_inflect,
                    ymin = lower,
                    ymax = upper),
                col = 'black',
                alpha = 0.2) +
  custom_theme +
  labs(x = "True Inflection",
       y = "TABI Estimated Inflection",
       title = "TABI Inferred Inflection vs Simulated Inflection Value")


Fig_1B = TABI_fit_slopes_df %>%
  dplyr::filter(par == "beta[1,1]") %>%
  dplyr::rename(lower := "2.5%", upper := "97.5%") %>%
  #dplyr::mutate(across(.cols = c(lower, upper, mean), ~.x/3.102418)) %>% 
  dplyr::mutate(slope = slope*3.102418) %>% 
  dplyr::select(Gene_ref, mean, lower, upper, slope) %>%
  rowwise() %>%
  #mutate(mean = log(mean), lower = log(lower), upper = log(upper)) %>%
  #filter(slope > 0) %>%
  ggplot(aes(x=slope,
             y= mean)) +
  geom_point(col = "dodgerblue") +
  geom_abline(slope = 1, intercept = 0, col = "red", lty = 2) +
  geom_errorbar(aes(x= slope,
                    ymin = lower,
                    ymax = upper),
                col = 'black',
                alpha = 0.2) +
  custom_theme +
  labs(x = "True Slope",
       y = "TABI Estimated Slope",
       title = "TABI Inferred vs Simulated Slope Value")

################################################################################
