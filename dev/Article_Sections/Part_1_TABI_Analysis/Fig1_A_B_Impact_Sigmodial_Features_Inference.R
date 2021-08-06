

# n cores / for running in parallel (set up as forked, needs to be changed if running on windows)
n_cores = 8

#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Intention of Script



## Show TABI is able to retrieve slope and inflection for realistic simulated data

# Requires sigmoidal_sim_df function (see Functions_Simulate_RNAseq script)



################################################################################



# Section 1 : Slope / Date for Figure 1A
# Simulate data set over range of slopes
# run TABI on data set
# show TABI is able to retrieve realistic slope values
# show the relationship between slope and correctly identifying DE transcripts



# Section 2: Inflection / Date for Figure 1B
# Simulate data set over range of inflection values
# run TABI on data set
# show TABI is able to retrieve realistic inflection values
# show the relationship between inflection and correctly identifying DE transcripts



################################################################################

# Load required packages 
box::use(dplyr[...], 
         roxygen2[...],
         furrr[...],
         future[...]) 


# If package box is not installed, then install package 
if (
     !("box" %in% installed.packages()[ , "Package"])) {
  install.packages("box")
}



# Load functions from Functions_Simulate_RNA_seq_df using box 
options(box.path = setwd("~/TABI/dev/Article_Sections"))   


#' @export
box::use(./TABI_Article_Functions) 

# If you want details and example usage 
# of the functions used with module TABI_Article_Functions
# use box::help(TABI_Article_Functions$function_name)

# e.g. 

box::help(TABI_Article_Functions$sigmoidal_sim_df)

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

slopes = setdiff( 
                   seq(
                       from=-5, 
                               to =5, 
                                    by =0.1), 
                                             c(0)) #Exclude slope of = 0 (this is not a DE gene!)



# Number of different slopes simulated 

slopes %>% 
      length() # 100

slopes = rep(slopes, 10)


# All other parameters remain constant with: 
# (reasonable parameters given TCGA data, see Figure 1C)

A = 6 #lower plateau

k = 6 + 7 #upper plateau


# n workers (for furrr parallelisation)

plan(multisession, 
     workers = 15)

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
                               seed_n  = FALSE),
                               .progress = TRUE) 

# Check: Structure of data frame 

slopes_sim_df %>% 
  head()

slopes_sim_df %>% 
  dim()

# Total number of simulated genes - should be 10*number of elements in slopes vector
slopes_sim_df %>% 
  select(Gene_ref) %>% 
  distinct() %>% 
  nrow()

# Check: plot some genes for realism

TABI_Article_Functions$plot_rsim_df(slopes_sim_df)


saveRDS(slopes_sim_df, 
        file = "~/TABI/dev/Article_Sections/Part_1_TABI_Analysis/slopes_sim_df.rds")


# Run this data set on TABI and extract fit data 

slopes_df_by_gene_ref = slopes_sim_df %>% 
                        rename(sample = sample_id) %>% 
                        split(.$Gene_ref) 
  
  
  TABI_fit_slopes = slopes_df_by_gene_ref %>% 
  future_map(
    ~TABI::TABI_glm(
      .data = .,
      ~CAPRA_S,
      .abundance = value, # simulated RNA seq value
      .sample = sample, #sample column is sample_id
      .transcript = Gene_ref, # column of transcript / gene ids 
      chains = 4,
      cores = 1,
      iter = 5000, # total number of iterations
      warmup = 2000
    ) %$% 
      fit) 
  

saveRDS(TABI_fit_slopes_df, 
        file = "~/TABI/dev/Article_Sections/Part_1_TABI_Analysis/TABI_fit_slopes_df.rds")
  
  summary %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "par") %>% 
  dplyr::filter(par == "inflection[1]") %>% 
  rename(upper = '97.5%', lower = '2.5%')
  
  
  # Plot true slope data vs inferred slope by TABI
  ggplot(aes(x=inflec, 
             y= mean)) + 
  geom_point(col = "dodgerblue") + 
  geom_smooth(method = lm, se = FALSE, col = "red", lty = 2) + 
  geom_errorbar(aes(x= inflec,
                    ymin = lower,
                    ymax = upper),
                col = 'black',
                alpha = 0.2) 





##################################################################################


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


####################################################################################

# Figure 1B
# Section 2: TABI retrieval of inflection, and impact of inflection on DE calls


# # # # # # # Set up Parameters for simulation


# Inflection  = - alpha / beta, so for simplicity in this example let beta = 1

slopes = 1 

# alpha values to simulate over / -inflection values to simulate over
# 1 value either side of the range of covariates

alpha = seq(from = -6,
                      to = 6,
                              by = 0.1)

alpha = rep(alpha, 10)

# Number of different alpha simulated 

alpha %>% 
          length()


# All other parameters remain constant with 

k  = 6 + 7 # upper plateau 

A = 6 # lower plateau 


# Set up furrr parrallelisation

plan(multicore, 
               workers = n)


# Simulate data sets of DE transcripts over alpha values

inflection_sim_df = future_map_dfr(.x = alpha,  
                               .f =  ~sigmoidal_sim_df(
                                 n_true_tests = 1, #Number of tests per alpha value
                                 n_false_tests = 0, #No null tests
                                 beta = slope, 
                                 k = k,
                                 A = A,
                                 sample_size =  sample_size,
                                 disp_size = 0.425, 
                                 alpha = c(.x,.x), #Keep infection constant - and centered in the middle
                                 covar =  seq(from = -5, to =5, by =0.5)),
                               .progress = TRUE) 

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

# Check: plot some genes for realism

plot_rsim_df(slopes_sim_df)


saveRDS(inflection_sim_df, 
        file = "/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Article_Sections/Part_1_TABI_Analysis/inflection_sim_df.rds")


# Run this data set on TABI and extract fit data 

inflection_df_by_gene_ref = inflection_sim_df %>% 
  rename(sample = sample_id) %>% 
  split(.$Gene_ref)
  
  
  TABI_fit_inflection_df = inflection_df_by_gene_ref %>% 
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


saveRDS(TABI_fit_inflection_df, 
        "/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Article_Sections/Part_1_TABI_Analysis/TABI_fit_inflection_df.rds")


################################################################################


