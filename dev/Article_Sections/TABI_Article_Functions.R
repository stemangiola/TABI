#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Intention of Script 

# Set up  /define functions used in all parts of analysis 


#(Section 1)
# Function to simulate a dataset - sigmoidal_sim_df
# of a mixture of sigmoidal based distributions (true associations)
# and null associations 

#(Section 2)
# Function to plot random gene from simulated data frame


#(Section 3) make_groups function
# Function which is able to make a data frame to match CAPRA_S values to group labels
# Given a desired number of groups to sort the simulated samples into 


#(Section 4)
# Function which performs differential expression by statifying samples into
# groups using make_groups function and tidybulk
# Can be edgeR, DESEq2, limma + voom

# (Section 5)
# Function which performs differential expression using Bayseq




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Section 1 


#Create simulation function

# For example simulated data frames 
# And sigmoid functions see script for Pt 1_1 



##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 


#Function for simulating Data with both null and true (sigmoidal curve) tests


#Function returns table with simulated data from
#specified number of true positive curves - 
#of a single equation but varying inflection value 
#And specified number of null curves


# Function is designed to simulate DE genes with a single set of A, k, beta, sample size, and dispersion size
# over a set of alpha values (or potentially a single alpha value)
# To simulate of other parameter values - need to map function over a set of parameter values

#' Simulate Dataset of Negative Bionomial Distributed Transcripts
#' 
#' @description 
#' \code{sigmoidal_sim_df} simulates a dataset (returned as a tibble) with both null and truly 
#' differentially abundant transcripts, with respect to a single continous or pseudo-continous variable.
#' The log mean at each covariate value (\code{log_y_hat}) for 
#' differentially abundant transcripts is set from a Richard's curve. 
#' 
#' @section Details: 
#' 
#' Function is designed to simulate a data set of a mixture of null (not differentially abundant)
#' and true (differentially abundant) transcripts. All transcripts are assumed to be negative bionomial distributed (parametised as mean and overdispersion).
#' The mean for null transcripts is constant across the covariate and 
#' is simulated from a vector of values (each
#' value in the vector a possible mean). By default, this vector is of transcript-wise
#' means from the TCGA prostate cancer dataset. 
#' 
#' Differentially abundant transcripts have a log mean which is dependent on the covariate. For a single data set returned by a single use of this function, set a 
#' single set of A, k, beta, sample size, and dispersion size. 
#' Differentially abundant transcripts are simulated over the full range of the covariate (default), or a given range of alpha values (vector of two values) or a single alpha value (if alpha is set to a single number). 
#' 
#' If you need to simulate a dataset with differentially abundant transcripts over a range of
#' Richard curve parameter values, map function over a set of parameter values 
#' keeping all other function parameters the same 
#' (recommended if combining these datasets afterwards: make sure to set \code{set_seed} to a single value)
#' 
#' If you only want to simulate null (not differentially abundant transcripts) then set \code{n_true_tests = 0}. 
#' 
#' If you only want to simulate differentially abundant transcripts, set \code{n_false_tests = 0}. 
#' 
#' 
#' 
#' 
#' @section Detailed Description of Returned Tibble: 
#' 
#' The returned data frame is in the tidy format. 
#' That is, each row is an observation, and each column is a variable.
#' 
#' 
#' \strong{Columns:}  
#' \tabular{lll}{
#' \code{sample_id} \tab \code{<int>} \tab Reference ID for the biological replicate. Can consider this to be the 'Sample' or 'Patient' ID column.  \cr
#' \code{Gene_ref} \tab  \code{<chr>} \tab Reference ID for the transcript. Can consider this to be the unique name of the gene / transcript.  \cr
#' \code{Gene_number} \tab \code{<chr>} \tab Reference ID for the transcript. If you are simulating multiple data sets, this value is not unique. (Instead use Gene_Ref Column) \cr
#' \code{Null_test} \tab \code{<log>} \tab Is this a null (TRUE) or differentially abundant transcript (FALSE)? \cr
#' \code{CAPRA_S} \tab \code{<dbl>} \tab Value of the covariate for this sample / patient. \cr
#' \code{value} \tab \code{<dbl>} \tab Read count / Transcript abundance for this observation. \cr
#' \code{log_y_hat}  \tab \code{<dbl>} \tab True log mean transcript abundance for this value of the covariate. \cr
#' \code{slope} \tab \code{<dbl>} \tab Value for the slope. If 'NA' this is not a differentially abundant transcript. \cr
#' \code{sample_size} \tab \code{<int>} \tab Number of biological replicates (samples) per transcript (single value per dataset) \cr
#' \code{inflect} \tab \code{<dbl>} \tab Value for the Inflection. If 'NA' this is not a differentially abundant transcript. \cr
#' \code{disp} \tab \code{<dbl>} \tab Value of the dispersion \cr
#' \code{A} \tab \code{<dbl>} \tab Value for A (the lower plateau). If 'NA' this is not a differentially abundant transcript. \cr
#' \code{k} \tab \code{<dbl>} \tab Value for k (the upper plateau). If 'NA' this is not a differentially abundant transcript. \cr
#' \code{seed_number} \tab \code{<int>} \tab  Seed number used for random generation of transcript abundance from mean transcript abundance and dispersion size. \cr
#' }
#' 
#' 
#' @param n_true_tests The number of true positive transcripts 
#' (i.e. truly differentially abundant transcripts) 
#' \code{ (<int> >= 0) } 
#' @param n_false_tests The number of null transcripts (i.e. not differentially abundant transcripts) 
#' \code{(<int> >= 0)} 
#' @param beta Value of beta parameter in the Richard's curve equation. 
#' (inversely proportional to inflection). Applicable to differentially abundant transcripts only.
#' \code{(<dbl> != 0)}
#' @param k Value for k parameter (the upper plateau) in Richard's curve sigmoid equation. 
#' Applicable to differentially abundant transcripts only.
#' \code{(<dbl> > 0)} 
#' @param A  Value for A parameter (the lower plateau) in Richard's curve sigmoid equation. 
#' Applicable to differentially abundant transcripts only. 
#' \code{(<dbl> >=0)}
#' @param sample_size Number of samples per transcript. Can consider this the number of 
#' patients or samples in this dataset. \code{(<int> > 0)} 
#' @param fixed_disp Is the overdispersion a single fixed value for this dataset (default, TRUE), 
#' or linearly related to the mean transcript abundance (FALSE)
#' @param disp_size Value for the overdispersion, 
#' used to simulate from the negative binomial distribution. 
#' If \code{fixed_disp = TRUE} (default), then this is a single value \code{(<dbl> > 0)}
#' Else if \code{fixed_disp = FALSE}, then this should be a vector of two values. The first value is 
#' a value for the slope, and the second is a value for the y-intercept, for the linear relationship between
#' overdispersion and transcript mean abundance.  \code{ ( vector,  c(<dbl> < 0, <dbl> > 0) )}
#' @param covar_discrete Is there a discrete number of covariate values (TRUE, default) or is the x-coordinate (covariate) continuous (FALSE)?
#' @param covar Covariate values to simulate over. 
#' If covariate is discrete then vector of all possible covariate values  
#' (by default it is \code{seq(from = -5, to = 5, by = 0.5)}). 
#' Else if the covariate is continous, it should be a two-element vector to 
#' be used as min and max values for a continuous uniform distribution.  \code{ ( vector,  c(<dbl>, <dbl>) )}
#' @param alpha The range of alpha parameter values in the Richard's curve sigmoid equation. 
#' (proportional to inflection). 
#' By default, it is set as all possible alpha values which would keep inflection within the range
#' of the covariate values. \code{( = "all")}
#' Otherwise, it can be set as a range (a two-element vector be used as min and max values 
#' for a continuous uniform distribution).  \code{ ( vector,  c(<dbl>, <dbl>) )}
#' Or it can be set as a single value. \code{ (<dbl>)}
#' @param null_mean_distribution Vector of mean transcript abundance to sample from for null test 
#' (By default; based on values of mean for a single transcript and CAPRA-S score from TCGA prostate cancer dataset)
#' @param seed_n Seed number for generating read count from a negative binomial distribution. 
#' Let = FALSE if you want to randomly select the seed used, else give numerical value. 
#' The actual seed used will be returned as a column in the final return data frame, so you can always replicate a previous dataset. 
#' 
#' @return Returns a data frame of simulated transcripts with a 
#' specified number of true positive curves (with varied inflection value) 
#' and specified number of null curves.
#' 
#' @examples 
#' 
#' # A dataset of 40 patients, with 100 differentially abundant transcripts. 
#' # Each patient is scored on a pseudo-continous (discrete) 
#' # covariate, with possible values, seq(from = 0, to = 7, by = 1)
#'                 
#' sigmoidal_sim_df(n_true_tests = 100,
#'                  n_false_tests = 0,
#'                  k = 8,
#'                  A = 3,
#'                  beta = 1,
#'                  covar = seq(from = 0, to = 7, by = 1)
#'                  disp_size = 0.85, 
#'                  sample_size = 40, # 5 patients will be assigned for each covariate value
#'                  seed_n = 50) 
#'  
#' # As above but this time a dataset of only null transcripts
#' 
#' sigmoidal_sim_df(n_true_tests = 0,
#'                  n_false_tests = 100,
#'                  covar = seq(from = 0, to = 7, by = 1),
#'                  disp_size = 0.85, 
#'                  sample_size = 40,
#'                  seed_n = 50)
#'            
#' # As above, but this time 41 patients
#' # Notice this time the number of patients is not a multiple of the 
#' # number of discrete covariate values. 
#' # In this case, patients will be randomly assigned
#' # covariate values using sampling with replacement. 
#'            
#' sigmoidal_sim_df(n_true_tests = 0,
#'                  n_false_tests = 100,
#'                  covar = seq(from = 0, to = 7, by = 1),
#'                  disp_size = 0.85, 
#'                  sample_size = 41,
#'                  seed_n = 50)               
#'                  
#' # A dataset of 50 patients, with both differentially abundant 
#' # and null transcripts. 
#' # Each patient has been scored on a continous covariate, 
#' # which ranges from 0 to 7.   
#'                
#' sigmoidal_sim_df(n_true_tests = 300,
#'                  n_false_tests = 300,
#'                  k = 7,
#'                  A = 2,
#'                  beta = 3,
#'                  covar_discrete = FALSE,
#'                  covar = c(0, 7),
#'                  disp_size = 0.85, 
#'                  sample_size = 40,
#'                  seed_n = 70)
#' 
#'
#' 
#' @export
sigmoidal_sim_df<-function(n_true_tests, #Number of True Positives (integer)
                          n_false_tests, #Number of Null Curves (integer)
                          beta = NA, #Value of beta for equation (double)
                          k = NA, #k value ( k + a is upper plateau) (double > 0)
                          A = NA, #A value (lower plateau level) (double >=0)
                          sample_size, #Number of samples per gene - needs to be a multiple of covar length (integer)
                          fixed_disp = TRUE, 
                          disp_size, #value of dispersion - used to simulate for neg binomial distribution (double)
                          covar_discrete = TRUE, # is there a discrete number of x-coordinates (TRUE) or is the x-coordinate continuous?
                          covar = seq(from = -5, to=5, by = 0.5), #If x-coordinate is discrete then: vector of x-coordinates of samples (simulated range of CAPRA_S values)
                          # If x-coordinate is continuous then: a vector of length two where the first element is the minimum possible value and the maximum is the second element
                          alpha = "all", # range of alpha values to simulate (proportional to inflection)
                          #by default all possible x-coordinates (- vector with min value, and max value - 
                          #otherwise a two-element vector be used as min and max values for a continuous uniform distribution)
                          null_mean_distribution = TCGA_mean_distribution, #distribution of means for a null test (based on values of mean for a single gene and CAPRA score from TCGA)
                          seed_n = 30 #Let = FALSE if you don't want to select a random seed used, else give numerical value 
                          # If simulating multiple datasets - and wanting different x-coordinate values
                          # then seed_n should be false - the seed number used will be returned in the data frame for reproducibility
                          ){ 
  
  
  
  ############################################################################
  
  # Set up 
  
  box::use(reshape2[...],
           dplyr[...],
           data.table[...],
           magrittr[...], 
           stats[...])
  
  
  # require(reshape2)
  # require(stats)
  # require(dplyr)
  # require(data.table)
  
  # Distribution of mean values for null tests in simulation 
  # Take at the distribution of means from TCGA - for gene and CAPRA_S i.e. 
  
  load("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Copy_Of_Article_Sections/TCGA_Prostate_Simple.rda")
  
  TCGA_mean_distribution = TCGA %>% 
    group_by(transcript, CAPRA_S) %>% 
    summarise(mean = mean(read_count_normalised)) %$%
    mean
  

  
  # Total number of genes / transcripts / tests to simulate 
  
  n_tests = (n_true_tests) + (n_false_tests)
  
  
  # If a random seed number is not set, randomly select seed number to set (from numbers 1 to 10000)
  
  if(seed_n == FALSE){ 
    seed_n = sample.int(10000, size = 1)
  }
  
  
  ###############################################################################
  
  # a. Set up / assign x-coordinate values 
  
  ##############################################################################
  
  
  
  
  #Number of distinct x-coordinate values ("simulation" CAPRA_S values)
  # covar<-seq(from = -5, to=5, by = 0.5) #21 different values of CAPRA_S
  n_cord<-length(covar)
  

  # Condition if: 
 
  # If the simulated dataset is over a continuous x-axis then simulate x-coord values using a uniform distribution
  # taking the original covar to be 
  #vector of length two where the first element is the minimum possible value and the maximum is the second element
  
  if(covar_discrete != TRUE){
    
  # Check that covar in is in the correct form - else stop and return error
     if(length(covar)!=2|covar[2]<covar[1]){
       stop(cat("The x-coordinate was given as continous, but the supplied covar value was not in the correct form. \n For continous x-coordinates, covar should be a vector of length two where \n the first element is the minimum possible value and the maximum is the second element. \n \n"))
     }
    
    # Warn that in the case of a continuous x-coordinate, 
    cat(" \n Continous x-coordinates are simulated using a uniform distribution. \n Samples could be unbalanced along the x-coordinate.\n")

# Simulate continous x-coordinates    
# Set of covar values which are going to correspond to the x values for each sample
    covar = stats::runif(n = sample_size, min = covar[1], max = covar[2])
  }
  

    # Condition if: 
  
   # If x-coordinate is discrete and sample size is not a multiple of the supplied distribution of x-coordinates
   # simulate sample covar values from the supplied discrete distribution of x-coordinates (giving a warning)
  
  else{
    if(sample_size %% n_cord != 0){
      
      # If sample size is not a multiple of the supplied distribution of x-coordinates
      # give warning 
      cat("Warning: Sample size is not a multiple of the simulated range of x-coordinates. \n There will likely be an uneven number of samples assigned to each x-coordinate.\n \n")
      
      cat("Samples were assigned x-coordinate values by random sampling with replacement the supplied discrete distribution of x-coordinates. \n \n")
      set.seed(seed_n)
      
      
      # Set of covar values which are going to correspond to the x values for each sample
      covar = sample(covar, 
                      size = sample_size, 
                      replace = TRUE)
      
    }
    
   # Condition if:  
    
   else{ # If x-coordinate is discrete and sample size is a multiple of the supplied distribution of x-coordinates
     
     # Then simulate an equal number of 
     # Number of samples per x-coordinate value 
     n_sam_cord = sample_size/n_cord
     
     # Replicate each x-coordinate value, so that each x-coordinate value is has n_sam_cord samples
     covar = sapply(covar, 
                     function(x) rep(x, n_sam_cord)) %>% c()
     
     } 
  }
  
 
  # Final covar values for samples in this data set
  # Arrange by value 
   
  covar = covar %>% 
           sort()
  
    #         #          #          #         #           #
  
  
  
  ###################################################################
  
  
  
  # Setting up alpha values 
  
  
  
  ##################################################################
  
  
  if(n_true_tests!=0) { #alpha is only applicable for DE genes - so if no DE genes ignore
    
  if(is.na(beta)|is.na(A)|is.na(k)){ #If any of the required parameter values for the sigmoid curve
    stop("One or more of the required parameter values (beta, A, k) for differentially abundant transcripts was not set.
         Either set these values, or don't simulate any differentially abundant transcripts (by setting n_true_tests = 0).")
    
  }
  #Creating the sigmoidal curve based genes (true tests - diff transcribed genes)
    
  
    
  #################################################################
  
  
    
  #Alpha is linearly dependent on the value of the inflection, inflection = - alpha / beta
  # so -inflection*beta = alpha
  
    
  #Alpha is set to simulate such that inflection values are within all possible values of the x-coordinate range
  
  
  # Condition if:    
    
  # If simulating over all possible values of alpha
  # then simulate from uniform distribution, with 
  # max = beta*max(x-coord)
  # min = beta*min(x-coord)
    
  if (alpha[1] == "all") {
  # simulate inflections for the full range of x coordinate values  
 
  set.seed(seed_n) 
    
  alpha = runif(min = min(covar)*abs(beta), 
                 max = max(covar)*abs(beta),
                 n = n_true_tests) 
  
  } 
    
  # Condition if 
  # Alpha is provided a range 
    
  else{ # Else if a range of alpha is specified use that range 
    
    # Give warning if range is not in desired format, i.e. not a length two vector c(alpha_min, alpha_max)
    
    if(length(alpha)!=2|alpha[2]<alpha[1]){
      stop(cat("The alpha was given as a range, but the supplied alpha vector was not in the correct form. \n For alpha as a range, alpha should be a vector of length two where \n the first element is the minimum possible value and the maximum is the second element \n \n"))
    }

    # Explain how alpha is being simulated
    else{cat("Alpha values were simulated using a uniform distribution with the supplied range \n \n")
      
    set.seed(seed_n) 
    
    alpha = runif(min = min(alpha),
                  max = max(alpha),
                  n = n_true_tests)
    }
  }
    
    
    
    #########################################################################################
    
    
    
    
    # Simulate true tests 
    
    
    
    
    #########################################################################################
    
    
    
    
    
  #Using the sigmoidal equation with 
  #values as defined in the function
  
  sig_eq<-function(x, Alpha) {
    
    #Lower plateau = A
    #Upper plateau = k 
    #Alpha adjusts inflection
    #beta adjusts slope
    
    k = k
    alpha_val = Alpha
    beta = beta
    A = A
    
    
   sig_eq_res = c(x, 
                  alpha_val, 
                  (k-A)/(1+exp(
                                -(alpha_val + x*beta))
                                                        ) + A
                                                               )  
   
    
   names(sig_eq_res) = c("CAPRA_S",
        "alpha",
        "log_y_hat")
      
      return(sig_eq_res) 
  }
  
  #Use the sigmodail eq  above to simulate log_y_hat values 
  #(i.e. mean values for each x coordinate values )
  
  #For log_y_hat list
  #Each list element is a gene
  #Top row is x-coordinate values / CAPRA-S values
  #Middle row is 
  #Bottom row is a simulated mean (log_y_hat) value
  

  
  log_y_hat_list = lapply(alpha, function(a)
    sapply(covar, #For better behaviour of simulated curve, centre values around 0 
                #Replacating normalisation process undertaken by TABI
                function(x)
                  sig_eq(x,a)
  )) 
  
  
  
  # Then convert list to a data frame with an id number column for each gene 
  # and rows converted to columns
  log_y_hat_df = log_y_hat_list %>% 
    lapply(., function(x) t(x) %>% 
                                   as_tibble()) %>% 
                                                data.table::rbindlist(.,idcol = "id")
  
  true_tests = log_y_hat_df %>% 
    mutate(Gene_number = paste0("X", id)) %>% 
    rowwise() %>% 
    mutate(value = rnbinom(mu = exp(log_y_hat), #From each mean value of each x coordinate (log_y_hat), Simulated Negative Binomial Distributed Values
                size = disp_size, #Precision/overdispersion values
                n = 1)) %>% 
    ungroup() %>% 
    mutate(sample_id = rep(1:(sample_size), n_true_tests)) %>% 
    arrange(id, sample_id) %>% 
    select(-id) %>% 
    mutate(Null_test = "FALSE")

  
  # Add parameter about the simulated curves information - 
  true_tests = true_tests %>% 
    mutate(sample_size = sample_size) %>% 
    mutate(A = A) %>% 
    mutate(slope = beta) %>% 
    mutate(k = k) 
  
  } 


  #         #          #          #         #           #
  
  
  #Simulation of null tests (negbinom distributed, with varying means  - same dispersion)
  #Similar procedure as above, but mean remains constant (as not diff transcribed)
  
  # n_null<-n_tests-n_true_tests #Number of null tests
  
  # if there are no null tests - then the simulated for the final dataframe is complete 
  if(n_false_tests == 0) {
    final_df<-true_tests
  }
  
  else{ # Else simulate null genes / transcripts as well 
    
  set.seed(seed_n) 
  #Simulation a n_null number of means - any integer between 0 and 10000
  
  null_means = data.frame(log_y_hat = sample(log(null_mean_distribution), 
                        size = n_false_tests, 
                        replace = T), 
                        id = 1:n_false_tests,
                         Gene_number = paste0("V",
                                              1:n_false_tests))
  
  
  null_df = lapply(1:sample_size,
                   function(x) {
                     null_means %>% 
                     mutate(sample_id = x) %>% 
                     mutate(CAPRA_S = covar[x]) 
                     }) %>% 
            bind_rows() 
  
  #Simulate null_values 
  
  null_tests = null_df %>% 
                mutate(value = rnbinom(mu = exp(log_y_hat), 
                                      size = disp_size,
                                      n =1)) %>% 
    mutate(Null_test = "TRUE") %>% 
    mutate(sample_size = sample_size) %>% 
    arrange(id, sample_id) %>% 
    select(-id)

  #Return a data frame with CAPRA_S (x-coordinate), 
  #Gene_number (Gene identifier column), value (simulated gene expression value)
  #Null_test (identifying true vs null tests), Alpha and Sample (Sample Identifcation Column)
  #All True Tests are Differentiated by having a "X" in Gene_number
  if(n_true_tests != 0) {
    final_df<-bind_rows(true_tests, 
                    null_tests)
  }
  
  else{
    final_df = null_tests
  }

  }
  
  # Add all the required information to the final data frame
  # Inflection, dispersion value, Gene ID column and seed number
  return(final_df %>% 
           rowwise() %>% 
           mutate(seed_number = seed_n) %>% 
           mutate(disp = disp_size) %>% 
           mutate(inflect = ifelse(Null_test == TRUE,
                                  NA,
                                  -1*as.numeric(alpha)/(slope))) %>%  #inflection is = -alpha/beta - only include for genes which are truly DE
           mutate(Gene_ref = ifelse(Null_test == TRUE, # Add gene reference / identification column
                                    paste0(Gene_number, #If gene is null - then gene id includes gene number sample size, dispersion, seed_number
                                           sample_size,
                                           disp_size,
                                           seed_number),
                                    paste0(Gene_number, # Else if DE gene - then gene id also includes parameter values 
                                    slope, 
                                    alpha, 
                                    k, 
                                    A, 
                                    sample_size, 
                                    disp_size, 
                                    seed_number)
                                    ) 
                  ) %>% 
           mutate(Null_test = as.logical(Null_test)) %>% 
           mutate(value = as.integer(value)) %>% 
           mutate(seed_number = as.integer(seed_number)) %>% 
           select(sample_id,
                  Gene_ref,
                  Gene_number,
                  Null_test,
                  CAPRA_S,
                  value,
                  log_y_hat, 
                  slope,
                  inflect, 
                  sample_size,
                  disp, 
                  A,
                  k,
                  seed_number
                  )
         ) 
  
}



##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 




#(Section 2)
# Create function to plot random gene from simulated data frame


# Create Plotting function 
# Which takes a simulated dataframe 
# Selects a random simulated gene
# And Plots it (listing relevant features)



##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



#' Function to plot random gene from simulated data frame (from sigmoid_sim_df)
#' 
#' 
#' @param sim_df Data frame of simulated read counts 
#' (outcome of \code{sigmoidal_sim_df})
#' 
#' @param type  What type of transcript do you want to randomly select - 
#' differentially abundant/expressed (let = "DE"), null (= "null"), or "any" 
#' 
#' @export
plot_rsim_df<-function(sim_df, # Dataframe of simulated read counts (outcome of sigmoidal_sim_df function)
                       type = "any") { # What type of gene do you want to randomly select - diff expressed (let = "DE"), null ("null"), or "any" 
  
  box::use(dplyr[...],
        magrittr[...],
        ggplot2[...])
  
  # If 
  if(type == "DE") {
    
  gene_name<-sim_df %>%  
    select(Gene_ref) %>% 
    filter(grepl("X", Gene_ref))  %>% 
    ungroup() %>% 
    distinct() %>% 
    sample_n(1) %$%
    Gene_ref 
  } 
  
  else if(type == "null") { 
    
    gene_name<-sim_df %>%  
      select(Gene_ref) %>% 
      filter(grepl("V", Gene_ref)) %>% 
               ungroup() %>% 
               distinct() %>% 
               sample_n(1) %$%
               Gene_ref
    }
  
  else {
    
      gene_name<-sim_df %>%  
        select(Gene_ref) %>% 
        ungroup() %>% 
        distinct() %>% 
        sample_n(1) %$%
        Gene_ref
    
  }
  
  # If the gene selected is a truly DE gene (with identifier X and not V)
  # take important parameter attributes for this e.g. slope, sample size etc. 

  if(grepl("X", paste0(gene_name))) {
    att= sim_df %>% 
       ungroup() %>% 
       filter(Gene_ref == gene_name) %>%
      select(sample_size, slope, inflect) %>%
      mutate(slope = paste(slope, "inflection  = ", inflect)) %>% 
      select(sample_size, slope) %>% 
      distinct() }

 # Otherwise if null gene, then label as such 
  # and only take sample size as an important parameter
  else {
    att= sim_df %>% 
      ungroup() %>% 
      filter(Gene_ref == gene_name) %>% 
      select(sample_size) %>% 
      distinct() %>% 
      mutate( slope = "No Slope (Null Test)")
  }
  

 plot = sim_df %>% 
    ungroup() %>% 
    filter(Gene_ref == gene_name) %>% 
    mutate(Test_type = ifelse(Null_test == TRUE, "Null", "DE")) %>% 
    ggplot(aes(x=CAPRA_S, 
               y=value + 1,
               Sample_ID = sample_id,
               Gene_number = Gene_number,
               Test_type = Test_type)) + 
    geom_point(col = "dodgerblue", size = 2) + 
    scale_y_log10() + 
    labs(subtitle = paste("Sample size = ", att$sample_size,
                            "Slope = ", att$slope), 
         y = "Simulated normalised RNA seq counts + 1",
         x = "Simulated CAPRA-S values (x-coordinate)",
         title = paste0("Example Simulated Gene,",
                        gene_name)) + 
   theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
         panel.background = element_blank(),
         panel.grid.major = element_line(colour = "grey84"), 
         panel.grid.minor = element_line(colour = "grey84", linetype = 2)) + 
   list(
     #scale_fill_manual(values = friendly_cols,   na.value = "black"),
     #scale_color_manual(values = friendly_cols,   na.value = "black"),
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
 
  
  print(plot)
  
}








##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



#(Section 3) make_groups function
# Function which is able to make a data frame to match CAPRA_S values to group labels
# Given a desired number of groups to sort the simulated samples into 

# (Section 1)

# Make groups for using tidybulk differential expression


# Input simulated df table from function sigmoidal_sim_df 

# Output: list consisting of $group (df matching CAPRA-S / x coordinate value to group letter)
# and $contrast_values (vector of values which correspond to the %>% ddle value between each group  - 
# i.e. if highest CAPRA-S for group A is 4 and lowest for group B is a 5 then the contrast value is 4.5)


#' Stratify samples into groups, by covariate value, 
#' 
#' 
#' To use for tidybulk to test differential transcript abundance along a continous or pseudo-continous value
#' 
#' @description 
#' \code{make_groups} takes a dataset produced from function \code{sigmoidal_sim_df} and returns a list  
#' differentially abundant transcripts, with respect to a single continous or pseudo-continous variable.
#' The log mean at each covariate value (\code{log_y_hat}) for 
#' differentially abundant transcripts is set from a Richard's curve. 
#' 
#' 
#' 
#' @param data Data frame of of simulated RNA seq data counts (from \code{sigmoidal_sim_df})
#' @param .sample Name of the column which identifies the sample (e.g. Patient / Sample / Biological Replicate ID)
#' By default, this is \code{sample_id} (this is the typical output column from \code{sigmoidal_sim_df})
#' @param .transcript Name of the column which identifies the transcript. 
#' By default, this is \code{Gene_ref}, the unique reference column to a simulated transcript (from \code{sigmoidal_sim_df})
#' @param .abundance Name of the column 
#' @param .covar Name of the column of interest, i.e. the covariate. By default it is \code{CAPRA_S}
#' @param n_groups The number of groups to stratify samples into 
#' @param covar_range
#' 
#' @return A list consisting of: 
#'
#' - $group (A data frame which matches covariate value (e.g. CAPRA-S) to a letter
#' representing a group)
#' 
#' - $contrast_values (A vector of values which correspond to the middle value between each group  - 
#' i.e. if highest CAPRA-S for group A is 4 and lowest for group B is a 5 then the contrast value is 4.5)
#' 
#' 
#' 
#' @export
make_groups<-function(data, #Table of simulated RNA seq data counts
                      .sample = sample_id, #Sample ID column 
                      .transcript = Gene_ref, #transcript / gene ID column
                      .abundance = value, # read count columns
                      .covar = CAPRA_S,# x - coordinate column / factor of interest, in this case CAPRA_S
                      covar_discrete = TRUE, 
                      n_groups, # number of groups 
                      covar_range = seq(from = -5, to =5, by = 0.5) #For grouping, vector of discrete x-coordinates values / simulated CAPRA values
){
  
  
  # Set up columns
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  .covar = enquo(.covar)
  
  # Number of distinct simulated x-coordinate values 
  n_cord = length(x_cord)
  
  # Reduce df of simulated data to only required samples
  # sample id column (sample)
  # name of transcript 
  # RNAseq count column
  data = data %>% 
    select(!!.sample,
           !!.transcript,
           !!.abundance)
  
  if (!covar_discrete) {
    
    
    
  }
  
  
  if((n_cord %% n_groups) ==0){ #If number of groups to make is a multiple of possible x_cord values
    
    # Then number of different CAPRA S value or x-coordinate values with per groups is : 
    n_per_group<-n_cord/n_groups # Number of different x-cord per group 
    
    # Vector of number of samples in each group (assuming one sample per CAPRA-S value)
    groups<-rep(n_per_group, 
                n_groups)
    
    x_cord_id<-cumsum(groups)[-length(groups)]
    
    # Contrast values are the point estimate of the 'jump' 
    # i.e. if there is a change detected between the two groups, what is the estimate
    contrast_values<-(x_cord[x_cord_id]+x_cord[x_cord_id+1])/2
    # # Assign each sample a group letter as a group id
    
    
    # Assign each x-coordinate value a group letter based on how many x-coordinate values there are per group
    group_let<-sapply(LETTERS[1:n_groups], 
                      function(x) rep(x, n_per_group)) %>% 
      c()
    
    
    
    # Combine group id and corresponding x-cord value into a data frame
    df<-data.frame(group = group_let,
                   CAPRA_S = x_cord)
    
    
  } 
  else{ #if the number of groups is not a multiple of possible x_cord values
    # Then group 'lower' x-coordinate values together, such that there are bigger groups
    # with earlier CAPRA-S values
    
    n_per_group<-n_cord %/% n_groups # How many samples per group, if only accept whole numbers
    
    n_left_over<-n_cord %% n_groups # How many samples would be left over if we used the above number
    
    # Vector of number of samples in each group (assuming one sample per CAPRA-S value)
    # Give early / left most groups add an additional sample
    # than late / right most groups
    groups<-c(rep(n_per_group+1,n_left_over), 
              rep(n_per_group, n_groups-n_left_over))
    
    
    
    group_let<-sapply(1:n_groups, 
                      function(x) rep(LETTERS[x], groups[x])) %>% 
      Reduce(c, .)
    
    
    x_cord_id<-cumsum(groups)[-length(groups)]
    
    
    # Contrast values are the point estimate of the 'jump' 
    # i.e. if there is a change detected between the two groups, what is the estimate
    contrast_values<-(x_cord[x_cord_id]+x_cord[x_cord_id+1])/2 
    
    
    # Combine group id and corresponding x-cord value into a data frame
    df<-data.frame(group = group_let,
                   CAPRA_S = x_cord)
    
  } 
  
  
  # Return data frame of group id and corresponding x-cord value
  # And point estimate the point estimate of the 'jump' between two groups
  
  return(list(group = df,
              contrast_values = contrast_values)) 
  
  
}



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


