#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Intention of Script 

# Set up  /define functions used in all parts of analysis 


#(Section 1)
#Create function to simulate a dataset - sigmoidal_sim_df
# of a mixture of sigmoidal based distributions (true associations)
# and null associations 

#(Section 2)
# Create function to plot random gene from simulated data frame



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

#' @export
box::use(./Functions_Simulated_Sigmodial_df_RNAseq)



#Section 1 


#Create simulation function

# For example simulated data frames 
# And sigmoid functions see script for Pt 1_1 



##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 


box::use(magrittr[...], 
         dplyr[...])


# Distribution of mean values for null tests in simulation 
# Take at the distribution of means from TCGA - for gene and CAPRA_S i.e. 

load("/stornext/Home/data/allstaff/b/beasley.i/TABI/dev/Article_Sections/TCGA_Prostate_Simple.rda")

TCGA_mean_distribution = TCGA %>% 
  group_by(transcript, CAPRA_S) %>% 
  summarise(mean = mean(read_count_normalised)) %$%
  mean
  


#Function for simulating Data with both null and true (sigmoidal curve) tests


#Function returns table with simulated data from
#specified number of true positive curves - 
#of a single equation but varying inflection value 
#And specified number of null curves


# Function is designed to simulate DE genes with a single set of A, k, beta, sample size, and dispersion size
# over a set of alpha values (or potentially a single alpha value)
# To simulate of other parameter values - need to map function over a set of parameter values

#' @export
sigmoidal_sim_df<-function(n_true_tests, #Number of True Positives (integer)
                          n_false_tests, #Number of Null Curves (integer)
                          beta, #Value of beta for equation (double)
                          k, #k value ( k + a is upper plateau) (double > 0)
                          A, #A value (lower plateau level) (double >=0)
                          sample_size, #Number of samples per gene - needs to be a multiple of x_cord length (integer)
                          disp_size, #value of dispersion - used to simulate for neg binomial distribution (double)
                          x_cord_discrete = TRUE, # is there a discrete number of x-coordinates (TRUE) or is the x-coordinate continuous?
                          x_cord = seq(from = -5, to=5, by = 0.5), #If x-coordinate is discrete then: vector of x-coordinates of samples (simulated range of CAPRA_S values)
                          # If x-coordinate is continuous then: a vector of length two where the first element is the minimum possible value and the maximum is the second element
                          alpha = "all", # range of alpha values to simulate (proportional to inflection)
                          #by default all possible x-coordinates (- vector with min value, and max value - 
                          #otherwise a two-element vector be used as min and max values for a continuous uniform distribution)
                          null_mean_distribution = TCGA_mean_distribution, #distribution of means for a null test (based on values of mean for a single gene and CAPRA score from TCGA)
                          seed_n = 30 #Let = FALSE if you don't want to select a random seed used, else give numerical value 
                          # If simulating multiple datasets - and wanting different x-coordinate values
                          # then seed_n should be false - the seed number used will be returned in the data frame for reproducibility
                          ){ 
  
  box::use(reshape2[...], 
           dplyr[...],
           data.table[...],
           stats[...]) 
  
  # Total number of genes / transcripts / tests to simulate 
  
  n_tests = (n_true_tests) + (n_false_tests)
  
  
  # If a random seed number is not set, randomly select seed number to set (from numbers 1 to 1000)
  
  if(seed_n == FALSE){ 
    seed_n = sample.int(1000, size = 1)
  }
  
  
  
  ###############################################################################
  
  # a. Set up / assign x-coordinate values 
  
  ##############################################################################
  
  
  
  
  #Number of distinct x-coordinate values ("simulation" CAPRA_S values)
  # x_cord<-seq(from = -5, to=5, by = 0.5) #21 different values of CAPRA_S
  n_cord<-length(x_cord)
  

  # Condition if: 
 
  # If the simulated dataset is over a continuous x-axis then simulate x-coord values using a uniform distribution
  # taking the original x_cord to be 
  #vector of length two where the first element is the minimum possible value and the maximum is the second element
  
  if(x_cord_discrete != TRUE){
    
  # Check that x-cord in is in the correct form - else stop and return error
     if(length(x_cord)!=2|x_cord[2]<x_cord[1]){
       stop(cat("The x-coordinate was given as continous, but the supplied x_cord value was not in the correct form. \n For continous x-coordinates, x_cord should be a vector of length two where \n the first element is the minimum possible value and the maximum is the second element. \n \n"))
     }
    
    # Warn that in the case of a continuous x-coordinate, 
    cat(" \n Continous x-coordinates are simulated using a uniform distribution. \n Samples could be unbalanced along the x-coordinate.\n")

# Simulate continous x-coordinates    
# Set of x_cord values which are going to correspond to the x values for each sample
    x_cord = stats::runif(n = sample_size, min = x_cord[1], max = x_cord[2])
  }
  
  
  # Condition if: 
  
   # If x-coordinate is discrete and sample size is not a multiple of the supplied distribution of x-coordinates
   # simulate sample x-cord values from the supplied discrete distribution of x-coordinates (giving a warning)
  
  else{
    if(sample_size %% n_cord != 0){
      
      # If sample size is not a multiple of the supplied distribution of x-coordinates
      # give warning 
      cat("Warning: Sample size is not a multiple of the simulated range of x-coordinates. \n There will likely be an uneven number of samples assigned to each x-coordinate.\n \n")
      
      cat("Samples were assigned x-coordinate values by random sampling with replacement the supplied discrete distribution of x-coordinates. \n \n")
      set.seed(seed_n)
      
      
      # Set of x_cord values which are going to correspond to the x values for each sample
      x_cord = sample(x_cord, 
                      size = sample_size, 
                      replace = TRUE)
      
    }
    
   # Condition if:  
    
   else{ # If x-coordinate is discrete and sample size is a multiple of the supplied distribution of x-coordinates
     
     # Then simulate an equal number of 
     # Number of samples per x-coordinate value 
     n_sam_cord = sample_size/n_cord
     
     # Replicate each x-coordinate value, so that each x-coordinate value is has n_sam_cord samples
     x_cord = sapply(x_cord, 
                     function(x) rep(x, n_sam_cord)) %>% c()
     
     } 
  }
 
 
  # Final x-cord values for samples in this data set
  # Arrange by value 
   
  x_cord = x_cord %>% 
           sort()
  
    #         #          #          #         #           #
  
  
  
  ###################################################################
  
  
  
  # Setting up alpha values 
  
  
  
  ##################################################################
  
  
  if(n_true_tests!=0) { #alpha is only applicable for DE genes - so if no DE genes ignore
    
  
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
    
  alpha = runif(min = min(x_cord)*beta, 
                 max = max(x_cord)*beta,
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
    #Upper plateau = k + A
    #Alpha adjusts inflection
    #beta adjusts slope
    
    k = k
    alpha_val = Alpha
    beta = beta
    A = A
    
    
   sig_eq_res = c(x, 
                  alpha_val, 
                  k/(1+exp(-(alpha_val + x*beta))) + A)  
   
    
   names(sig_eq_res) = c("CAPRA_S",
        "alpha",
        "y_hat")
      
      return(sig_eq_res) 
  }
  
  #Use the sigmodail eq  above to simulate y_hat values 
  #(i.e. mean values for each x coordinate values )
  
  #For y_hat list
  #Each list element is a gene
  #Top row is x-coordinate values / CAPRA-S values
  #Middle row is 
  #Bottom row is a simulated mean (y_hat) value
  

  
  y_hat_list = lapply(alpha, function(a)
    sapply(x_cord, #For better behaviour of simulated curve, centre values around 0 
                #Replacating normalisation process undertaken by TABI
                function(x)
                  sig_eq(x,a)
  )) 
  
  
  
  # Then convert list to a data frame with an id number column for each gene 
  # and rows converted to columns
  y_hat_df = y_hat_list %>% 
    lapply(., function(x) t(x) %>% 
                                   as_tibble()) %>% 
                                                data.table::rbindlist(.,idcol = "id")
  
  true_tests = y_hat_df %>% 
    mutate(Gene_number = paste0("X", id)) %>% 
    rowwise() %>% 
    mutate(value = rnbinom(mu = exp(y_hat), #From each mean value of each x coordinate (y_hat), Simulated Negative Binomial Distributed Values
                size = disp_size, #Precision/overdispersion values
                n = 1)) %>% 
    ungroup() %>% 
    mutate(sample_id = rep(1:(sample_size), n_true_tests)) %>% 
    arrange(id, sample_id) %>% 
    select(-id) %>% 
    mutate(Null_test = "FALSE")

  #From each mean value of each x coordinate
  #Simulated Negative Binomial Distributed Values
  # sample<-lapply(y_hat, function(x) rnbinom(mu = exp(x), 
  #                                            size = disp_size, 
  #                                            n = 1))
                
  
  #In Sample
  #Each column is a simulated curve / gene 
  #Each row is a sample
  
  
  # #Combine respective CAPRA_S values with their respective simulated gene counts
  # True_tests<-data.frame(
  #   CAPRA_S= c(
  #     sapply(x_cord, 
  #            function(x) rep(x, 1))), #
  #   sample
  # )
  

  
  #Reshape the above such that it mimics TCGA data set (column labeling gene count, column labeling gene name,
  #column labbeling CAPRA_S value)
  # 
  # true_tests<-True_tests %>% 
  #   reshape2::melt(id.vars = 1) %>% #Melt with respect to CAPRA_S 
  #   mutate(Null_test = "FALSE") %>% #Add a column which indicates this is not a null test (i.e. it is DE)
  #   dplyr::rename(Gene_number = variable) #Make Gene_number the test identification column
  # 
  # #Add inflection data to Gene_number
  # inflection<-cbind(true_tests %>% 
  #   select(Gene_number) %>% 
  #   distinct(), alpha)
  # 
  # #Add information about alpha and Sample number to final table
  # true_tests=inner_join(true_tests,
  #                        inflection, by = "Gene_number") %>% 
  #   mutate(sample_id = rep(1:(sample_size), n_true_tests))
  
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
  
  null_means = data.frame(y_hat = sample(log(null_mean_distribution), 
                        size = n_false_tests, 
                        replace = T), 
                        id = 1:n_false_tests,
                         Gene_number = paste0("V",
                                              1:n_false_tests))
  
  
  null_df = lapply(1:sample_size,
                   function(x) {
                     null_means %>% 
                     mutate(sample_id = x) %>% 
                     mutate(CAPRA_S = x_cord[x]) 
                     }) %>% 
            bind_rows() 
  
  #Simulate null_values 
  
  null_tests = null_df %>% 
                mutate(value = rnbinom(mu = exp(y_hat), 
                                      size = disp_size,
                                      n =1)) %>% 
    mutate(Null_test = "TRUE") %>% 
    mutate(sample_size = sample_size) %>% 
    arrange(id, sample_id) %>% 
    select(-id)
  # null_values<-sapply(null_means, 
  #                     function(x) 
  #                       y_hat = (rnbinom(n=sample_size,
  #                                mu = x, 
  #                                size = disp_size)))
  # 
  
  # null_tests<-data.frame(
  #   CAPRA_S= c(
  #     sapply(x_cord, 
  #            function(x) rep(x, 
  #                            1))), #
  #   null_values
  # ) %>% 
  #   reshape2::melt(id.vars = 1) %>% #Melt with respect to CAPRA_S 
  #   mutate(Null_test = "TRUE") %>% #Add a column which indicates this is a null 
  #   dplyr::rename(Gene_number = variable) 
  
  
  # null_tests$Gene_number<-sub("X", "V", 
  #                             as.factor(null_tests$Gene_number)) #Convert X test recognition to V to further
  #Distinguish from true tests (now true tests have an X prefix, and null have a V prefix)
  
  #Add column of alpha (not relevant so NA), and Sample Identification
  # null_tests<-null_tests %>% 
  #   mutate(alpha = "NA") %>% 
  #   mutate(sample_id = rep(1:(sample_size), 
  #                          (n_tests-n_true_tests)))
  
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
           mutate(inflec = ifelse(Null_test == TRUE,
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

#' @export
plot_rsim_df<-function(sim_df, # Dataframe of simulated read counts (outcome of sigmoidal_sim_df function)
                       type = "any") { # What type of gene do you want to randomly select - diff expressed (let = "DE"), null ("null"), or "any" 
  
  require(dplyr)
  require(magrittr)
  require(ggplot2)
  
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
      select(sample_size, slope, inflec) %>%
      mutate(slope = paste(slope, "Inflection  = ", inflec)) %>% 
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
         panel.grid.minor = element_line(colour = "grey84", linetype = 2))
  
  print(plot)
  
}








##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

