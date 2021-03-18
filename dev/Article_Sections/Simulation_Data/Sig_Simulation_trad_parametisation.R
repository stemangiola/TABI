#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #



#Intention of Script

#(Section 1)
#Test function to create sigmodial curve data 

#(Section 2)
#Create functions to simulate a dataset
# of a mixture of sigmoidal based distributions (true associations)
# and null associations 

#(Section 3)
# Use this function to simulate data 



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #
#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #




##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

# Section 1: 

#Create Generalised Sigmoid Function Equation

#sig_eq is GLA eqaution which log y hat (log of the mean of gene read count) should be equiv to

#Plot examples as sanity check

#Note: 
#Typical parametisation of sigmodial function, 
#i.e. in form k/(1+exp(xb-a)) + A is used to simulate sigmoidal equations
# (re-parametisation form which used used for inference by TABI
#is less useful for simulation as changing eta, (inflection x-coordinate)
# with other parameters remaining the same causes wildly different upper and lower plateaus,
# which unintentionally leads to unrealistic simulation data) 


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 


#Show example: 

sig_eq<- function(x) {
  
  #Lower plateau = A
  #Upper plateaa = k + A
  #Alpha adjusts inflection
  #beta adjusts slope
  
  k = 3
  alpha = 0
  beta = 2
  A = 5
    
    
  return(k/(1+exp(-(alpha + x*beta))) + A) 
}

library(dplyr)

library(ggplot2)


#Plot of sigmoid equation
data.frame(x = seq(from = -5, to = 5, by = 0.5)) %>% 
  ggplot(aes(x=x)) + 
  stat_function(fun=sig_eq, geom="line") +
  labs(y = "Log Mean Gene Count",
       x = "CAPRA-S")


# Example of simulated data 

#Use the eq above to simulate log_y_hat values 
log_y_hat<-sapply(seq(from = -5, to=5, by = 0.5),
                  function(x)
                    sig_eq(x)  
)


#Simulated Negative Binomial Distributed Values
sample<-sapply(log_y_hat,
               function(y) rnbinom(mu = exp(y), #Mean is exp of log y hat (GLA eq)
                                   size = 0.25, #Precision/overdispersion value of 0.3 seems approp. level of noise
                                   n = 10))

plot_sig<-function(x) {
  exp(sig_eq(x))
}

f<-data.frame(x = seq(from = -5, to = 5, by = 0.5), value = log_y_hat %>% exp())

#Plot simulated data 
data.frame(
  CAPRA_S= c(
    sapply(seq(from = -5, to=5, by = 0.5), 
           function(x) rep(x, 10))), 
  sample = c(sample)
) %>% 
  ggplot() + 
  geom_point(aes(x=CAPRA_S, y=sample+1), col="dodgerblue")  +
  #labs(title="True Positive Curve", x="CAPRA-S", y="Gene Count") + 
  labs(x = "", y= "") + 
  geom_line(aes(x=x, y=value), data = f) + 
  stat_function(fun=plot_sig, geom="line") + 
  scale_y_log10()  + 
  theme(panel.grid.major = element_line(colour='white', linetype = 1),
        panel.grid.minor= element_line(colour='white', linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) %>% 
  labs(main = "Example Simulated Gene ")

##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 


#Section 2 


#Create simulation function


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 


#Function for simulating Data with both null and true (sigmodial curve) tests


#Function returns table with simulated data from
#specified number of true positive curves - 
#of a single equation but varying inflection value 
#And specified number of null curves



sigmodial_tests<-function(n_true_tests, #Number of True Positives
                          n_false_tests, #Number of Null Curves
                          beta, #Value of beta for equation
                          k, #k value for slopes 
                          A, #A value (lower plateau level)
                          sample_size, #Number of samples per gene - sample size needs to be a multiple of x_cord length
                          disp_size, #Size of dispersion - used to simulate for negbinomial distribution 
                          x_cord = seq(from = -5, to=5, by = 0.5) #Calculation of x-coordinates ("simulation" CAPRA_S values)
                          #(for rnbiom(mu, size = disp_size))
                          ) { #With the way this function is coded Sample size needs to be multiple of 21
  
  require(reshape2)
  require(dplyr)
  
  
  #Calculation of x-coordinates ("simulation" CAPRA_S values)
  # x_cord<-seq(from = -5, to=5, by = 0.5) #21 different values of CAPRA_S
  n_cord<-length(x_cord)
  
  
  n_tests = (n_true_tests) + (n_false_tests)
  

    #         #          #          #         #           #
  
  
  
  #Creating the Sigmodial curve based genes (true tests - diff transcribed genes)
  
  
  #Simulation all possible values of the inflection (alpha)
  #Alpha is linearly dependent on the value of the inflection
  
  
  #If the number of True Positive Curves is less than 22
  #Then only simulate for less extreme inflections (between -2/2)
  #Otherwise simulate inflections for the full range of x coordinate values
  
  
  if (n_true_tests<=21) {
  set.seed(35)
  alpha = sample(seq(from = -2, to = 2, by =1), 
                   size = n_true_tests, replace = T)
  }
  
  else{
  set.seed(35)
  alpha = sample(seq(from = -5, to = 5, by =0.5), 
               size = n_true_tests, replace = T) }
  
  
  #Using the sigmodial equation with 
  #values as defined in the function
  
  sig_eq<-function(x) {
    
    #Lower plateau = A
    #Upper plateau = k + A
    #Alpha adjusts inflection
    #beta adjusts slope
    
    k = k
    alpha = alpha
    beta = beta
      A = A
      
      
      return(k/(1+exp(-(alpha + x*beta))) + A) 
  }
  
  #Use the sigmodail eq  above to simulate y_hat values 
  #(i.e. mean values for each x coordinate values )
  
  y_hat<-sapply(seq(from = -5, to = 5, by =0.5), #For better behaviour of simulated curve, centre values around 0 
                #Replacating normalisation process undertaken by TABI
                function(x)
                  sig_eq(x)  
  )
  
  #For y_hat
  #Each column is a for a value of CAPRA_S 
  #Each row is a simulated curve
  
  
  #From each mean value of each x coordinate
  #Simulated Negative Binomial Distributed Values
  sample<-apply(y_hat, 
                1, #For each row (a curve)
                function(x) 
                  sapply(x, #Simulate sample size/15 
                         function(y) rnbinom(mu = exp(y), 
                                             size = disp_size, #Precision/overdispersion values
                                             n = sample_size/n_cord))
                
  )
  
  #In Sample
  #Each column is a simulated curve / gene 
  #Each row is a sample
  
  
  #Combine respective CAPRA_S values with their respective simulated gene counts
  True_tests<-data.frame(
    CAPRA_S= c(
      sapply(seq(from = -5, to=5, by = 0.5), 
             function(x) rep(x, sample_size/n_cord))), #
    sample
  )
  
  #Reshape the above such that it mimics TCGA data set (column labeling gene count, column labeling gene name,
  #column labbeling CAPRA_S value)
  
  true_tests<-True_tests %>% 
    melt(id.vars = 1) %>% #Melt with respect to CAPRA_S 
    mutate(Null_test = "FALSE") %>% #Add a column which indicates this is a null 
    rename(Gene_number = variable) #Make Gene_number the test identification column
  
  #Add inflection data to Gene_number
  inflection<-cbind(true_tests %>% 
    select(Gene_number) %>% 
    distinct(), alpha)
  
  #Add information about alpha and Sample number to final table
  true_tests=inner_join(true_tests,
                         inflection) %>% 
    mutate(Sample = rep(1:(sample_size), n_true_tests))
  
  #         #          #          #         #           #
  
  
  #Simulation of null tests (negbinom distributed, with varying means  - same dispersion)
  #Similar proceedure as above, but mean remains constant (as not diff transcribed)
  
  n_null<-n_tests-n_true_tests #Number of null tests
  
  #Simulation a n_null number of means - any integer between 0 and 15000
  null_means<-sample(0:10000, 
                     size = n_null, 
                     replace = T)
  
  #Simulate null_values 
  null_values<-sapply(null_means, 
                      function(x) 
                        (rnbinom(n=sample_size,
                                 mu = x, 
                                 size = disp_size)))
  
  null_tests<-data.frame(
    CAPRA_S= c(
      sapply(seq(from = -5, to=5, by = 0.5), 
             function(x) rep(x, sample_size/n_cord))), #
    null_values
  ) %>% 
    melt(id.vars = 1) %>% #Melt with respect to CAPRA_S 
    mutate(Null_test = "TRUE") %>% #Add a column which indicates this is a null 
    rename(Gene_number = variable) 
  
  
  null_tests$Gene_number<-sub("X", "V", as.factor(null_tests$Gene_number)) #Convert X test recognition to V to further
  #Distinguish from true tests (now true tests have an X prefix, and null have a V prefix)
  
  #Add column of alpha (not relevant so NA), and Sample Identification
  null_tests<-null_tests %>% 
    mutate(alpha = "NA") %>% 
    mutate(Sample = rep(1:(sample_size), (n_tests-n_true_tests)))
  
  #Return a data frame with CAPRA_S (x-coordinate), 
  #Gene_number (Gene identifier column), value (simulated gene expression value)
  #Null_test (identifying true vs null tests), Alpha and Sample (Sample Identifcation Column)
  #All True Tests are Differentiated by having a "X" in Gene_number
  
  return(rbind(true_tests, 
               null_tests) %>% 
           mutate(disp = disp_size) %>% 
           mutate(sample_size = sample_size) %>% 
           mutate(A = A) %>% 
           mutate(True_Slope = beta) %>% 
           mutate(k = k))
}


##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



#Section 3


#Create simulation dataset from function



##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 



library(foreach)
library(dplyr)



#Values of sigmodial equation to simulate over
slopes = setdiff(seq(from=-2, to =2, by =0.25), c(0)) #Set of slopes (beta values) to simulate
ks = seq(from =2, to = 4, by=1) #Set of k values to simulate
As = seq(from = 0, to =5, by =1) #Set of A values to simulate
disp = seq(from = 0.25, to = 4, by = 0.25)
sample_sizes = seq(from = 21, to = 105, by = 21)



#Create Simulation Dataset 
TABI_df<-foreach(k = ks, 
                 .combine = 'rbind') %do% { 
                       foreach(j=slopes,  
                               .combine = 'rbind') %do% {
          foreach(i = As,
                  .combine = 'rbind') %do% {
                    foreach(d= disp, 
                            .combine = 'rbind') %do% {
                              foreach(s = sample_sizes,
                                      .combine = 'rbind') %do% {
      
  sigmodial_tests(100, #Number true tests (per single set of slope, k and A values)
                  10, #Number of null tests / gene (per single set of slope, k and A values)
                  j, #Slope value
                  k, #k value
                  i, #A value
                  s, #Number of samples per simulated gene (default x_cord means multiple of 21)
                  d #Dispersion value (size in rnbinom)
                  ) %>% 
    filter(grepl("X", Gene_number))
                            }
                            }     
                  } 
        } 
        } 


#Total number of genes / tests in final simulated dataset
n_total_tests<-Sig_multiple_value_data %>% 
  select(Gene_ref) %>% 
  distinct() %>% 
  nrow()



##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ## 

