

#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


#Intention of Script

#Create functions to simulate a dataset
# of a mixture of sigmoidal based distributions (true associations)
# and null associations (and situations of all tsrue and all false)
# and use this function to simulate data (for comparing TABI / other methods in another script) 



#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


# Section 1

#Testing potetional sigmodial curves (for realism)




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

#Create Generalised Sigmoid Function Equation

#eq is GLA euqation which log y hat (log of the mean of gene read count) should be equiv to

eq <- function(x) {
  #Testing Parameter values
  #Values of beta =1, y_0 = 2, A = 1 seem reasonable for eta from-2 to 2
  
  eta =  -1.19 #All values from -2/2 seem reasonable - eta is X coord inflection
  beta= -7.34 #Slope value
  y_0=0.7 #y_cross / y_0
  A = 2.86 #vertical translation
  
  top<- ((y_0-A)*(1+exp(-eta*beta))) #Numerator of the GLA function
  bottom<-(1+exp(-eta*beta+x*beta)) #Denominator of the GLA function
  
  return(exp(A+top/bottom)) #Evaluate the GLA function 
}

peq<-function(x) {
  #Testing Parameter values
  #Values of beta =1, y_0 = 2, A = 1 seem reasonable for eta from-2 to 2
  
  eta =  -0.25 #All values from -2/2 seem reasonable - eta is X coord inflection
  beta= 0.99 #Slope value
  y_0=0.05 #y_cross / y_0
  A = 7.35 #vertical translation
  
  top<- ((y_0-A)*(1+exp(-eta*beta))) #Numerator of the GLA function
  bottom<-(1+exp(-eta*beta+x*beta)) #Denominator of the GLA function
  
  return(exp(A+top/bottom)) #Evaluate the GLA function 
}




#Plot curve across all values of CAPRA_S (0-7) (in log space)

#To test that mean values of log_y_hat are reasonable 

library(ggplot2)

ggplot(#CAPRA_S values
  data.frame(x=c(-2,2)), #Shift CAPRA_S values of 0 - 7 by 3.5 to better simulate/understand behavouir
  #(TABI normalises CAPRA_S, such that values are centered around 0)
  aes(x=x)) + 
  stat_function(fun=peq, geom="line")	+ 
  labs(title="True Positive Curve", x="CAPRA-S", y="Gene Count")



#Use the eq above to simulate log_y_hat values 
log_y_hat<-sapply(seq(from = -2, to=2, by = 0.5),
              function(x)
                eq(x)  
)

#Simulated Negative Binomial Distributed Values
sample<-sapply(log_y_hat,
                       function(y) rnbinom(mu = exp(y), #Mean is exp of log y hat (GLA eq)
                                           size = 4, #Precision/overdispersion value of 0.5 seems approp. level of noise
                                           n = 10))

#Plot simulated data 

data.frame(
  CAPRA_S= c(#Now using "real" un-normalised CAPRA_S
    sapply(seq(from = -2, to=2, by = 0.5), 
           function(x) rep(x, 10))), 
  sample = c(sample)
)  %>% 
  ggplot() + 
  geom_point(aes(x=CAPRA_S, y=sample+1), col="dodgerblue")  +
  scale_y_log10() 	+ 
  labs(title="True Positive Curve", x="CAPRA-S", y="Gene Count")


Sig_big_subset %>% 
  filter(Gen_ref == "X1_Slope_-5_k_3") %>% 
  ggplot(aes(x=CAPRA_S, y=value)) + 
  geom_point()




#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


#Section 2 


#Data simulation Functions


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #


#Function for simulating Data with both null and true (sigmodial curve) tests



sigmodial_tests<-function(n_true_tests,
                          n_tests, 
                          sample_size) { #With the way this function is code Sample size needs to be multiple of 15
  
  require(reshape2)
  require(dplyr)
  
  #Calculation of x-coordinates ("simulation" CAPRA_S values)
  x_cord<-seq(from = -2, to=2, by = 0.5) #9 different values of CAPRA_S
  
  
  
  #         #          #          #         #           #
  
  
  
  #Creating the Sigmodial curve based genes (true tests - dif transcribed genes)
  
  
  #Simulation all possible values of the inflection (eta)
  set.seed(20)
  eta = sample(seq(from = -2, to = 2, by =0.5), 
               size = n_true_tests, replace = T)
  
  #As above
  #Create Generalised Sigmoid Function Equation
  #Simulates log ( mean ) values 
  eq <- function(x) {
  #Paramter values
  beta=1.4
  y_0=2
  A = 1
  
  top<- ((y_0-A)*(1+exp(eta*beta))) #Numerator of the GLA function
  bottom<-(1+exp(eta*beta-x*beta)) #Denominator of the GLA function
  
  return((A+top/bottom)) #Evaluate the function 
  }
  
  #Use the eq  above to simulate y_hat values 
  y_hat<-sapply(seq(from = -2, to = 2, by =0.5), #For better behavouir of simulated curve, centre values around 0 
  #Replacating normalisation process undertaken by TABI
                function(x)
                  eq(x)  
  )
  
  #For y_hat
  #Each column is a for a value of CAPRA_S 
  #Each row is a simulated curve
  
  #Simulated Negative Binomial Distributed Values
  sample<-apply(y_hat, 
                1, #For each row (a curve)
                function(x) 
                  sapply(x, #Simulate sample size/15 
                         function(y) rnbinom(mu = exp(y), 
                                                size = 0.5, #Precision/overdispersion values
                                                n = sample_size/9))
                
  )
  
  #In Sample
  #Each column is a simulated curve / gene 
  #Each row is a sample
  
  
  #Combine respective CAPRA_S values with their respective simulated gene counts
  True_tests<-data.frame(
    CAPRA_S= c(
      sapply(seq(from = -2, to=2, by = 0.5), 
             function(x) rep(x, sample_size/9))), #
             sample
  )
  
  #Reshape the above such that it mimics TCGA data set (column labeling gene count, column labeling gene name,
  #column labbeling CAPRA_S value)
  
  true_tests<-True_tests %>% 
    melt(id.vars = 1) %>% #Melt with respect to CAPRA_S 
    mutate(Null_test = "FALSE") %>% #Add a column which indicates this is a null 
    rename(Gene_number = variable) #Make Gene_number the test identification column
  
  
  
  #         #          #          #         #           #
  
  
  
  #Simulation of null tests (negbinom distributed, with varying means  - same dispersion)
  
  n_null<-n_tests-n_true_tests #Number of null tests
  
  #Simulation a n_null number of means - any integer between 0 and 15000
  null_means<-sample(0:5000, 
                     size = n_null, 
                     replace = T)
  
  #Simulate null_values 
  null_values<-sapply(null_means, 
                      function(x) 
                        (rnbinom(n=sample_size,
                                 mu = x, 
                                 size = 0.5)))
  
  null_tests<-data.frame(
    CAPRA_S= c(
      sapply(seq(from = -2, to=2, by = 0.5), 
             function(x) rep(x, sample_size/9))), #
    null_values
  ) %>% 
    melt(id.vars = 1) %>% #Melt with respect to CAPRA_S 
    mutate(Null_test = "TRUE") %>% #Add a column which indicates this is a null 
    rename(Gene_number = variable) 
  
  
  null_tests$Gene_number<-sub("X", "V", as.factor(null_tests$Gene_number)) #Convert X test recognition to V to further
  #Distinguish from true tests (now true tests have an X prefix, and null have a V prefix)

  #Return a data frame with CAPRA_S, Gene_number (Gene identifier colum), value (simulated gene expression value)
  #and Null_test (identifying true vs null tests)
  return(rbind(true_tests, 
               null_tests))
  
}


#Simulate a Mixture of null and true tests

sigmodial_test<-sigmodial_tests(100, 1000, 150) #with size of precision/dispersion size = 0.5, sample size= 150, 
#100 true associations, 1000 tests in total


#Double check some simulated curves seem reasonable by plotting 

library(ggplot2)

#Plotting 
#True associations 

sigmodial_test %>% 
  filter(Gene_number == "X40") %>% 
ggplot(aes(x=CAPRA_S, y=value+1)) +
  geom_point(col = "dodgerblue") +
  scale_y_log10() + 
  labs(title="True Positive Curve", x="CAPRA-S", y="Gene Count")

sigmodial_test %>% 
  filter(Gene_number == "X50") %>% 
  ggplot(aes(x=CAPRA_S, y=value+1)) +
  geom_point(col = "dodgerblue") +
  scale_y_log10() + 
  labs(title="True Positive Curve", x="CAPRA-S", y="Gene Count")

sigmodial_test %>% 
  filter(Gene_number == "X1") %>% 
  ggplot(aes(x=CAPRA_S, y=value+1)) +
  geom_point(col = "dodgerblue") +
  scale_y_log10() + 
  labs(title="True Positive Curve", x="CAPRA-S", y="Gene Count")

sigmodial_test %>% 
  filter(Gene_number == "X15") %>% 
  ggplot(aes(x=CAPRA_S, y=value+1)) +
  geom_point(col = "dodgerblue") +
  scale_y_log10() + 
  labs(title="True Positive Curve", x="CAPRA-S", y="Gene Count")

#Null associations check


sigmodial_test %>% 
  filter(Gene_number == "V40") %>% 
  ggplot(aes(x=CAPRA_S, y=value+1)) +
  geom_point(col = "dodgerblue") +
  scale_y_log10() + 
  labs(title="No underlying Curve", x="CAPRA-S", y="Gene Count")

sigmodial_test %>% 
  filter(Gene_number == "V500") %>% 
  ggplot(aes(x=CAPRA_S, y=value+1)) +
  geom_point(col = "dodgerblue") +
  scale_y_log10() + 
  labs(title="No underlying Curve", x="CAPRA-S", y="Gene Count")

sigmodial_test %>% 
  filter(Gene_number == "V1") %>% 
  ggplot(aes(x=CAPRA_S, y=value+1)) +
  geom_point(col = "dodgerblue") +
  scale_y_log10() + 
  labs(title="No underlying Curve", x="CAPRA-S", y="Gene Count")

sigmodial_test %>% 
  filter(Gene_number == "V666") %>% 
  ggplot(aes(x=CAPRA_S, y=value+1)) +
  geom_point(col = "dodgerblue") +
  scale_y_log10() + 
  labs(title="No underlying Curve", x="CAPRA-S", y="Gene Count")


#Save for later use (testing on TABI / other methods)

saveRDS(sigmodial_test,
        compress = "gzip",
        file = "sigmodial_test.rds")

#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

#Simulate all null results 

#Simulation a n_null number of means - any integer between 0 and 7000

n_null<-1000
sample_size<-150
CAPRA_S<-seq(from = 0, to=7, by = 0.5)

null_means<-sample(0:7000, 
                   size = n_null, 
                   replace = T)

#Simulate null_values 
null_values<-sapply(null_means, 
                    function(x) 
                      (rnbinom(n=sample_size,
                               mu = x, 
                               size = 4)))

null_result<-as.data.frame(
  cbind(null_values,
        CAPRA_S)) %>%
  as_tibble() %>%
  melt(id.vars = c(n_null+1)) %>%
  rename(Gene_number = variable, CAPRA_S = CAPRA_S)


saveRDS(null_result, 
        compress = "gzip",
        file = "null_result.rds")


#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

#Simulation all results as true

#Function only producing true results
only_sigmodial_tests<-function(n_true_tests, n_tests, sample_size) { #Sample size needs to be multiple of 15
  
  require(reshape2)
  require(dplyr)
  
  #Calculation of x-coordinates ("simulation" CAPRA_S values)
  x_cord<-seq(from = 0, to=7, by = 0.5) #15 different values of CAPRA_S
  
  #Sigmodial curve (true tests)
  #Simulation all possible values of the inflection (eta)
  eta = sample(c(0:7), size = n_true_tests, replace = T)
  
  #Create Generalised Sigmoid Function Equation
  #Simulates log ( mean ) values 
  eq <- function(x) {
    
    beta=0.6
    y_0=1.5
    A = 1
    
    top<- ((y_0-A)*(1+exp(eta*beta)))
    bottom<-(1+exp(eta*beta-x*beta))
    
    return((A+top/bottom))
  }
  
  #Use the eq  above to simulate y_hat values 
  y_hat<-sapply(x_cord,
                function(x)
                  eq(x)  
  )
  
  #For y_hat
  #Each column is a for a value of CAPRA_S 
  #Each row is a simulated curve
  
  #Simulated Negative Binomial Distributed Values
  sample<-apply(y_hat, 
                1, #For each row (a curve)
                function(x) 
                  sapply(x, #Simulate sample size/15 
                         function(y) rnbinom(mu = exp(y), 
                                             size = 0.5, #Precision/overdispersion values
                                             n = sample_size/15))
                
  )
  
  #In Sample
  #Each column is a simulated curve / gene 
  #Each row is a sample
  
  
  #Combine the 
  True_tests<-data.frame(
    CAPRA_S= c(
      sapply(seq(from = 0, to=7, by = 0.5), 
             function(x) rep(x, sample_size/15))), #
    sample
  )
  
  #Reshape the 
  true_tests<-True_tests %>% 
    melt(id.vars = 1) %>% #Melt with respect to CAPRA_S 
    mutate(Null_test = "FALSE") %>% #Add a column which indicates this is a null 
    rename(Gene_number = variable) #

  return(true_tests)
}



true_results<-only_sigmodial_tests(1000, 1000, 150)

saveRDS(true_results, 
        compress = "gzip",
        file = "true_results.rds")

#   #    #   #   #   #   #   #   #    #   #   #   #   #  #   #    #   #   #   #   #

#Mixture of null and true tests

sigmodial_test<-sigmodial_tests(100, 1000, 150) #with size of precision/dispersion size = 4

saveRDS(sigmodial_test,
        compress = "gzip",
        file = "sigmodial_test.rds")