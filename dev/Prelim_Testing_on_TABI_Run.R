
#Preliminary Testing TABI - Unadjusted CAPRA_S

library(dplyr)
library(ggplot2)
library(tidyr)


#First 50000 Genes 

fit_100<-rbind_all(fit_100)

fit_p1<-rbind(fit_100, fit_5000)

# (1) Genes for which slope and inflection are different to zero

Sig_trends<- fit_p1 %>% 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  filter(`2.5%`<0&`97.5%`<0|`2.5%`>0&`97.5%`>0) 

# Sort by furthest from zero (slope, 2.5% lower margin)

#Infletion 
Sig_trends %>% 
  select(parameters, `2.5%`, Gene_name) %>%
  pivot_wider(names_from = parameters, values_from = `2.5%`) %>%
  rename(Inflection="inflection[1]", Slope ="beta[2,1]") %>%
  arrange(desc(abs(Infletion)))

#Slope
Sig_trends %>% 
  select(parameters, `2.5%`, Gene_name) %>%
  pivot_wider(names_from = parameters, values_from = `2.5%`) %>%
  rename(Inflection="inflection[1]", Slope ="beta[2,1]") %>%
  arrange(desc(abs(Infletion)))

# (2) Scatterplot of relationship between inflection and slope 

# Pick Unbiased Sample of Genes 

lm_test<-fit_p1 %>% #Data Table 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  select(parameters, mean, Gene_name) %>% 
  pivot_wider(names_from = parameters, values_from = mean)%>%
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]") 


fit_p1 %>% #Data Table 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  select(parameters, "mean", Gene_name) %>% 
  pivot_wider(names_from = parameters, values_from = "mean")%>%
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]") %>% 
  ggplot(aes(x=Inflection_mean, y= Slope_mean)) +
  geom_point() + stat_smooth(se = FALSE, method = 'lm', col="red")


(fit_p1 %>% #Data Table 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  select(parameters, "mean", Gene_name) %>% 
  pivot_wider(names_from = parameters, values_from = "mean")%>%
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]"))[sample.int(5000, 1500, replace=FALSE),] %>% #Take random sample of rows i.e. different genes
  ggplot(aes(x=Inflection_mean, y= Slope_mean)) +
  geom_point() + stat_smooth(se = FALSE, method = 'lm', col="red")

  
fit_p1 %>% 
  filter(Gene_name == c(sample))
class(c(sample))
summary(lm(Inflection_mean ~ Slope_mean, data = lm_test))

# (2.1) Look at genes for which up/low on graph

#Bottom Right Corner Genes,
#Inflection toward the right, negative slope

fit_100 %>% 
  select(parameters, Gene_name, mean) %>% 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  pivot_wider(names_from = parameters, values_from = mean) %>% 
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]") %>% 
  filter(Inflection_mean>1&(Slope_mean<1*-1))


fit_p1 %>% 
  select(parameters, Gene_name, mean) %>% 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  pivot_wider(names_from = parameters, values_from = mean) %>% 
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]") %>% 
  filter(Inflection_mean>2&(Slope_mean<1*-1.5))



normalised_PC_TCGA %>% 
  filter(transcript == c("AKT3")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

normalised_PC_TCGA %>% 
  filter(transcript == c("ARF1")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

normalised_PC_TCGA %>% 
  filter(transcript == c("ATP6V1E1")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()


#Fit Seems Reasonable

normalised_PC_TCGA %>% 
  filter(transcript == c("AARS")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +geom_point() +scale_y_log10()


normalised_PC_TCGA %>% 
  filter(transcript == c("ABCA5")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

#?Questionable 

normalised_PC_TCGA %>% 
  filter(transcript == c("ABCB8")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

fit_100 %>% 
  filter(Gene_name =="ABCB8")

#Top Left Corner Genes
#Inflection toward the left, positive slope

fit_p1 %>% 
  select(parameters, Gene_name, mean) %>% 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  pivot_wider(names_from = parameters, values_from = mean) %>% 
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]") %>% 
  filter(Inflection_mean<1*-2&(Slope_mean)>2)

normalised_PC_TCGA %>% 
  filter(transcript == c("BAP1")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

normalised_PC_TCGA %>% 
  filter(transcript == c("C12orf49")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

normalised_PC_TCGA %>% 
  filter(transcript == c("CALCOCO2")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()


#Reasonable Fit (first 100, <1.5, >2)

normalised_PC_TCGA %>% 
  filter(transcript == c("A2M")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()


normalised_PC_TCGA %>% 
  filter(transcript == c("ABCF3")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()


# Genes that fit in top left corner

fit_p1 %>% 
  select(parameters, Gene_name, mean) %>% 
  filter(parameters == "beta[2,1]"|parameters == "inflection[1]") %>% 
  pivot_wider(names_from = parameters, values_from = mean) %>% 
  rename(Inflection_mean="inflection[1]", Slope_mean ="beta[2,1]") %>% 
  filter(Inflection_mean>2&(Slope_mean)>0.5)

normalised_PC_TCGA %>% 
  filter(transcript == c("CALHM3")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()
  

# Genes that fit in top right corner


normalised_PC_TCGA %>% 
  filter(transcript == c("CDH20")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

normalised_PC_TCGA %>% 
  filter(transcript == c("CCL22")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

normalised_PC_TCGA %>% 
  filter(transcript == c("CD5L")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

normalised_PC_TCGA %>% 
  filter(transcript == c("CD69")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

normalised_PC_TCGA %>% 
  filter(transcript == c("CLDN4")) %>% 
  filter(is.na(CAPRA_S)==F, is.na(`read count normalised`)==F) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10()

# Look at Genes which would fit in top right (inflection toward the right - late change, positive slope)
# Or in bottom left (inflection toward the left - early change, negative slope)

#Plotting Function 
plot <-function(x) {
  normalised_PC_TCGA %>% 
  filter(transcript == levels(normalised_PC_TCGA$transcript)[x]) %>% 
  ggplot(aes(x=CAPRA_S, y=`read count normalised`)) +
  geom_point() +
  scale_y_log10() +
  labs(title =levels(normalised_PC_TCGA$transcript)[x])}

#Plot Random Genes
plot(sample.int(length(levels(normalised_PC_TCGA$transcript)),1))


#Which genes appeared to fit in top right / bottom left (to be tested in TABI)
normalised_PC_TCGA %>% 
  filter(transcript == c("PIGQ", "DNAI1", "SERPINE2", "SSTR3", "ISM2", "GPAT4", "PLD3", "TAOK3", "RAPGEF3", "MTND5P28", "NSL1", "STARD6", "CST7", "NDUFA1", "OSGEPL1-AS1", "LRRC17", "ERICD"))



# (2.2) Heatmap / Test on Simulated Data

#Inflection_mean > 2, Positive Slope > 2 


sim_eq<- function(x,eta, beta, y_0, A) {
  top<-(y_0-A)*(1+exp(eta*beta))
  bottom<-(1+exp(-x*beta+eta*beta))
  return(A + top/bottom)
}

CAPRA_S = seq(from=0, to =8, by=0.01)
sim_eq(CAPRA_S,3, 3, 2, 2)



# (3) Statistics and Histograms about Slopes and Inflections

# Slope Summary Statistics 
fit_100 %>% #Data Table
  filter(parameters == "beta[2,1]") %>% #Only Slope Information
  select(-parameters, -Gene_name) %>% 
  summary() 


# Histogram of Mean Slope
# Spread is larger than standard normal / maybe t 2 degrees freedom / slightly left centre
fit_100 %>%  #Data Table
  filter(parameters == "beta[2,1]") %>%
  select(parameters, "mean", Gene_name) %>% 
  ggplot() + geom_histogram(aes(x=mean))  #Mean name here, "X2.5"


# Inflection Summary Statistics
fit_100 %>% #Data Table
  filter(parameters == "inflection[1]") %>% #Only Slope Information
  select(-parameters, -Gene_name) %>% 
  summary() 


# Histogram of Mean Inflection 
# Appears to follow standard normal
fit_100 %>%  #Data Table
  filter(parameters == "inflection[1]") %>%
  select(parameters, "mean", Gene_name) %>%  #Mean name here, "X2.5"
  ggplot() +geom_histogram(aes(x=mean))  #Mean name here, "X2.5"
