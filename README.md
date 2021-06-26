![TABI](man/TABI_logo.png)

Transcriptomic Analyses through Bayesian Inference

# Usage 

To apply TABI, read count data should be a table (tibble, data.frame) in tidy format (i.e. each row an observation, each column a variable). Additionally, the data should have: 

- a column which identifies the sample (biological replicate) 
- a column which identifies the transcript 
- a column of unnormalised read counts (*Importantly* as TABI is designed to work with unnormalised read counts, if you supply TABI with normalised read counts the model will behave unpredictablly)
- a column for each factor of interest 

# Example usage on a subset of 100 transcripts from the TCGA Dataset

We provide a example below of running TABI on a subset of transcripts from the TCGA prostate cancer data. In this example we have undertaken RNA sequencing on tumours from a cohort of patients with prostate cancer. We are interested in using TABI to measure how read count varies across a pseudo-continuous factor of interest, CAPRA-S (a clinical measurement of likelihood of prostate cancer reoccurrence).

```R

if("package:TABI" %in% search()) detach("package:TABI", unload=TRUE, force=TRUE)
library(TABI)


# Load sample TCGA data (is in data folder)

TCGA_example = readRDS("./data/example_TCGA_data.rds")


# Data should be in tidy format 

head(TCGA_example) # For example see format of TCGA



# We need to filter out NAs in the covariate before applying TABI

TCGA_example_na_rm = TCGA_example %>% dplyr::filter(!is.na(CAPRA_S))


# If we are interested in 

# Run TABI_glm 

TABI::TABI_glm(
 .data = TCGA_example_na_rm,
 .transcript = transcript, #The column which identifies the transcript 
 .sample = sample.ID, # The column which identifies the sample
 .abundance = read_count, # The column of read count 
	formula = ~ CAPRA_S+purity.score,
	
)


```

 We also found previously using PCA the purity score of the tumour sample , hence we 



# Test




```R
# Import  
if("package:TABI" %in% search()) detach("package:TABI", unload=TRUE, force=TRUE); library(TABI)  

# Simulate data  
d = simulate_from_sigmoid(n_samples = 13,precision_NB = 100)  

# Run model  
res = TABI_glm(
	~ cov,
	tibble( cov = d$X[,2] ) %>% bind_cols(as_tibble(d$y) ) %>% setNames(gsub("^V", "", colnames(.)))
)

# Analyse results  
source("https://raw.githubusercontent.com/betanalpha/knitr_case_studies/master/qr_regression/stan_utility.R")
check_all_diagnostics(res$fit)
plot_posterior(res, CI=0.95)
plot_generated_gene(res, "99")

```
# Test real data

```R

# Run TABI_glm
if("package:TABI" %in% search()) detach("package:TABI", unload=TRUE, force=TRUE)
library(TABI)

tabi_res = TABI_glm(
	formula = ~ CAPRA_TOTAL + batch,
	data = prostate_df, 
	prop_DE =0.1,
	scale_DE = 5,
	model=rstan::stan_model("~/PhD/TABI/src/stan_files/DE_sigmoid.stan")
)

ggplotly(plot_posterior(tabi_res, CI=0.95, covariate = "CAPRA_TOTAL"))


```

