![TABI](man/TABI_logo.png)

Transcriptomic Analyses through Bayesian Inference

# Test

```R
# Import  
if("package:TABI" %in% search()) detach("package:TABI", unload=TRUE, force=TRUE); library(TABI)  

# Simulate data  
d = simulate_from_sigmoid(n_samples = 13,precision_NB = 100)  

# Run model  
res = TABI_glm(
	~ cov,
	tibble( cov = d$X[,2] ) %>% bind_cols(as_tibble(d$y) ) %>% setNames(gsub("^V", "", colnames(.))),
	prior = list(prop_DE = 0.05 ,	scale_DE = 5) 
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
	prior = list( 
		prop_DE =0.1,
		scale_DE = 5
	),
	model=rstan::stan_model("~/PhD/TABI/src/stan_files/DE_sigmoid.stan")
)

ggplotly(plot_posterior(tabi_res, CI=0.95, covariate = "CAPRA_TOTAL"))


```
