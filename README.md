![TABI](man/TABI_logo.png)

Transcriptomic Analyses through Bayesian Inference

# Test

if("package:TABI" %in% search()) detach("package:TABI", unload=TRUE, force=TRUE)
library(TABI)

d = simulate_from_sigmoid(n_samples = 13,precision_NB = 100)

res = TABI_glm(
	~ cov,
	tibble( cov = d$X[,2] ) %>% bind_cols(as_tibble(d$y) ) %>% setNames(gsub("^V", "", colnames(.))),
	prior = list(prop_DE = 0.05 ,	scale_DE = 5) 
)

source("https://raw.githubusercontent.com/betanalpha/knitr_case_studies/master/qr_regression/stan_utility.R")
check_all_diagnostics(res$fit)

plot_posterior(res, CI=0.95)
plot_generated_gene(res, "99")
