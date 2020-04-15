library(tidyverse)
library(TABI)
library(tidybulk)

TABI_TP <- 
  TABI::test_df %>%
  mutate(count = rnorm(n(), 50, 10) %>% as.integer) %>%
  TABI_glm(
    ~ CAPRA_S,
    sample, transcript, count,
    model = rstan::stan_model("inst/stan/DE_sigmoid_hierarchical.stan"),
    control=list(
      adapt_delta=0.9,
      stepsize = 0.01,
      max_treedepth =10
    ),iter = 2000, warmup = 1000,
  )

# Censoring
TABI_TP <- 
  TABI::test_df %>%
  mutate(alive = sample(0:1, size = n(), replace = T)) %>%
  mutate(months = rgamma(n(), shape = 0.884, rate = 0.884 * exp(-4.72))) %>%
  mutate(count = rnorm(n(), 50, 10) %>% as.integer) %>%
  TABI_glm(
    ~ censored(months, alive),
    sample, transcript, count,
    model = rstan::stan_model("inst/stan/DE_sigmoid_hierarchical.stan"),
    control=list(
      adapt_delta=0.9,
      stepsize = 0.01,
      max_treedepth =10
    ),iter = 2000, warmup = 1000,
  )
mcmc_parcoord(as.array(TABI_TP$fit, pars = c("beta", "inflection", "y_cross_raw", "od", "A", "unseen")), np = nuts_params(TABI_TP$fit)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

gglines = 
  TABI_TP$fit %>%
  tidybayes::spread_draws(inflection[G], y_cross_raw[G], beta1_z[G], A[G]) %>%
  ungroup() %>%
  nest(data = .draw) %>%
  mutate(line = map(data, ~ stat_function(fun=eq, geom="line", args=c(eta=.x$inflection,   beta=.x$beta1_z,  y_0=.x$y_cross_raw,  A=.x$A)) )) %>%
  pull(line)

eq <- function(x, y_0, beta, eta, A) {
  # eta=4
  # beta=2
  # y_0=20
  # A=5000
  top<- ((y_0+A)*(1+exp(eta*beta)))
  bottom<-(1+exp(eta*beta-x*beta))
  return(exp(top/bottom + A)) 
  }


#Plot CUrve

test_df %>%
  ggplot(aes(x=scale(CAPRA_S), y=count)) +
  geom_point() +
  stat_function(fun=eq, geom="line", args=c(eta=2.60,   beta=-7.50,  y_0=0.564,  A=4.36))
  #stat_function(fun=eq, geom="line", args=c(eta=-0.860,   beta=8.31,  y_0=1.29,  A=4.13))
  
ggplot(data.frame(x=c(-2, 2)), aes(x=x)) + 
  labs(title="True Positive Curve", x="CAPRA-S", y="Normalised Gene Count")

