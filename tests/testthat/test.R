library(tidyverse)
library(TABI)
library(tidybulk)

TABI_TP <- 
  test_df %>%
  mutate(count = rnbinom(n(), mu=50, size = 30) %>% as.integer) %>%
  TABI_glm(
    ~ CAPRA_S,
    sample, transcript, count,
    model = rstan::stan_model("inst/stan/DE_sigmoid_hierarchical.stan"),
    control=list(
      adapt_delta=0.9,
      stepsize = 0.01,
      max_treedepth =10
    ),
    iter = 2000,
    warmup = 1000,
  )

TABI_TP$fit %>% pairs(pars=c("beta", "inflection", "A", "od", "y_cross"))

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

