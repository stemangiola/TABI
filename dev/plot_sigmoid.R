log1p_exp = function(x){ log(1+exp(x)) }




plot(-10:10, sigmoid_4_param(
  matrix(-10:10, ncol = 1),
  A = 2,
  y_cross = 1,
  slope = matrix(1, nrow = 1),
  inflection= 1
  )
)




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

