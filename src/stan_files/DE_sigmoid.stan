functions{
	vector log_gen_inv_logit(row_vector y, row_vector k) {
		return to_vector( k - log1p_exp(- ( y ) ) );
		}
}
data {
      int<lower = 0> G;                   // all genes
      int<lower = 0> T;                   // tube
			int<lower=0> R;
      int<lower = 0> y[T, G];             // RNA-seq counts
			matrix[T,R] X;

			int exposure[T];                 // How many reads have been sequenced for each sample
}
parameters {
	matrix[R, G] beta;
	row_vector[G] k;

	//real<lower=0> beta_sigma[R-1];

}
model {

	for (t in 1:T)
		y[t] ~ multinomial( softmax( log_gen_inv_logit(X[t] * beta, k) ) );

	beta[1] ~normal(0,1);
	for(r in 2:R) beta[r] ~ double_exponential(0, 0.2);

	k ~ normal(0,1);
	//beta_sigma ~ normal(0,1);

	sum(k) ~ normal(0, 0.01 * G);

}
generated quantities{
  int y_gen[T,G];          // RNA-seq counts

	for (t in 1:T)
			y_gen[t] = multinomial_rng( softmax( log_gen_inv_logit(X[t] * beta, k)  ), exposure[t] );
}
