functions{
	vector log_gen_inv_logit(row_vector y, row_vector k) {
		return to_vector( log(k) - log1p_exp(- ( y ) ) );
		}
}
data {
      int<lower = 0> F;                   // fixed genes
      int<lower = 0> G;                   // all genes
      int<lower = 0> T;                   // tube
			int<lower=0> R;
      int<lower = 0> y[T, G];             // RNA-seq counts
			matrix[T,R] X;

			int exposure[T];                 // How many reads have been sequenced for each sample
}
transformed data{
	matrix[R, F] zeros;
	zeros =  rep_matrix(0, R, F);
}
parameters {
	matrix[R, G-F] beta_changing;
	row_vector<lower=0>[G] k;

	real<lower=0> beta_sigma[R-1];

}
transformed parameters{
	matrix[R, G] beta = append_col(zeros,	beta_changing);
}
model {

	for (t in 1:T)
		y[t] ~ multinomial( softmax( log_gen_inv_logit(X[t] * beta, k) ) );

	beta_changing[1] ~normal(0,1);
	for(r in 2:R) beta_changing[r] ~ normal(0, beta_sigma[r-1]);

	k ~ normal(0,5);
	beta_sigma ~ normal(0,1);

}
generated quantities{
  int y_gen[T,G];          // RNA-seq counts

	for (t in 1:T)
			y_gen[t] = multinomial_rng( softmax( to_vector( log(k) - log1p_exp(- ( X[t] * beta ) ) ) ), exposure[t] );
}
