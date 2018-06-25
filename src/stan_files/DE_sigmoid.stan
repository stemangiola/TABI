functions{

  vector log_gen_inv_logit(row_vector y, vector inversion, vector intercept) {
    return  intercept + log1p_exp(-inversion) - log1p_exp(- to_vector(y)  ) ;
  }

  vector reg_horseshoe(
					vector zb,
					real aux1_global ,
					real aux2_global,
					vector  aux1_local ,
					vector aux2_local ,
					real  caux,
					real scale_global ,
					real slab_scale
					) {
    int K = rows(zb);

    // Horseshoe variables
		real tau ; // global shrinkage parameter
		vector [ K] lambda ; // local shrinkage parameter
		vector [ K] lambda_tilde ; // ’ truncated ’ local shrinkage parameter
		real c; // slab scale

		// Horseshoe calculation
		lambda = aux1_local .* sqrt ( aux2_local );
		tau = aux1_global * sqrt ( aux2_global ) * scale_global * 1 ;
		c = slab_scale * sqrt ( caux );
		lambda_tilde = sqrt ( c ^2 * square ( lambda ) ./ (c ^2 + tau ^2* square ( lambda )) );
		return  zb .* lambda_tilde * tau ;
  }
}

data {
	int <lower=0, upper = 1> prior_only; // For testing purpose
	int<lower = 0> G;                   // all genes
	int<lower = 0> T;                   // tube
	int<lower=0> R_1;
	int<lower = 0> y[T, G];             // RNA-seq counts
	matrix[T,R_1+1] X;
	int exposure[T];                 // How many reads have been sequenced for each sample

	// Horseshoe
	real < lower =0 > par_ratio ; // proportion of 0s
	real < lower =1 > nu_global ; // degrees of freedom for the half -t prior
	real < lower =1 > nu_local ; // degrees of freedom for the half - t priors
	real < lower =0 > slab_scale ; // slab scale for the regularized horseshoe
	real < lower =0 > slab_df; // slab degrees of freedom for the regularized

}

transformed data{
	real < lower =0 > scale_global = par_ratio / sqrt(1.0 * T); // scale for the half -t prior for tau
}

parameters {
	// Linear model
	vector[G] inversion;
	vector[G] intercept;
	real intercept_mu;
	real<lower=0> intercept_sigma;
	vector[G] beta1_z[R_1];
	vector[G] normalization;

	// Overdispersion of Dirichlet-multinomial
	real overdispersion;

	// Horseshoe
	real < lower =0 > aux1_global ;
	real < lower =0 > aux2_global ;
	vector < lower =0 >[ G] aux1_local ;
	vector < lower =0 >[ G] aux2_local ;
	real < lower =0 > caux ;

	// Non sparse sigma
	vector<lower=0>[R_1-1] non_sparse_sigma;

}

transformed parameters {

	vector[G] beta1[R_1];
	matrix[R_1+1, G] beta;
	matrix[T, G] X_beta;
	vector[G] y_hat[T];
	real<lower=0> exp_overdispersion;

	// Horseshoe calculation
	beta1[1] =
		reg_horseshoe(
			beta1_z[1],
			aux1_global ,
			aux2_global,
			aux1_local ,
			aux2_local ,
			caux,
			scale_global,
			slab_scale
		);

	// Other non sparse priors
	if(R_1 > 1)	for(r in 2:R_1)
		beta1[r] = beta1_z[r] * non_sparse_sigma[r-1];

	// make beta
	beta[1] = to_row_vector(inversion);
	for(r in 1:R_1) beta[r+1] = to_row_vector(beta1[r]);

	// Matrix multiplication for speed up
	X_beta = X * beta;
	for(t in 1:T) y_hat[t] = log_gen_inv_logit(X_beta[t], inversion, intercept) ;

	// Overdispersion
	exp_overdispersion = exp(overdispersion);

}
model {

	// Linear system
	for(r in 1:R_1) beta1_z[r] ~ normal (0 , 1);
	inversion ~ normal(0 ,1);
	intercept ~ normal(intercept_mu, intercept_sigma);
	intercept_mu ~ normal(0,1);
	intercept_sigma ~ cauchy(0, 2.5);
	normalization ~ normal(0,1);
	sum(normalization) ~ normal(0, 0.01*T);

	// Horseshoe
	aux1_local ~ normal (0 , 1);
	aux2_local ~ inv_gamma (0.5* nu_local , 0.5* nu_local );
	aux1_global ~ normal (0 , 1);
	aux2_global ~ inv_gamma (0.5* nu_global , 0.5* nu_global );
	caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df );

	// Non sparse sigma
	if(R_1 > 1) non_sparse_sigma ~ normal(0, 1);

	// Overdispersion
	overdispersion ~ normal(0, 1);

	// Likelihood
	if(prior_only == 0) for(t in 1:T) y[t] ~ neg_binomial_2_log	( log(exposure[t]) + normalization[t] + y_hat[t],  rep_vector(exp_overdispersion, G));

}

generated quantities{
  int y_gen[T,G];          // RNA-seq counts

	// Generate the data
	for (t in 1:T) for(g in 1:G)
		y_gen[t,g] = neg_binomial_2_log_rng( log(exposure[t]) + normalization[t] + y_hat[t,g],  exp_overdispersion );


	}
