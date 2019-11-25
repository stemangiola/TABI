functions{

  vector log_gen_inv_logit(row_vector y_log, row_vector b0, vector log_y_cross) {
    return  log_y_cross + log1p_exp(-to_vector(b0)) - log1p_exp(- to_vector(y_log)  ) ;
  }

	real gamma_log_lpdf(vector x_log, real a, real b){

		vector[rows(x_log)] jacob = x_log; //jacobian
		real norm_constant = a * log(b) -lgamma(a);
		real a_minus_1 = a-1;
		return sum( jacob ) + norm_constant * rows(x_log) + sum(  x_log * a_minus_1 - exp(x_log) * b ) ;

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

	vector gen_inv_logit_overdispersion(vector y, real k) {
			// It uses the empirical mean of x axes instead of 0
	    return k ./ (1 + exp(-y));
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
	// real < lower =0 > aux1_global = 2;
	// real < lower =0 > aux2_global = 1;
	// real < lower =0 > caux = 1;

	// Overdispersion of Dirichlet-multinomial
	// real<lower=0> od_inflection = 5.775609e+00;
	// real<lower=0> od1 = 4.659860e+00;
	// real<lower=0> od_k = 2.571095e-01;
}

parameters {
	// Linear model
	row_vector[G] inflection;
	vector[G] log_y_cross;
	real<lower=0> log_y_cross_prior[2];
	vector[G] beta1_z[R_1];
	vector[T] normalization;

	// // Overdispersion of Dirichlet-multinomial
	real<lower=0> od_inflection;
	real<lower=0> od1;
	real<lower=0> od_k;

	// Horseshoe
	vector < lower =0 >[ G] aux1_local ;
	vector < lower =0 >[ G] aux2_local ;
	real < lower =0 > aux1_global;
	real < lower =0 > aux2_global;
	real < lower =0 > caux;

	// Non sparse sigma
	vector<lower=0>[R_1-1] non_sparse_sigma;

}

transformed parameters {

	matrix[R_1+1, G] beta;
	vector[G] y_hat[T];
	vector[G] overdispersion[T];
	real od0 = -od_inflection * od1;

	// Building matrix factors of interest
	beta[2] = to_row_vector(
		reg_horseshoe(
			beta1_z[1],
			aux1_global ,
			aux2_global,
			aux1_local ,
			aux2_local ,
			caux,
			scale_global,
			slab_scale
		)
	);
		if(R_1 > 1)	for(r in 2:R_1) beta[r+1] = to_row_vector( beta1_z[r] * non_sparse_sigma[r-1]);

	# Inflection point
	beta[1] = to_row_vector(-inflection .* beta[2]);

	// Calculation of generalised logit
	for(t in 1:T) y_hat[t] = log_gen_inv_logit(X[t] * beta, beta[1], log_y_cross) ;

	// Overdispersion for negative binomial
	for(t in 1:T) overdispersion[t] = 1 + gen_inv_logit_overdispersion(od0 + od1 * y_hat[t],  inv(od_k) ) ;

}
model {

	// Linear system
	for(r in 1:R_1) beta1_z[r] ~ normal (0 , 1);
	inflection ~ normal(0 ,5);
	log_y_cross ~ gamma_log(exp(log_y_cross_prior[1]) * inv(exp(log_y_cross_prior[2])), inv(exp(log_y_cross_prior[2])) );
	log_y_cross_prior ~ normal(0,5);

	normalization ~ normal(0,1);
	sum(normalization) ~ normal(0, 0.01*T);

	// Horseshoe
	aux1_local ~ normal (0 , 1);
	aux2_local ~ inv_gamma (0.5* nu_local , 0.5* nu_local );
	// aux1_global ~ normal (0 , 1);
	// aux2_global ~ inv_gamma (0.5* nu_global , 0.5* nu_global );
	// caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df );

	// Non sparse sigma
	if(R_1 > 1) non_sparse_sigma ~ normal(0, 1);

	// overdispersion
	od_inflection ~ normal(0,10);
	od1 ~ normal(0,1);
	od_k ~ normal(0,1);

	// Likelihood
	if(prior_only == 0) for(t in 1:T) y[t] ~ neg_binomial_2_log	(  normalization[t] + y_hat[t],  overdispersion[t]);

}

generated quantities{
  int y_gen[T,G];          // RNA-seq counts
		// vector[G] gamma_log_sampling;

	// Generate the data
	for (t in 1:T) for(g in 1:G)
		y_gen[t,g] = neg_binomial_2_log_rng(  normalization[t] + y_hat[t,g],  overdispersion[t,g] );

		// for(g in 1:G) gamma_log_sampling[g] = log( gamma_rng(log_y_cross_prior[1] + 1, inv(log_y_cross_prior[2] * 100)  ) );

	}
