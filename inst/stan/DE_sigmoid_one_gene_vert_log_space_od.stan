functions{

  vector gen_inv_logit(row_vector y_lin, row_vector b0, vector y_cross, vector A) {
    return  A + (y_cross -A) .* (1+exp(-to_vector(b0))) ./ (1+exp(- to_vector(y_lin) )) ;
  }

	real gamma_log_lpdf(vector x_log, real a, real b){

		vector[rows(x_log)] jacob = x_log; //jacobian
		real norm_constant = a * log(b) -lgamma(a);
		real a_minus_1 = a-1;
		return sum( jacob ) + norm_constant * rows(x_log) + sum(  x_log * a_minus_1 - exp(x_log) * b ) ;

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
	real < lower =0 > aux1_global = 2;
	real < lower =0 > aux2_global = 1;
	real < lower =0 > caux = 1;

	// Overdispersion of Dirichlet-multinomial
	// real<lower=0> od_inflection = 5.775609e+00;
	// real<lower=0> od1 = 4.659860e+00;
	// real<lower=0> od_k = 2.571095e-01;
}

parameters {
	// Linear model
	row_vector[G] inflection;
	vector<lower=0>[G] y_cross_raw;
	//real<lower=0> log_y_cross_prior[2];
	vector[G] beta1_z[R_1];
	//vector[T] normalization;

// Overdispersion
	real od;
	//Vertical Translatioon
	vector[G] A;
	

	// Non sparse sigma
	vector<lower=0>[R_1-1] non_sparse_sigma;

}

transformed parameters {

	matrix[R_1+1, G] beta;
	vector[G] y_hat[T];
	vector[G] phi[T];
  vector[G] y_cross = y_cross_raw + A;
  
	// Building matrix factors of interest
  beta[2] = to_row_vector(beta1_z[1]);
	if(R_1 > 1)	for(r in 2:R_1) beta[r+1] = to_row_vector( beta1_z[r] * non_sparse_sigma[r-1]);

	// Inflection point
	beta[1] = to_row_vector(-inflection .* beta[2]);

	// Calculation of generalised logit
	for(t in 1:T) y_hat[t] = gen_inv_logit(X[t] * beta, beta[1], y_cross, A);
	
	//Calculation of Overdispersion 
	for(t in 1:T) phi[t] = -0.3186 * y_hat[t] + od;
	//phi[,1] = y_hat[,1] * -0.3186 + 0.7847; 
}
model {

	// Linear system
	for(r in 1:R_1) beta1_z[r] ~ normal (0 , 2);
	inflection ~ normal(0 ,2);
	y_cross_raw ~ normal(0,1); //gamma_log(explog_y_cross_prior[1]) * inv(exp(log_y_cross_prior[2])), inv(exp(log_y_cross_prior[2])) );
	//log_y_cross_prior ~ normal(0,5);

	// normalization ~ normal(0,1);
	// sum(normalization) ~ normal(0, 0.001*T);
	
	//Vertical Translation
	A ~ normal(0,2);
	
	od ~ normal(0,2);

	// Non sparse sigma
	if(R_1 > 1) non_sparse_sigma ~ normal(0, 1);

	// Likelihood
	if(prior_only == 0) for(t in 1:T) y[t] ~ neg_binomial_2_log	(y_hat[t], 1 ./ exp(phi[t]));


}

generated quantities{
  int y_gen[T,G];          // RNA-seq counts
		// vector[G] gamma_log_sampling;

	// Generate the data
	for (t in 1:T) for(g in 1:G)
		y_gen[t,g] = neg_binomial_2_log_rng(   y_hat[t,g], 1 ./ exp(phi[t,g]));

		// for(g in 1:G) gamma_log_sampling[g] = log( gamma_rng(log_y_cross_prior[1] + 1, inv(log_y_cross_prior[2] * 100)  ) );

	}
