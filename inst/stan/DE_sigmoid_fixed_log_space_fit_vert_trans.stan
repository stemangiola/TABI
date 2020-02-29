functions{
  //Generalised sigmoid function for fitting to the mean of y in log space (i.e. to log_y_hat)
  vector gla_eq_2(row_vector y_lin, row_vector neg_eta_beta_1, vector y_cross, vector A) {
    // From paper - 
    //y_lin is linear equation in gla_eq (eta*beta - X*beta),
    //neg_eta_beta_1 is -eta*beta_1, 
    //y_cross is y_0 
    return A + (y_cross -A) .*  exp(  log1p_exp(-to_vector(neg_eta_beta_1)) - log1p_exp(- to_vector(y_lin)) );
  }
}

data {
	int <lower=0, upper = 1> prior_only; // For testing purpose (if ==0 then use data to generate quantities)

	int<lower = 0> G;                   // All Genes
	int<lower = 0> T;                   // All Samples (e.g. tube/individual)
	int<lower=0> R_1;                   // All Covariates (Intercept, Chosen Covariate 1 ... ) 
	int<lower = 0> y[T, G];             // RNA-seq counts
	matrix[T,R_1+1] X;                 // Design Matrix
	int exposure[T];                 // How many reads have been sequenced for each sample
	int<lower=0> TMM[T,G];           // TMM algorithim factor

	// Horseshoe (for multiple genes)
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
}

parameters {
	// Linear model
	
	row_vector[G] inflection; //Value of the inflection point on the x axis
	vector<lower=0>[G] y_cross_raw; 
	
	vector[G] beta1_z[R_1]; //Vector of slopes 
	
	//vector[T] normalization;

  // Overdispersion
	real od;
	
	//Vertical Translatioon
	vector[G] A;
	

	// Non sparse sigma
	vector<lower=0>[R_1-1] non_sparse_sigma;

}

transformed parameters {

	matrix[R_1+1, G] beta; //matrix of coefficents
	vector[G] log_y_hat[T];  //log of the mean of y
	vector[G] phi[T]; //log of the precision paramter - i.e dispersion in neg binomial is 1/exp(phi)
  vector[G] y_cross = y_cross_raw + A; //Restricted/defined y_cross to prevent problems with alterating signs in y_0 and A giving same result
  //hence preventing cases of multiple solutions
  
	// Building matrix factors of interest
  beta[2] = to_row_vector(beta1_z[1]);
	if(R_1 > 1)	for(r in 2:R_1) beta[r+1] = to_row_vector( beta1_z[r] * non_sparse_sigma[r-1]);

	// Eta*beta 
	beta[1] = to_row_vector(-inflection .* beta[2]);

	// Calculation of generalised logit - fitting in log space (i.e. log of the means follows gla eq)
	for(t in 1:T) log_y_hat[t] = gla_eq_2(X[t] * beta + beta[1], beta[1], y_cross, A);
	
	//old eq for(t in 1:T) log_y_hat[t] = gla_eq_2(X[t] * beta, beta[1], y_cross, A);
	
	//Calculation of Overdispersion 
	for(t in 1:T) phi[t] = -0.3186 * log_y_hat[t] + od;
	
}
model {

	// Linear system
	//Restricted priors on beta1_z[r], and inflection (were originally n(0,2)), preventing larger generated values
	//As these dramatically increase log_y_hat - which causes problems with neg_binomial_2_log / neg_binomial_2_log_rng
	for(r in 1:R_1) beta1_z[r] ~ normal (0,1);
	
	inflection ~ normal(0,1);
	
	y_cross_raw ~ normal(0,1); 
	
	//gamma_log(explog_y_cross_prior[1]) * inv(exp(log_y_cross_prior[2])), inv(exp(log_y_cross_prior[2])) );
	//log_y_cross_prior ~ normal(0,5);

	// normalization ~ normal(0,1);
	// sum(normalization) ~ normal(0, 0.001*T);
	
	//Vertical Translation
	A ~ normal(0,2);
	
	//overdispersion 
	od ~ normal(0,1);

	// Non sparse sigma
	if(R_1 > 1) non_sparse_sigma ~ normal(0, 1);
	

	// Likelihood - fiting 
	if(prior_only == 0) for(t in 1:T) y[t] ~ neg_binomial_2_log(to_vector(TMM[t]).*log_y_hat[t], 1 ./ exp(phi[t]));

}

generated quantities{
  int y_gen[T,G];          // RNA-seq counts

	//Generate the data 
	for (t in 1:T) {
	  for(g in 1:G) {
	  y_gen[t,g] = neg_binomial_2_log_rng(log_y_hat[t,g], 1 ./ exp(phi[t,g])); 
	  } 
	  }
	  
	}
