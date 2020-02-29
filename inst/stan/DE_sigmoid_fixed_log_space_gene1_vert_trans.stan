functions{

  vector gla_eq(row_vector y_lin, row_vector b0, vector y_cross, vector A) {
    return  A + (y_cross -A) .* (1+exp(-to_vector(b0))) ./ (1+exp(- to_vector(y_lin) )) ;
  }
  
  vector gla_eq_2(row_vector y_lin, row_vector b0, vector y_cross, vector A) {
    return A + (y_cross -A) .*  exp(  log1p_exp(-to_vector(b0)) - log1p_exp(- to_vector(y_lin)) );
  }
}

data {
	int <lower=0, upper = 1> prior_only; // For testing purpose

	int<lower = 0> G;                   // All Genes
	int<lower = 0> T;                   // All Samples (e.g. tube/individual)
	int<lower=0> R_1;                   // All Covariates (Intercept, Chosen Covariate 1 ... )  
	int<lower = 0> y[T, G];             // RNA-seq counts
	matrix[T,R_1+1] X;                 // Design Matrix
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
}

parameters {
	// Linear model
	
	row_vector[G] inflection;
	vector<lower=0>[G] y_cross_raw;
	
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
	vector[G] log_y_hat[T];
	vector[G] phi[T];
  vector[G] y_cross = y_cross_raw + A;
  
	// Building matrix factors of interest
  beta[2] = to_row_vector(beta1_z[1]);
	if(R_1 > 1)	for(r in 2:R_1) beta[r+1] = to_row_vector( beta1_z[r] * non_sparse_sigma[r-1]);

	// Inflection point
	beta[1] = to_row_vector(-inflection .* beta[2]);

	// Calculation of generalised logit
	for(t in 1:T) log_y_hat[t] = gla_eq_2(X[t] * beta, beta[1], y_cross, A);
	
	//Calculation of Overdispersion 
	for(t in 1:T) phi[t] = -0.3186 * log_y_hat[t] + od;
	//phi[,1] = y_hat[,1] * -0.3186 + 0.7847; 
	
	
	//Print comparing gla_eq and gla_eq_2
	
	print(gla_eq([7.0], [12.0], [-1.0]', [1.0]'), " gla_eq(1,1,1,1) ");
	print(gla_eq_2([7.0], [12.0], [-1.0]', [1.0]'), " gla_eq_2(1,1,1,1) ");
	
	
	print(gla_eq([2.0], [-15.0], [21.0]', [16.0]'), " gla_eq(1,1,1,1) ");
	print(gla_eq_2([2.0], [-15.0], [21.0]', [16.0]'), " gla_eq_2(1,1,1,1) ");
	
	print(gla_eq([2.0], [54.0], [-3.0]', [1.0]'), " gla_eq(1,1,1,1) ");
	print(gla_eq_2([2.0], [54.0], [-3.0]', [1.0]'), " gla_eq_2(1,1,1,1) ");
	
	print(gla_eq([-4.0], [-2.0], [-3.0]', [-1.0]'), " gla_eq(1,1,1,1) ");
	print(gla_eq_2([-4.0], [-2.0], [-3.0]', [-1.0]'), " gla_eq_2(1,1,1,1) ");
}
model {

	// Linear system
	for(r in 1:R_1) beta1_z[r] ~ normal (0 , 2);
	inflection ~ normal(0 ,2);
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

	// Likelihood
	if(prior_only == 0) for(t in 1:T) y[t] ~ neg_binomial_2_log(log_y_hat[t], 1 ./ exp(phi[t]));


}

generated quantities{
  int y_gen[T,G];          // RNA-seq counts
		// vector[G] gamma_log_sampling;
	//int gamma_rate[T,G];   //For working problem cases

 //print(phi[1,1], " phi ", 
 //log_y_hat[1,1], " log_y_hat " , 
 //A[1], " A ", 
 //od, "od",
 //inflection[1], " inflection ", 
 //y_cross_raw[1], " y_cross_raw ", 
 //non_sparse_sigma, " non_sparse_sigma", 
 //beta1_z, " beta1_z");
 
 

 
	//Generate the data
	for (t in 1:T) {
	  for(g in 1:G) {
	  // If non-valid neg binomial values return parameter values instead
	  //if(is_nan(log_y_hat[t,g]))	
	  //print(phi[t], " phi ", log_y_hat[t], " log_y_hat " , A[g], " A ", inflection[g], " inflection ", y_cross_raw[g], " y_cross_raw ")
	  //else if( gamma_rng( 1 ./ exp(phi[t,g]),  1 ./ (exp(phi[t,g]) .* exp(log_y_hat[t,g]) ))>= exp(20.79)) //Gamma rate is greater than or equal to 
	  //print(phi[t], " phi ", log_y_hat[t], " log_y_hat" , A[g], " A ", inflection[g], "inflection", y_cross_raw[g], "y_cross_raw")
	  //else 
	  
	  y_gen[t,g] = neg_binomial_2_log_rng(log_y_hat[t,g], 1 ./ exp(phi[t,g])); 
	  } 
	  }
	  
	  
	

//*/
  //int neg_binomial_2_log_safe_rng(real eta, real phi) {
    //real gamma_rate = gamma_rng(phi, phi / exp(eta));
    //if (gamma_rate >= exp(20.79))
      //return -9;      
    //return poisson_rng(gamma_rate);
		// for(g in 1:G) gamma_log_sampling[g] = log( gamma_rng(log_y_cross_prior[1] + 1, inv(log_y_cross_prior[2] * 100)  ) );

	}