data {
	int<lower = 0> G;                   // all genes
	int<lower = 0> T;                   // tube

  int<lower=0, upper=1> hyperprior_type; //0 - no hyperprior, 1 - normal horseshoe	
	int<lower=0,upper=1> link_type; //0 - linear (only with normal likelihood), 1 - exp
	int<lower=0,upper=0> likelihood_type; //0 - normal
	
	//NOTE: Using different data/params for different likelihood/link types is done following
	//http://www.martinmodrak.cz/2018/04/24/optional-parameters/data-in-stan/
	
	//int<lower = 0> y[T, G];             // RNA-seq counts
	
	// stand-in for counts for normal likelihood only
	matrix[likelihood_type == 0 ? T : 0, likelihood_type == 0 ? G : 0] y_real; 
	vector[T] X;    
	
	real intercept_mean;
	real<lower=0> intercept_sigma;
	
	real<lower=0> beta_sigma[hyperprior_type == 0 ? 1 : 0];
	real<lower=0> horseshoe_tau[hyperprior_type == 1 ? 1 : 0];
	
	real<lower=0> xi_sigma;
}

transformed data {
  matrix[T, 2] X_with_intercept;
  
  if(link_type == 0 && likelihood_type != 0) {
    reject("Linear link is allowed only for normal likelihood");
  }
  
  X_with_intercept[,1] = rep_vector(1, T);
  X_with_intercept[,2] = X;
}

parameters {
  vector[G] beta1_z;
	// Linear model
	vector[G] intercept_z;
	
	real<lower=0> horseshoe_sigma_z[hyperprior_type == 1 ? 1 : 0];
	
	real<lower=0> xi_z;
}

transformed parameters {

	real<lower=0> horseshoe_sigma[hyperprior_type == 1 ? 1 : 0];
	vector[G] intercept;
  vector[G] beta1;
	matrix[2, G] beta;
	matrix[T, G] X_beta;
	matrix[T, G] y_hat;
	real<lower=0> xi;

  //------------ Hyperprior -------------
  if(hyperprior_type == 0) {
    beta1 = beta1_z * beta_sigma[1];
  } else if(hyperprior_type == 1) {
    horseshoe_sigma[1] = horseshoe_sigma_z[1] * horseshoe_tau[1];
    beta1 = beta1_z * horseshoe_sigma[1];
  }

	//------------ Latent linear model ----
	intercept = intercept_z * intercept_sigma + intercept_mean;
	beta[1] = to_row_vector(intercept);
	beta[2] = to_row_vector(beta1);

	// Matrix multiplication for speed up
	X_beta = X_with_intercept * beta;
	
	//------------ Link function ----------
	if(link_type == 0) {
	  //linear link
	  y_hat = X_beta;
	} else if(link_type == 1) {
	  //exp link
	  y_hat = exp(X_beta);
	}

  //----------- Likelihood noise -------
  xi = xi_z * xi_sigma;

}

model {

	//---------Priors for linear model---
	beta1_z ~ normal (0, 1);
	intercept_z ~ normal(0, 1);

  //---------Hyperprior----------------
  if(hyperprior_type == 0) {
    //No hyperprior
  }
  else if(hyperprior_type == 1) {
    //Normal horseshow
    horseshoe_sigma_z ~ normal(0, 1);
  }

  xi_z ~ normal(0, 1);

	//---------Likelihood----------------
	to_vector(y_real) ~ normal(to_vector(y_hat), xi);
}

