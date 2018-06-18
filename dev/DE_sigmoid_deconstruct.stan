functions {
 	real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
  	real alpha_plus = sum(alpha);

    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }

}

data {
	int<lower = 0> G;                   // all genes
	int<lower = 0> T;                   // tube

  // Indicator for hyperprior type 
  // 0 - no hyperprior
  // 1 - normal horseshoe	
  int<lower=0, upper=1> hyperprior_type; 
  
  // Indicator for link function type
  // 0 - linear (allowed only with normal and lognormal likelihood), 
  // 1 - exp
  // 2 - Log. Sigmoid (inv. logit), plateau parametrized via value at x=0 (inversion).
  // 3 - Sigmoid (inv. logit) reparametrization (TODO describe)
	int<lower=0,upper=2> link_type; 
	
	// Indicator for likelihood type
	// 0 - normal (for basic tests)
	// 1 - lognormal (for other tests)
	// 2 - negative binomial
	// 3 - dirichelt_multinomial + softmax
	int<lower=0,upper=3> likelihood_type; 
	
	//NOTE: Using different data/params for different likelihood/link types is done following
	//http://www.martinmodrak.cz/2018/04/24/optional-parameters/data-in-stan/
	
	//int<lower = 0> y[T, G];             // RNA-seq counts

  // Gene counts (for integer likelihoods only)	
	int y[likelihood_type > 1 ? T : 0, likelihood_type > 1 ? G : 0]; 
	
	// stand-in for counts for normal/lognormal likelihood only
	matrix[likelihood_type <= 1 ? T : 0, likelihood_type <= 1 ? G : 0] y_real; 
	
	vector[T] X;    
	
	real intercept_mean;
	real<lower=0> intercept_sigma;
	
	real<lower=0> beta_sigma[hyperprior_type == 0 ? 1 : 0];
	real<lower=0> horseshoe_tau[hyperprior_type == 1 ? 1 : 0];
	
  real inversion_mean[link_type == 2 ? 1 : 0];
  real<lower=0> inversion_sigma[link_type == 2 ? 1 : 0];
	
	//Total number of sequenced samples for multinomial likelihood
	int exposure[likelihood_type == 3 ? T : 0];
	
	real<lower=0> xi_sigma;
}

transformed data {
  matrix[T, 2] X_with_intercept;
  
  if(link_type == 0 && likelihood_type > 1) {
    reject("Linear link is allowed only for normal and lognormal likelihoods");
  }
  
  X_with_intercept[,1] = rep_vector(1, T);
  X_with_intercept[,2] = X;
}

parameters {
  vector[G] beta1_z;
	vector[G] intercept_z;
	
	real<lower=0> horseshoe_sigma_z[hyperprior_type == 1 ? 1 : 0];

  vector[link_type == 2 ? G : 0] inversion_z;
	
	real<lower=0> xi_z;
}

transformed parameters {

	real<lower=0> horseshoe_sigma[hyperprior_type == 1 ? 1 : 0];
  vector<lower=0>[link_type == 2 ? G : 0] inversion;
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
  } else {
    reject("Unknown hyperprior");
  }

	//------------ Latent linear model ----
	intercept = intercept_z * intercept_sigma + intercept_mean;
	if(link_type == 3) {
	  //Interpret intercept as the point where X_beta = 0
	  intercept = -intercept .* beta1;
	}
	
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
	} else if(link_type == 2) {
	  //Log Sigmoid via value at x=0 (inversion) point
    inversion = inversion_z * inversion_sigma[1] + inversion_mean[1];
	  for(g in 1:G) {
	    //iterating over genes should improve performance due to memory locality
	    //Original code:
	    //y_hat[t,] = log(to_row_vector(inversion)) + log1p_exp(to_row_vector(-intercept)) - log1p_exp(-X_beta[t,]);
	    //New code:
	    //TODO: there probably should be log(inversion[g]) and inversion constrained to positive
	    y_hat[,g] = inversion[g] + log1p_exp(-intercept[g]) - log1p_exp(-X_beta[,g]);

	  }
  } else {
    reject("Unknown link");
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
  } else {
    reject("Unknown hyperprior");
  }

  xi_z ~ normal(0, 1);

  //--------- Link function params ----
  inversion_z ~ normal(0, 1);

	//---------Likelihood----------------
	if(likelihood_type == 0) {
	  to_vector(y_real) ~ normal(to_vector(y_hat), xi);
	} else if(likelihood_type == 1) {
	  to_vector(y_real) ~ lognormal(to_vector(y_hat), xi);
	} else if(likelihood_type == 2) {
	  to_array_1d(y) ~ neg_binomial_2(to_array_1d(y_hat), xi);
	} else if(likelihood_type == 3) {
	  for(t in 1:T) {
	    y[t] ~ dirichlet_multinomial( xi * softmax(to_vector(y_hat[t])));
	  }
	} 
	else {
	  reject("Unknown likelihood")
	}
}

