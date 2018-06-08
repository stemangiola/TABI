functions{
	vector log_gen_inv_logit(row_vector y, row_vector k) {
		return to_vector( k - log1p_exp(- ( y ) ) );
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
    int<lower = 0> G;                   // all genes
    int<lower = 0> T;                   // tube
		int<lower=0> R;
    int<lower = 0> y[T, G];             // RNA-seq counts
		matrix[T,R] X;
		int exposure[T];                 // How many reads have been sequenced for each sample

		// Horseshoe
		real < lower =0 > par_ratio ; // proportion of 0s
		real < lower =1 > nu_global ; // degrees of freedom for the half -t prior
		real < lower =1 > nu_local ; // degrees of freedom for the half - t priors
		real < lower =0 > slab_scale ; // slab scale for the regularized horseshoe
		real < lower =0 > slab_df ; // slab degrees of freedom for the regularized
	}

transformed data{
		real < lower =0 > scale_global = par_ratio / sqrt(1.0 * T); // scale for the half -t prior for tau
}

parameters {

		// Linear model
		row_vector[G] beta0;
		row_vector[G] k;
		vector<lower=0>[G] sigma_trick;

		// Horseshoe
		vector [ G] beta1_z;
		real < lower =0 > aux1_global ;
		real < lower =0 > aux2_global ;
		vector < lower =0 >[ G] aux1_local ;
		vector < lower =0 >[ G] aux2_local ;
		real < lower =0 > caux ;
}

transformed parameters {
	row_vector[G] beta1 = to_row_vector(reg_horseshoe(beta1_z, aux1_global , aux2_global, aux1_local , aux2_local , caux, 		scale_global, slab_scale ));
	matrix[R, G] beta = append_row(beta0, beta1);
}

model {

		// Horseshoe
		beta1_z ~ normal (0 , 1);
		aux1_local ~ normal (0 , 1);
		aux2_local ~ inv_gamma (0.5* nu_local , 0.5* nu_local );
		aux1_global ~ normal (0 , 1);
		aux2_global ~ inv_gamma (0.5* nu_global , 0.5* nu_global );
		caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df );

	// Lnear system
	beta0 ~ normal(0,1);

	// trick
	for(g in 1:G) beta[,g] ~ normal(0,sigma_trick[g]);
	sigma_trick ~ normal(0,.2);

	k ~ normal(0,2);
	sum(k) ~ normal(0, 0.01 * G);

	for (t in 1:T) y[t] ~ multinomial( softmax( log_gen_inv_logit(X[t] * beta, k) ) );
}

generated quantities{
  int y_gen[T,G];          // RNA-seq counts

	for (t in 1:T)
			y_gen[t] = multinomial_rng( softmax( log_gen_inv_logit(X[t] * beta, k)  ), exposure[t] );
}
