functions{
  //Generalised sigmoid function for fitting to the mean of y in log space (i.e. to log_y_hat)
  row_vector gla_eq_2(row_vector y_lin, row_vector neg_eta_beta_1, vector y_cross, vector A) {
    // From paper - 
    //y_lin is linear equation in gla_eq (eta*beta - X*beta),
    //neg_eta_beta_1 is -eta*beta_1, 
    //y_cross is y_0 
    return to_row_vector(A + ((y_cross) .*  exp(  log1p_exp(-to_vector(neg_eta_beta_1)) - log1p_exp(- to_vector(y_lin)) )));
  }
  
  real gla_eq(row_vector x, real inflection, vector slope, real y_cross, real A) {

  real inf_slop = inflection * slope[1];
  return(
    
       A + 
      (
        (y_cross) * 
        exp(   
          log1p_exp(-inf_slop) -
          log1p_exp(   -( x * slope ) - inf_slop  ) 
          
        )
      ) 
    
  );
    
  }
  
  real gla_eq_3(row_vector x, real inflection, vector slope, real y_cross, real A_y) {

  return(
    
       A_y + 
      (
        (A_y*y_cross) * 
        exp(   
          log1p_exp(inflection * slope[1]) -
          log1p_exp(   -( x * slope ) + inflection * slope[1]  ) 
          
        )
      ) 
    
  );
    
  }
  
  int[] get_elements_per_shard(int lenth_v, int shards){

	// Returned integer(max_size, last_element_size)
	int tentative_size = lenth_v / shards;
	int tentative_remaining = lenth_v - (tentative_size * shards);
	int elements_per_shard = tentative_remaining > 0 ? tentative_size + 1 : tentative_size;
	int remaining =  (elements_per_shard * shards) - lenth_v;

	int length_obj[shards];

	for(s in 1:shards) {
		length_obj[s] =
			s != shards ?
			elements_per_shard :
			elements_per_shard - remaining;  // Actual elements in it for last object
	}

 	return length_obj;

}

  vector[] get_mu_sigma_vector_MPI(vector mus, vector sigmas, int shards){

		int elements_per_shard[shards] = get_elements_per_shard(rows(mus), shards); // Length of the returned object
		int size_MPI_obj = elements_per_shard[1]; // the first element is always the full size of the object
		vector[size_MPI_obj * 2] v_MPI[shards] ; // Set values to -999 for the ones that are not filled

		int i = 0; // Index sweeping the vector

		for(s in 1:shards){

			// If last shard fill in
			if(s == shards) v_MPI[s] = rep_vector(-999.0, size_MPI_obj * 2);

			v_MPI[s, 1:elements_per_shard[s]] = mus[ (i + 1) : i + elements_per_shard[s] ];
			v_MPI[s, (elements_per_shard[s]+1):(elements_per_shard[s]+elements_per_shard[s])] = sigmas[ (i + 1) : i + elements_per_shard[s] ];

			i += elements_per_shard[s];
		}

		return v_MPI;
	}

  int get_real_buffer_size(vector v, real threshold){
  	// This function finds how may fake indexes -1 there are in a vector, added for map_rect needs
  
  	real i = threshold; // Value of the index
  	int n = 0; // Length of the buffer
  	int s = rows(v); // Size of the whole vector
  
  	while(i == threshold){
  		i = v[s-n];
  		if(i==threshold) n += 1;
  	}
  
  	return n;
  }

	vector lp_reduce_simple( vector global_parameters , vector mus_sigmas , real[] real_data , int[] int_data ) {

		real lp;
		real threshold = -999;
		int size_buffer = get_real_buffer_size(mus_sigmas, threshold);
		int size_vector = (rows(mus_sigmas)-size_buffer)/2;

		if(min(mus_sigmas[1:(size_vector*2)]) == threshold) print("ERROR! The MPI implmentation is buggy");

		// Reference / exposure rate
		lp = neg_binomial_2_lpmf(
			int_data[1:size_vector] |
			mus_sigmas[1:size_vector],
			mus_sigmas[size_vector+1:size_vector+size_vector]
		);

	 return [lp]';

	}
	
	int[] append_int(int[] i1, int[] i2){
  	int i3[size(i1)+size(i2)];

  	i3[1:size(i1)] = i1;
  	i3[size(i1)+1:size(i1)+size(i2)] = i2;

  	return i3;
  }

	int[] int_2D_to_1D(int[,] elems);
  int[] int_2D_to_1D(int[,] elems) {
    int num_elems = size(elems[1]);

    if (num_elems == 1) return(elems[,1]);

    if (num_elems == 2) return(append_int(elems[,1], elems[,2]));

    // else
 	  return(append_int(elems[,1],  int_2D_to_1D(elems[,2:num_elems]) ));
  }
  
	int[,] get_int_MPI(int[] v, int shards){
    // Simple MPI for int vector
  
  	int elements_per_shard[shards] = get_elements_per_shard(size(v), shards); // Length of the returned object
  	int size_MPI_obj = elements_per_shard[1]; // the first element is always the full size of the object
  	int v_MPI[shards,size_MPI_obj] = rep_array(-1, shards,size_MPI_obj); // Set values to -1 for the ones that are not filled
  
  	int i = 0; // Index sweeping the vector
  
  	for(s in 1:shards){
  		v_MPI[s, 1:elements_per_shard[s]] = v[ (i + 1) : i + elements_per_shard[s] ];
  		i += elements_per_shard[s];
  	}
  
  	return v_MPI;
  }

  	 vector horseshoe_get_tp(vector zb, vector[] local, real[] global, real scale_global, real c2) {
  	    int K = rows(zb);
  	    vector[K] lambda = local[1] .* sqrt(local[2]);
  	    vector[K] lambda2 = square(lambda);
  	    real tau = global[1] * sqrt(global[2]) * scale_global;
  	    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
  	    return zb .* lambda_tilde * tau;
  	  }
  
    real horseshoe_get_lp(vector zb, vector[] local, real df, real[] global, real df_global, real c2, real df_slab){
  
      // real<lower=0> hs_df; // == 1  // If divergencies increase this
  	  // real<lower=0> hs_df_global; // == 1
  	  // real<lower=0> hs_df_slab; // == 4 // df of the outliers
  
    	vector[6] lp;
  
    	lp[1] = normal_lpdf(zb | 0, 2);
  	  lp[2] = normal_lpdf(local[1] | 0, 1) - 101 * log(0.5);
  	  lp[3] = inv_gamma_lpdf(local[2] | 0.5 * df, 0.5 * df);
  	  lp[4] = normal_lpdf(global[1] | 0, 1)  - 1 * log(0.5);
  	  lp[5] = inv_gamma_lpdf(global[2] | 0.5 * df_global, 0.5 * df_global);
  	  lp[6] = inv_gamma_lpdf(c2 | 0.5 * df_slab, 0.5 * df_slab);
  
  	  return(sum(lp));
    }
//   real neg_binomial_2_MPI_lpmf(int[] y, vector mus, vector sigmas, int shards){
//     real real_data[shards,1];
//     	
//     return(
//       sum(map_rect(
// 	    	lp_reduce_simple,
// 		    [0]', // global parameters
//     		get_mu_sigma_vector_MPI(mus,	sigmas,	shards),
//     		real_data,
//     		get_int_MPI( y, shards)
//     	))
//     );
//   }

}
data {
	int <lower=0, upper = 1> prior_only; // For testing purpose (if ==0 then use data to generate quantities)

	int<lower = 0> G;                   // All Genes
	int<lower = 0> T;                   // All Samples (e.g. tube/individual)
	int<lower=0> R_1;                   // All Covariates (Intercept, Chosen Covariate 1 ... ) 
	int<lower = 0> y[T, G];             // RNA-seq counts
	matrix[T,R_1+1] X;                 // Design Matrix
	matrix<lower=0>[T,G] multiplier;           // Scale factor 
	
	int shards;

}
transformed data{

	/// Horseshoe tuning
  real<lower=0> hs_df = 1; // == 1  // If divergencies increase this
  real<lower=0, upper=1> par_ratio = 0.01; // real<lower=0> hs_scale_global; // from par ratio // 0.1 mean 1 outliers every 10 non outliers  !! KEY PARAMETER
  real<lower=0> hs_scale_slab = 0.5; // == 2 // regularisation/scale of outliers                !! KEY PARAMETER

  // MPI
	real real_data[shards,1];

}
parameters {
	// Linear model
	
	row_vector<lower=-3, upper=3>[G] inflection; //Value of the inflection point on the x axis
	
	vector<lower=0>[G] y_cross_z; 
	
	//Vector of slopes 
	matrix[R_1,G] beta; 
	
	//Vertical Translatioon
	vector[G] A;

  // Overdispersion
	vector[G] od;
	
  // Horseshoe
  vector<lower=0>[1] hs_local[2]; // local parameters for horseshoe prior
  real<lower=0> hs_global[2]; // horseshoe shrinkage parameters
  real<lower=0> hs_c2; // horseshoe shrinkage parameters

}
transformed parameters {
  
  // Horseshoe
  vector[G] y_cross= horseshoe_get_tp(y_cross_z, hs_local, hs_global, par_ratio / sqrt(T), hs_scale_slab^2 * hs_c2);


	matrix[T, G] log_y_hat;  //log of the mean of y
	matrix[T,G] phi; //log of the precision paramter - i.e dispersion in neg binomial is 1/exp(phi)
  //hence preventing cases of multiple solutions
  
	// Calculation of generalised logit - fitting in log space (i.e. log of the means follows gla eq)
	for(t in 1:T) for(g in 1:G) log_y_hat[t,g] = gla_eq(X[t,2:(R_1+1)], inflection[g], beta[,g], y_cross[g], A[g]);
	
	//Calculation of Overdispersion 
	//for(g in 1:G) phi[,g] = -0.3186 * log_y_hat[,g] + od[g];
	for(g in 1:G) phi[,g] = 0 * log_y_hat[,g] + od[g];



	
}
model {
real lp  = 0;
	// Linear system
	//Restricted priors on beta1_z[r], and inflection (were originally n(0,2)), preventing larger generated values
	//As these dramatically increase log_y_hat - which causes problems with neg_binomial_2_log / neg_binomial_2_log_rng
	for(r in 1:R_1) beta[r] ~ normal (0,4);
	
	inflection ~ normal(0,2);
	//y_cross ~ double_exponential(0,0.01); 
	
	//Vertical Translation. DO NOT SET BECAUSE IT CREATES UNDETERMINABILITY
	//A ~ normal(0,2);
	
	//overdispersion 
	od ~ normal(0,1);
	
	
  if(prior_only == 0)
	target += sum(map_rect(
		lp_reduce_simple,
		[0]', // global parameters
		get_mu_sigma_vector_MPI(
			exp(to_vector(log_y_hat)), //to_vector(multiplier).*exp(to_vector(log_y_hat)),
			1 ./ exp(to_vector(phi)),
			shards
		),
		real_data,
		get_int_MPI( int_2D_to_1D(y), shards)
	));
	
  // Horseshoe
  target += horseshoe_get_lp(y_cross_z, hs_local, hs_df, hs_global, 1, hs_c2, 4);
  
  
}
generated quantities{
  int y_gen[T,G];          // RNA-seq counts

	//Generate the data 
	for (t in 1:T) {
	  for(g in 1:G) {
	  y_gen[t,g] = neg_binomial_2_rng(multiplier[t,g]*exp(log_y_hat[t,g]), 1 ./ exp(phi[t,g])); 
	  } 
	  }
	  
	}
