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
int<lower = 0> G;                   // all genes
int<lower = 0> T;                   // tube
int<lower=0> R_1;
int<lower = 0> y[T, G];             // RNA-seq counts
matrix[T,R_1+1] X;
int exposure[T];                 // How many reads have been sequenced for each sample

// Horseshoe
vector < lower =0 >[R_1] par_ratio ; // proportion of 0s
vector < lower =1 >[R_1] nu_global ; // degrees of freedom for the half -t prior
vector < lower =1 >[R_1] nu_local ; // degrees of freedom for the half - t priors
vector < lower =0 >[R_1] slab_scale ; // slab scale for the regularized horseshoe
vector < lower =0 >[R_1] slab_df; // slab degrees of freedom for the regularized
}

transformed data{
vector < lower =0 >[R_1] scale_global = par_ratio / sqrt(1.0 * T); // scale for the half -t prior for tau
}

parameters {
// Linear model
vector[G] inversion_z;
vector[G] intercept;
vector<lower=0>[G] sigma_trick;

// Horseshoe
vector [ G] beta1_z[R_1];
real < lower =0 > aux1_global[R_1] ;
real < lower =0 > aux2_global[R_1] ;
vector < lower =0 >[ G] aux1_local[R_1] ;
vector < lower =0 >[ G] aux2_local[R_1] ;
real < lower =0 > caux[R_1] ;
}

transformed parameters {

vector[G] beta1[R_1];
vector[G] beta1_trick[R_1];
vector[G] inversion;
matrix[R_1+1, G] beta;
matrix[T, G] y_hat;

// Horseshoe calculation
for(r in 1:R_1)
	beta1[r] =
		reg_horseshoe(
			beta1_z[r],
			aux1_global[r] ,
			aux2_global[r],
			aux1_local[r] ,
			aux2_local[r] ,
			caux[r],
			scale_global[r],
			slab_scale[r]
		);

	// trick
	for(r in 1:R_1) beta1_trick[r] = beta1[r] .* sigma_trick * 0.5;
	inversion = inversion_z .* sigma_trick * 0.5;

	// make beta
	beta[1] = to_row_vector(inversion);
	for(r in 1:R_1) beta[r+1] = to_row_vector(beta1_trick[r]);

	// Matrix multiplication for speed up
	y_hat = X * beta;

}
model {

// Horseshoe
for(r in 1:R_1) beta1_z[r] ~ normal (0 , 1);
for(r in 1:R_1) aux1_local[r] ~ normal (0 , 1);
for(r in 1:R_1) aux2_local[r] ~ inv_gamma (0.5* nu_local , 0.5* nu_local );
aux1_global ~ normal (0 , 1);
aux2_global ~ inv_gamma (0.5* nu_global , 0.5* nu_global );
caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df );

// Trick
sigma_trick ~ normal(0,1);
for(r in 1:R_1) beta1[r] ~ normal(0,10);

// Linear system
inversion_z ~ normal(0 ,1);
intercept ~ normal(0,5);
sum(intercept) ~ normal(0, 0.01 * G);

// Likelihood
for (t in 1:T) y[t] ~ multinomial( softmax( log_gen_inv_logit(y_hat[t], inversion, intercept) ) );

}

generated quantities{
  int y_gen[T,G];          // RNA-seq counts

	for (t in 1:T)
			y_gen[t] = multinomial_rng( softmax( log_gen_inv_logit(y_hat[t], inversion, intercept) ), exposure[t] );
}
