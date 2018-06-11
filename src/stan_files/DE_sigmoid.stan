functions{

  vector log_gen_inv_logit(row_vector y, row_vector inversion, row_vector intercept) {
    return to_vector( intercept + log1p_exp(-inversion) - log1p_exp(- y  ) );
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
row_vector[G] inversion_z;
row_vector[G] intercept;
row_vector<lower=0>[G] sigma_trick;

// Horseshoe
vector [ G] beta1_z;
real < lower =0 > aux1_global ;
real < lower =0 > aux2_global ;
vector < lower =0 >[ G] aux1_local ;
vector < lower =0 >[ G] aux2_local ;
real < lower =0 > caux ;
}
transformed parameters {

// Horseshoe variables
real < lower =0 > tau ; // global shrinkage parameter
vector < lower =0 >[ G] lambda ; // local shrinkage parameter
vector < lower =0 >[ G] lambda_tilde ; // ’ truncated ’ local shrinkage parameter
real < lower =0 > c; // slab scale

matrix[R, G] beta;
matrix[R-1, G] beta1;
row_vector[G] inversion;
row_vector[G] beta1_trick;

// Horseshoe calculation
lambda = aux1_local .* sqrt ( aux2_local );
tau = aux1_global * sqrt ( aux2_global ) * scale_global * 1 ;
c = slab_scale * sqrt ( caux );
lambda_tilde = sqrt ( c ^2 * square ( lambda ) ./ (c ^2 + tau ^2* square ( lambda )) );
beta1[1] =  to_row_vector(beta1_z .* lambda_tilde * tau)   ;

// trick
beta1_trick = beta1[1] .* sigma_trick * 0.5;
inversion = inversion_z .* sigma_trick * 0.5;

// make beta
beta = append_row(inversion, beta1_trick);
}
model {

// Horseshoe
beta1_z ~ normal (0 , 1);
aux1_local ~ normal (0 , 1);
aux2_local ~ inv_gamma (0.5* nu_local , 0.5* nu_local );
aux1_global ~ normal (0 , 1);
aux2_global ~ inv_gamma (0.5* nu_global , 0.5* nu_global );
caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df );

// Trick
sigma_trick ~ normal(0,1);
for(r in 2:R) beta1[r-1] ~ normal(0,10);

// Linear system
inversion_z ~ normal(0 , 1);
intercept ~ normal(0,5);
sum(intercept) ~ normal(0, 0.01 * G);

for (t in 1:T)
y[t] ~ multinomial( softmax( log_gen_inv_logit(X[t] * beta, inversion, intercept) ) );
}

generated quantities{
  int y_gen[T,G];          // RNA-seq counts

	for (t in 1:T)
			y_gen[t] = multinomial_rng( softmax( log_gen_inv_logit(X[t] * beta, inversion, intercept) ), exposure[t] );
}
