//several custom function
functions {
  real cburr(real x, real k, real c, real lambda) {
    return ((c*k)/lambda) * (x/lambda)^(c-1) * (1+(x/lambda)^c)^(-k-1);
  }
    
  real pburr(real x, real k, real c, real lambda) {
    return 1 - (1+(x/lambda)^c)^-k;
  }
    
  real cburr_truncated(real x, real k, real c, real lambda, real min, real max){
    real y = cburr(x, k, c, lambda)/(pburr(max, k, c, lambda) - pburr(min, k, c, lambda));
    if(x < min) y = 0.0;
    if(x > max) y = 0.0;
    return y;
  }
    
  real pburr_truncated(real x, real k, real c, real lambda, real min, real max){
    real y = (pburr(x, k, c, lambda) - pburr(min, k, c, lambda))/(pburr(max, k, c, lambda) - pburr(min, k, c, lambda));
    if(x < min) return 0.0;
    if(x > max) return 1.0;
    return y;
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int Ndat; //total number of observations per dataset
  real max_age; //maximum age per dataset
  real min_age; //maximum age per dataset
  int obs[Ndat]; //observed number of cases
  real age_low[Ndat]; //lower age_band of age_group
  real age_high[Ndat]; //upper age_band of age_group
}

//doesnt allow me to use this as transformed parameters due to the rng
transformed data {
  real cum_obs;
  real o = 0;
  for(a in 1:Ndat){
    o += obs[a];
  }
  cum_obs = o;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu_k;
  real mu_c;
  real mu_lambda;
  
  //real<lower=0, upper=50> k;
  //real<lower=1, upper=50> c;
  //real<lower=1, upper=50> lambda;
}

transformed parameters {
  real k;
  real c;
  real lambda;
  real expected_cases[Ndat];
  
  k = exp(mu_k);
  c = 1 + exp(mu_c);
  lambda = 1 + exp(mu_lambda);
  
  for(a in 1:Ndat){
    real p = pburr_truncated(age_high[a], k, c, lambda, min_age, max_age) - pburr_truncated(age_low[a], k, c, lambda, min_age, max_age);
    //calculate expected number of cases in this age group
    expected_cases[a] = cum_obs * p;
    //add a small value to prevent 0
    expected_cases[a] += 1e-20;  // Add a small value to avoid zero
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  mu_k ~ normal(0, 2.5);
  mu_c ~ normal(0, 2.5);
  mu_lambda ~ normal(0, 2.5);
  
  real y = pburr(max_age, k, c, lambda);
  if(y == 0.0) reject("CDF 0 at maxage; k: ", k, "; c: ", c, "; lambda: ", lambda);
    
  //loop through every age group in this dataset
  for(a in 1:Ndat){
    obs[a] ~ poisson(expected_cases[a]);
  }
}

generated quantities {
  vector[Ndat] log_lik;
  for (a in 1:Ndat) {
    log_lik[a] = poisson_lpmf(obs[a] | expected_cases[a]);
  }
}
