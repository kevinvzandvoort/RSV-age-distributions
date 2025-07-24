//several custom function
functions {
  real cdf_truncated(real x, real alpha, real beta, real min, real max){
    real y = (loglogistic_cdf(x, alpha, beta) - loglogistic_cdf(min, alpha, beta))/(loglogistic_cdf(max, alpha, beta) - loglogistic_cdf(min, alpha, beta));
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
  real cum_obs = 0;
  for(a in 1:Ndat){
    cum_obs += obs[a];
  }
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu_alpha;
  real mu_beta;
}

transformed parameters {
  real alpha;
  real beta;
  real expected_cases[Ndat];
  
  alpha = exp(mu_alpha);
  beta = exp(mu_beta);
  
  //loop through every age group in this dataset
  for(a in 1:Ndat){
    real p = cdf_truncated(age_high[a], alpha, beta, min_age, max_age) - cdf_truncated(age_low[a], alpha, beta, min_age, max_age);
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
  mu_alpha ~ normal(0, 2.5);
  mu_beta ~ normal(0, 2.5);
  
  real y = loglogistic_cdf(max_age, alpha, beta);
  if(y == 0.0) reject("CDF 0 at maxage; alpha: ", alpha, "; beta: ", beta);
    
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
