//
// Stan program for RSV pooled meta-regression
// Fit Burr age-distributions to case-counts by weeks of age
//

//several custom function
functions {
  real dburr(real x, real k, real c, real lambda) {
    return ( (c*k)/lambda ) * ( x/lambda )^(c - 1) * (1 + ( x/lambda )^c)^(-k - 1);
  }
    
  real pburr(real x, real k, real c, real lambda) {
    return 1 - (1 + (x/lambda)^c)^(-k);
  }
    
  real dburr_truncated(real x, real k, real c, real lambda, real min, real max){
    real y = dburr(x, k, c, lambda)/(pburr(max, k, c, lambda) - pburr(min, k, c, lambda));
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
  int<lower=0> D; //total number of datasets
  int<lower=0> Ntotal; //total number of observations in all datasets
  int Ndat[D]; //total number of observations per dataset
  real max_ages[D]; //maximum age per dataset
  real min_ages[D]; //maximum age per dataset
  int obs[Ntotal]; //observed number of cases
  real age_low[Ntotal]; //lower age_band of age_group
  real age_high[Ntotal]; //upper age_band of age_group
}

transformed data {
  real cum_obs[D];
  int i = 1;
  
  for(d in 1:D){
    real o = 0;
    
    for(a in 1:Ndat[d]){
      o += obs[i];
      i += 1;
    }
    
    cum_obs[d] = o;
  }
}

// The parameters accepted by the model.
// mu, c, and lambda
parameters {
  real<lower=0> mu_k;
  real<lower=1> mu_c;
  real<lower=1> mu_lambda;
  real<lower=0> sigma_k;
  real<lower=0> sigma_c;
  real<lower=0> sigma_lambda;
  vector<lower=-mu_k/sigma_k>[D] k_raw; //k should be >0
  vector<lower=(1-mu_c)/sigma_c>[D] c_raw; //c should be >1
  vector<lower=(1-mu_lambda)/sigma_lambda>[D] lambda_raw; //lambda should be >1
}

transformed parameters {
  real expected_cases[Ntotal];
  vector[D] k;
  vector[D] c;
  vector[D] lambda;
  
  {
    int n = 1;
    //loop through every dataset
    for(d in 1:D){
      k[d] = mu_k + sigma_k * k_raw[d];
      c[d] = mu_c + sigma_c * c_raw[d];
      lambda[d] = mu_lambda + sigma_lambda * lambda_raw[d];
      
      //loop through every age group in this dataset
      for(a in 1:Ndat[d]){
        //get proportion of cases between age_high and age_low
        real p = pburr_truncated(age_high[n], k[d], c[d], lambda[d], min_ages[d], max_ages[d]) - pburr_truncated(age_low[n], k[d], c[d], lambda[d], min_ages[d], max_ages[d]);
      
        //calculate expected number of cases in this age group
        expected_cases[n] = cum_obs[d] * p;
      
        //add a small value to avoid zero values for poisson probability
        expected_cases[n] += 1e-2;
      
        n += 1;
      }
    } 
  }
}

// The model to be estimated.
model {
  int n = 1;
  
  //sd for global k, c, and lambda
  sigma_k ~ cauchy(0, 1) T[0, ];
  sigma_c ~ normal(0, 1) T[0, ];
  sigma_lambda ~ normal(0, 1) T[0, ];
  
  //mean for global k, c, and lambda
  mu_k ~ normal(0, 5) T[0, ];
  mu_c ~ normal(1, 5) T[1, ];
  mu_lambda ~ normal(1, 5) T[1, ];
  
  //loop through every dataset
  for(d in 1:D){
    
    // study specific k, d, and lambda
    k_raw[d] ~ normal(0, 1);
    c_raw[d] ~ normal(0, 1);
    lambda_raw[d] ~ normal(0, 1);
    
    //check that CDF not 0 at max age
    real y = pburr(max_ages[d], k[d], c[d], lambda[d]);
    if(y == 0.0) reject("CDF 0 at maxage");
    
    //loop through every age group in this dataset
    for(a in 1:Ndat[d]){
      //assume poisson distribution of (rounded) observed cases
      obs[n] ~ poisson(expected_cases[n]);
      n += 1;
    }
  }
}

//log-likelihood for every observation
generated quantities {
  vector[Ntotal] log_lik;
  
  for (a in 1:Ntotal) {
    log_lik[a] = poisson_lpmf(obs[a] | expected_cases[a]);
  }
}

