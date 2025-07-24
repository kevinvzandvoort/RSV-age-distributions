#' setup distributions
#' specify Burr distributions
#' PDF
dburr = function(x, k, c, lambda=1, log=FALSE){
  ((c*k)/lambda) * (x/lambda)^(c-1) * (1+(x/lambda)^c)^(-k-1)
}
#' CDF
pburr = function(x, k, c, lambda){
  y = 1 - (1+(x/lambda)^c)^-k
  return(y)
}
#' Truncated PDF
dburr_truncated = function(x, k, c, lambda=1, log=FALSE, min=0, max=Inf){
  y = dburr(x, k, c, lambda)/(pburr(max, k, c, lambda) - pburr(min, k, c, lambda))
  y[x<min] = 0
  y[x>max] = 0
  if(log) y=log(y)
  return(y)
}
#' Truncated CDF
pburr_truncated = function(x, k, c, lambda=1, min=0, max=Inf){
  y = (pburr(x, k, c, lambda) - pburr(min, k, c, lambda))/
    (pburr(max, k, c, lambda) - pburr(min, k, c, lambda))
  y[x<min] = 0
  y[x>max] = 1
  return(y)
}

dgamma_truncated = function(x, alpha, beta, min=0, max=Inf){
  y = dgamma(x, shape = alpha, rate = beta)/(pgamma(max, shape = alpha, rate = beta) - pgamma(min, shape = alpha, rate = beta))
  y[x<min] = 0
  y[x>max] = 0
  return(y)
}
pgamma_truncated = function(x, alpha, beta, min=0, max=Inf){
  y = (pgamma(x, shape = alpha, rate = beta) - pgamma(min, shape = alpha, rate = beta))/
    (pgamma(max, shape = alpha, rate = beta) - pgamma(min, shape = alpha, rate = beta))
  y[x<min] = 0
  y[x>max] = 1
  return(y)
}

dllogis_truncated = function(x, alpha, beta, min=0, max=Inf){
  y = dllogis(x, scale = alpha, shape = beta)/(pllogis(max, scale = alpha, shape = beta) - pllogis(min, scale = alpha, shape = beta))
  y[x<min] = 0
  y[x>max] = 0
  return(y)
}
pllogis_truncated = function(x, alpha, beta, min=0, max=Inf){
  y = (pllogis(x, scale = alpha, shape = beta) - pllogis(min, scale = alpha, shape = beta))/
    (pllogis(max, scale = alpha, shape = beta) - pllogis(min, scale = alpha, shape = beta))
  y[x<min] = 0
  y[x>max] = 1
  return(y)
}

dlnorm_truncated = function(x, mu, sigma, min=0, max=Inf){
  y = dlnorm(x, meanlog = mu, sdlog = sigma)/(plnorm(max, meanlog = mu, sdlog = sigma) - plnorm(min, meanlog = mu, sdlog = sigma))
  y[x<min] = 0
  y[x>max] = 0
  return(y)
}
plnorm_truncated = function(x, mu, sigma, min=0, max=Inf){
  y = (plnorm(x, meanlog = mu, sdlog = sigma) - plnorm(min, meanlog = mu, sdlog = sigma))/
    (plnorm(max, meanlog = mu, sdlog = sigma) - plnorm(min, meanlog = mu, sdlog = sigma))
  y[x<min] = 0
  y[x>max] = 1
  return(y)
}