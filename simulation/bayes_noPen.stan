// one smooth with random effects
// no prior on finite differences in B-spline coefficients

data {
  int<lower=0> n; // total number of observations
  int<lower=0> N; // number of subjects
  vector[n] y; // outcome
  int<lower=0> p; // dimensions of B-spline basis
  int<lower=0> k; // order of finite difference matrix D
  matrix[n, p-1] F; // B-spline bases evaluated at x
  matrix[p-k-1, p-1] D; // finite difference matrix
  int<lower=1, upper=N> id[n]; // vector containing id numbers (1-N) for each observation
  matrix[n, N] Z;
}
parameters {
  real beta0; // population level intercept
  vector[p-1] beta; // population level parameters
  vector[N] b; // random intercepts

  // diffuse priors on sd's, not log(sd)
  real<lower=0> sigmaEpsilon; // sd term for overall errors
  real<lower=0> sigmaB; // sd term for random intercepts
  //real logSigmaLambda; // sd term for double exponential
  real<lower=0> sigmaLambda;
}

transformed parameters {
  vector[p-k-1] Dbeta;
  //real sigmaLambda;
    
  //sigmaLambda <- exp(logSigmaLambda);
  Dbeta <- D * beta;
}


model {
  b ~ normal(0, sigmaB);
  y ~ normal(beta0 + F * beta + Z * b, sigmaEpsilon);
}
