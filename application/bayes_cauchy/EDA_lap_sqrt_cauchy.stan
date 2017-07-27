// one smooth with random effects
// Laplace prior on finite differences in B-spline coefficients

data {
  int<lower=0> n; // total number of observations
  int<lower=0> q1; // number of columns of Zcheck1
  int<lower=0> q2; // number of columns of Ztilde2
  vector[n] y; // outcome
  int<lower=0> p1; // dimensions of B-spline basis
  int<lower=0> p2; // dimensions of B-spline basis
  int<lower=0> k1; // order of finite difference matrix D
  int<lower=0> k2; // order of finite difference matrix D
  matrix[n, p1 - 1] X1; // B-spline bases evaluated at x
  matrix[n, p2] X2; // B-spline bases evaluated at x
  matrix[p1-k1-1, p1 - 1] D1; // finite difference matrix
  matrix[p2-k2-1, p2] D2; // finite difference matrix
  matrix[n, q1] Zcheck1; // random effects design matrix -- penalized components
  matrix[n, q2] Ztilde2; // random effects design matrix -- unpenalized components
}

parameters {
  real beta0; // population level intercept
  vector[p1 - 1] beta1; // population level parameters
  vector[p2] beta2; // population level parameters
  vector[q1] btilde1; // random effects penalized
  vector[q2] btilde2; // random effects unpenalized

  // diffuse priors on sd's, not log(sd)
  real<lower=0> sigmaEpsilon; // sd term for overall errors
  real<lower=0> sigmaB; // 1/sd term for random intercepts "penalized"
  real<lower=0> sigmaB2; // 1/sd term for random intercepts "unpenalized"
  //real logSigmaLambda; // sd term for double exponential
  real<lower=0> sigmaLambda1;
  real logsigmaLambda2;
}

transformed parameters {
  vector[p1-k1-1] Dbeta1;
  vector[p2-k2-1] Dbeta2;
  real sigmaLambda2;
    
  Dbeta1 = D1 * beta1;
  Dbeta2 = D2 * beta2;
  sigmaLambda2 = exp(logsigmaLambda2);
}

model {
  btilde1 ~ normal(0, sigmaB);
  btilde2 ~ cauchy(0, sigmaB2);
  Dbeta1 ~ double_exponential(0, sigmaLambda1);
  Dbeta2 ~ double_exponential(0, sigmaLambda2);
  y ~ normal(beta0 + X1*beta1 + X2*beta2 + Zcheck1*btilde1 + Ztilde2*btilde2, sigmaEpsilon);
}
