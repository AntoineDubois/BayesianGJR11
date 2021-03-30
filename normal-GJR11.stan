data {
  //The main data
  int<lower=0> n;
  int<lower=0> p;
  matrix[n,p] X;
  vector[n] y; 
  //The hyperparameters of the priors
  vector[p] mu_gamma;
  cov_matrix[p] sigma_gamma;
  vector<lower=0>[3] mu_alpha;
  cov_matrix[3] sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> mu_beta;
}
parameters {
  vector<lower=0>[3] alpha;
  real<lower=0> beta;
  vector[p] gamma;
}

transformed parameters {
  real<lower=0> h[n];
  real<lower=0> sigma[n];
  vector[n] linpred;

  linpred = X*gamma; 
  //Constuction of the standard error of u_t for any t
  h[1] = alpha[1];
  sigma[1] = sqrt(h[1]);
  for (t in 2:n){
    if(y[t-1] - linpred[t-1] >= 0)
    {h[t] = alpha[1] + alpha[2]*pow(y[t-1] -linpred[t-1], 2)+ beta*h[t-1];}
    else
    {h[t] = alpha[1] + alpha[3]*pow(y[t-1] -linpred[t-1], 2)+ beta*h[t-1];}
    sigma[t] = sqrt(h[t]);}
}

model {
  y ~ normal(linpred, sigma);//normal inovation
  alpha ~ multi_normal(mu_alpha,sigma_alpha);
  beta ~ normal(mu_beta,sigma_beta);
  gamma ~ multi_normal(mu_gamma,sigma_gamma);
}
