functions{
  real rescale(real x, real lb, real ub){
    return((x-lb)/(ub-lb));
  }
  
  real rescale_inv(real x, real lb, real ub){
    return(lb+(ub-lb)*x);
  }
}

data{
  int n;
  int J;
  int H;
  int<lower=1> K;
  real lb;
  real ub;
  vector<lower=0,upper=K>[n] c; 
  vector<lower=0>[n] h;
  matrix[n,J] X;
  matrix[n,H] Z;
  vector[n] normaliz;
  vector[J] mu_beta;
  cov_matrix[J] Sigma_beta;
}

transformed data{
  vector<lower=0,upper=1>[n] cstar;
  
  for(i in 1:n){
    cstar[i] = fmin(1-1e-12,fmax(1e-12,rescale(c[i],lb,ub)));
  }
}

parameters{
  vector[J] betas;
  real<lower=1e-2> alpha_h;
  real<lower=1e-2> beta_h;
}

transformed parameters{
  vector<lower=0>[n] mu;
  mu =  exp(X*betas) .* normaliz;
}

model{
  
  betas ~ multi_normal(mu_beta,Sigma_beta);
  
  alpha_h ~ exponential(0.05);
  beta_h ~ exponential(0.1);
  
  for (i in 1:n){
    target += gamma_lpdf(h[i] | alpha_h, beta_h) + beta_lpdf(cstar[i] | 1e-6+(mu[i]/K)*h[i],1e-6+h[i]-h[i]*(mu[i]/K)); 
  }
}

generated quantities{
  
  array[n] real h_rep;
  array[n] real c_rep;
  
  for (i in 1:n){
    h_rep[i] = gamma_rng(alpha_h,beta_h);
    c_rep[i] = rescale_inv(beta_rng(1e-6+(1.0*mu[i]/K)*h_rep[i],1e-6+h_rep[i]-h_rep[i]*(1.0*mu[i]/K)),lb,ub);
  }
  
}




