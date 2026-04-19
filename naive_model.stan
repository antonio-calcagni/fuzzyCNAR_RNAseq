
data{
  int n;
  int J;
  int H;
  int K;
  array[n] int y; 
  matrix[n,J] X;
  matrix[n,H] Z;
  vector[n] normaliz;
  vector[J] mu_beta;
  vector[H] mu_gamma;
  cov_matrix[J] Sigma_beta;
  cov_matrix[H] Sigma_gamma;
}

parameters{
  vector[J] betas;
  vector[H] gammas;
}

transformed parameters{
  vector<lower=0>[n] mu;
  vector[n] phi;
  mu =  exp(X*betas) .* normaliz;
  phi = exp(Z*gammas);
}

model{
  betas ~ multi_normal(mu_beta,Sigma_beta);
  gammas ~ multi_normal(mu_gamma,Sigma_gamma);
  
  for (i in 1:n){
    y[i] ~ neg_binomial_2(mu[i], phi[i])T[,K];
  }
}

generated quantities{
  array[n] int y_rep;
  
  for (i in 1:n){
    y_rep[i] = K + 1;
    while (y_rep[i] > K){
      y_rep[i] = neg_binomial_2_rng(mu[i], phi[i]);
    }
  }
  
}





