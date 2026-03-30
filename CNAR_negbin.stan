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
  vector[H] mu_gamma;
  cov_matrix[J] Sigma_beta;
  cov_matrix[H] Sigma_gamma;
}

transformed data{
  vector<lower=0,upper=1>[n] cstar;
  
  for(i in 1:n){
    cstar[i] = fmin(1-1e-12,fmax(1e-12,rescale(c[i],lb,ub)));
  }
}

parameters{
  vector[J] betas;
  vector[H] gammas;
  real<lower=1e-2> alpha_h;
  real<lower=1e-2> beta_h;
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
  
  alpha_h ~ exponential(0.05);
  beta_h ~ exponential(0.1);
  
  // Marginalizing over y -- see: https://mc-stan.org/docs/stan-users-guide/latent-discrete.html
  //                              https://stats.stackexchange.com/questions/563129/marginalizing-out-discrete-response-variables-in-stan (Option 2)
  for (i in 1:n){
    //vector[K+1] lp;
    //real negbin2_Nconstant = neg_binomial_2_lcdf(K|mu[i],phi[i]);
    
    real joint_loglikel_unno = negative_infinity();
    real negbin_nConstant = negative_infinity();
    
    for (y in 0:K){ //marginalizing over y to get the pointwise joint loglikel for (c,h)_i
      real ystar = (y*1.0)/K;
      real negbin = neg_binomial_2_lpmf(y | mu[i], phi[i]); //p(y|mu,phi) --latent model for Y
      real be = beta_lpdf(cstar[i] | 1e-6+ystar*h[i],1e-6+h[i]-h[i]*ystar); //p(m|h,y) --observed model for C
      joint_loglikel_unno = log_sum_exp(joint_loglikel_unno,negbin+be); //partial sums for the unnormalized joint loglikel
      negbin_nConstant = log_sum_exp(negbin_nConstant,negbin); //partial sums for the normalization constant F(K|mu,phi) of the negbin model (truncated at K)
      //lp[y+1] = neg_binomial_2_lpmf(y | mu[i], phi[i]) + beta_lpdf(cstar[i] | 1e-6+ystar*h[i],1e-6+h[i]-h[i]*ystar);
    }
    target += gamma_lpdf(h[i] | alpha_h, beta_h) + joint_loglikel_unno - negbin_nConstant; //joint loglikel with normalization constant
    //target += gamma_lpdf(h[i] | alpha_h, beta_h) + log_sum_exp(lp); //joint likelihood pointwise after marginalizing over y (log_sum_exp for numerical stability)
  }
}

generated quantities{
  
  array[n] real log_lik_joint;
  array[n] int y_pred;
  array[n] int y_rep;
  array[n] real h_rep;
  array[n] real c_rep;
  
  for (i in 1:n){
    real joint_loglikel_unno = negative_infinity();
    real negbin_nConstant = negative_infinity();
    vector[K+1] lp;
    vector[K+1] negbin;
    
    for (y in 0:K){ 
      real ystar = (y*1.0)/K;
      negbin[y+1] = neg_binomial_2_lpmf(y | mu[i], phi[i]); //p(y|mu,phi) --latent model for Y
      real be = beta_lpdf(cstar[i] | 1e-6+ystar*h[i],1e-6+h[i]-h[i]*ystar); //p(m|h,y) --observed model for C
      lp[y+1] = negbin[y+1]+be;
      joint_loglikel_unno = log_sum_exp(joint_loglikel_unno,lp[y+1]); //partial sums for the unnormalized joint loglikel
      negbin_nConstant = log_sum_exp(negbin_nConstant,negbin[y+1]); //partial sums for the normalization constant F(K|mu,phi) of the negbin model (truncated at K)
    }
    log_lik_joint[i] = gamma_lpdf(h[i] | alpha_h, beta_h) + joint_loglikel_unno - negbin_nConstant; //joint loglikel with normalization constant (log_sum_exp for numerical stability)
    y_pred[i] = categorical_rng(softmax(lp)) - 1; //sample y from the discrete posterior
    
    y_rep[i] = categorical_rng(softmax(negbin)) - 1; //sample y from a truncated negbin
    h_rep[i] = gamma_rng(alpha_h,beta_h); //sample h
    c_rep[i] = rescale_inv(beta_rng(1e-6+(1.0*y_rep[i]/K)*h_rep[i],1e-6+h_rep[i]-h_rep[i]*(1.0*y_rep[i]/K)),lb,ub); //sample m|y,h
  }
  
}




