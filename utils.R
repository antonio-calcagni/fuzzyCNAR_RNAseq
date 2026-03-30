beta_type_gcount_approximation <- function(xsup=NULL,mux=NULL,alpha0=0.1,eps=0.05,h_max=100,dd=0.5){
  
  SK <- xsup; y <- mux
  alpha_grid <- seq(alpha0,0.99,length=101) #trimming at alpha0
  n_SK <- max(SK)
  
  L_Ya <- sapply(alpha_grid,function(a)A_L(alpha_cut(SK,y,a)))
  U_Ya <- sapply(alpha_grid,function(a)A_U(alpha_cut(SK,y,a)))
  W_Y <- W_AB(alpha_grid,SK,y)
  
  C_eps <- SK[y>(1-eps)]; K <- length(C_eps)
  parx <- matrix(nrow = K,ncol = 3)
  
  for(k in 1:K){
    m_k <- C_eps[k]
    
    h_low  <- tryCatch({uniroot(function(h) W_AB(x_alpha = alpha_grid,x = SK,mu_x = beta_type_FS(SK,m_k,h,n_SK)) - (W_Y+(W_Y*dd)), c(1e-9, h_max))$root},error=function(e){return(1e-9)})
    h_high  <- tryCatch({uniroot(function(h) W_AB(x_alpha = alpha_grid,x = SK,mu_x = beta_type_FS(SK,m_k,h,n_SK)) - (W_Y-(W_Y*dd)), c(1e-9, h_max))$root},error=function(e){return(h_max)})
    h_k <- optimize(function(x)D_AB_2(alpha_grid,SK,y,beta_type_FS(SK,m_k,x,n_SK)), interval=c(h_low, h_high))$minimum
    
    parx[k,] <- c(m_k,h_k,D_AB_2(alpha_grid,SK,y,beta_type_FS(SK,m_k,h_k,n_SK)))
  }
  parx_est <- parx[which.min(parx[,3]),1:2]
  
  return(list(c_hat=parx_est[1],h_hat=parx_est[2]))
}


beta_type_FS <- function(x,m,h,k){ #discrete beta-type FS
  px <- function(x,m,h,k,epsilon=0.5){return((x+epsilon)^((m+epsilon)/(h*k+2*h*epsilon)) * (k+epsilon-x)^((k+epsilon-m)/(h*k+2*h*epsilon)))}
  return(px(x,m,h,k) / px(m,m,h,k))
}

alpha_cut <- function(x=NULL,mu_x=NULL,alpha=0.01){return(x[mu_x>=alpha])}
A_L <- function(x){return(min(x))}; A_U <- function(x){return(max(x))}

D_AB_2 <- function(x_alpha=NULL,x,muA,muB){
  return(sum(sapply(x_alpha,function(a)a*((A_L(alpha_cut(x,muA,a))-A_L(alpha_cut(x,muB,a)))^2 + (A_U(alpha_cut(x,muA,a))-A_U(alpha_cut(x,muB,a)))^2)))^0.5)
}

W_AB <- function(x_alpha=NULL,x=NULL,mu_x=NULL){
  return(sum(sapply(x_alpha,function(a)A_U(alpha_cut(x,mu_x,a))-A_L(alpha_cut(x,mu_x,a)))))
}

I_AB <- function(mu_A,mu_B){
  ## Soruce: Fan, J., & Xie, W. (1999). Distance measure and induced fuzzy entropy. Fuzzy sets and systems, 104(2), 305-314.
  return(sum(1 - (1-mu_A)*exp(mu_A-mu_B) - mu_A*exp(mu_B-mu_A)))
}

De_AB <- function(mu_A,mu_B){
  ## Fan, J., & Xie, W. (1999). Distance measure and induced fuzzy entropy. Fuzzy sets and systems, 104(2), 305-314.
  return(I_AB(mu_A,mu_B)+I_AB(mu_B,mu_A))
}

bertoluzza_dist <- function(xsup,muy_a,muy_b,eps=1/2,num_alpha=101){
  alphas <- seq(1e-2,1-1e-2,length=num_alpha)
  fx <- sapply(alphas,function(a){
    xa <- range(alpha_cut(xsup,mu_x = muy_a,alpha = a))
    xb <- range(alpha_cut(xsup,mu_x = muy_b,alpha = a))
    eps*(sum(xa)/2 - sum(xb)/2)^2 + (1-eps)*(diff(xa) - diff(xb))^2
  })
  return(sqrt(pracma::trapz(x = alphas,y = fx)))
}

defuzzify_centroid = function(mi,si,what=c("mean","mode")){
  a1 = 1+si*mi; b1 = 1+si*(1-mi) 
  if(what=="mean"){(a1)/(a1+b1)}
  else if(what=="mode"){(a1-1)/(a1+b1-2)}
}
