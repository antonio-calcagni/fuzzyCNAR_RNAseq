defuzzify_centroid = function(mi,si,what=c("mean","mode")){
  a1 = 1+si*mi; b1 = 1+si*(1-mi) 
  if(what=="mean"){(a1)/(a1+b1)}
  else if(what=="mode"){(a1-1)/(a1+b1-2)}
}

curr_dir <- "rnaseq_fuzzy/M1_short_simulation/"
load("M1_short_simulation/data_common.rds")
load("M1_short_simulation/M1_Theta_post.rds")

M <- 500; K <- 1185 
set.seed(260326)
mmd <- sample(1:NROW(Theta_post),size = M,replace = FALSE)
for(m in 1:M){
  theta_m <- Theta_post[mmd[m],]
  y <- h <- c <- cd <- rep(NA,n)
  for(i in 1:n){
    y_i <- K+1; while(y_i>K){y_i <- MASS::rnegbin(n = 1,mu = exp(X[i,]%*%theta_m[1:3])*u[i],theta = exp(theta_m[4]))}; y[i] <- y_i
    h_i <- rgamma(n = 1,theta_m[5],theta_m[6]); h[i] <- h_i
    c_i <- rbeta(1,(y_i/K)*h_i,h_i-h_i*(y_i/K))*K; c[i] <- c_i
    cd_i <- defuzzify_centroid(c_i,si = h_i,what = "mean"); cd[i] <- cd_i
  }
  save(y,c,h,cd,file = paste0(curr_dir,"/simdata/simdata_",m,".rds"))
}

