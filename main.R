

# Set environment ---------------------------------------------------------
rm(list=ls())
source("utils.R")

require(cmdstanr)
set_cmdstan_path("/usr/local/cmdstan/"); cmdstan_version()

library(posterior); library(xtable)


# Beta-type approx of raw granular counts ---------------------------------
fls_list <- list.files(path = "data/raw_gcounts",full.names = TRUE,pattern = "SRR")
n <- length(fls_list)

set.seed(121)
xmin <- 0; xmax <- 1185 #based on the observed data
xsup <- seq(xmin,xmax)
sbjs <- rep(NA,n); Parx <- matrix(NA,n,2); colnames(Parx) <- c("c","h")
for(i in 1:n){
  cat(i,"\t")
  sbjs[i] <- paste0("SRR",strsplit(x = strsplit(x = fls_list[i],split = "/SRR")[[1]][2],split = "_gran")[[1]][[1]])
  mux <- as.vector(read.csv(file = fls_list[i],header = TRUE))[[1]]
  mux <- mux[xsup+1]; mux[is.na(mux)] <- 0.0
  c_hat <- mean(xsup[mux>=0.95]) 
  out <- tryCatch({optim(par = -2,fn = function(x,mu_A){De_AB(mu_A,beta_type_FS(xsup,c_hat,exp(x),max(xsup)))},method = "L-BFGS-B",lower = -Inf,mu_A=mux)},
                  error=function(e){NULL}); if(!is.null(out)){out$par <- exp(out$par)}else{out$par=exp(-4.5)} #if the fuzzy set is too narrow, h_hat becames fixed at its lower computable value
  Parx[i,] <- c(c_hat,out$par+runif(1,min = 1e-3,max = 2e-3)) #add a very small epsilon to avoid ties
}
dataout <- data.frame(sbjs,Parx,1/Parx[,2]); names(dataout)[4] <- "1/h"
write.csv2(dataout,file = "data/beta_approx_out.csv")



# CNAR models (Stan) ------------------------------------------------------

## Data and compiled files
mod <- cmdstanr::cmdstan_model(stan_file = "CNAR_negbin.stan")
dataout <- read.csv(file = "data/HAS3_complete.csv",header = TRUE,dec = ",")
str(dataout)

iid <- !is.na(dataout$hba1c_categ)
dataout <- dataout[iid,]
dataout$hba1c_categ <- factor(dataout$hba1c_categ,levels = c("normal","pre_diabetes","diabetes"))
dataout$hba1c_categ <- relevel(dataout$hba1c_categ,ref = "normal")
dataout$bio_sex <- as.factor(dataout$bio_sex)
dataout$bio_sex <- relevel(dataout$bio_sex,ref = "male")
n <- nrow(dataout)

## Fitting models: M0
J <- 1; H <- 1; X <- matrix(1,n,J); Z <- matrix(1,n,H)

mu.theta <- c(rep(0,J),rep(0,H)); sigma.2.theta <- c(rep(6.5,J),rep(6.5,H))
datastan <- list(n=n,J=J,H=H,K=1185,lb=0,ub=1185,c=dataout$y_c,h=dataout$y_1.h,X=X,Z=Z,normaliz=dataout$normaliz_factors,
                 mu_beta=mu.theta[1:J],Sigma_beta=sigma.2.theta[1:J]*diag(J),
                 mu_gamma=mu.theta[(J+1):(J+H)],Sigma_gamma=sigma.2.theta[(J+1):(J+H)]*diag(H))
mod_fit <- mod$sample(data = datastan,chains = 12,parallel_chains = 12,iter_sampling = 2.0e3,iter_warmup = 1.5e3)

mod_fit$diagnostic_summary()
mod_fit$print(variables = c("betas","gammas","alpha_h","beta_h"))

log_lik_joint <- mod_fit$draws("log_lik_joint", format = "matrix")
loo_out <- loo::loo(x = log_lik_joint) #PSIS-LOO-CV
waic_out <- loo::waic(log_lik_joint)
summary(loo_out$diagnostics$pareto_k)

save(mod_fit,log_lik_joint,loo_out,waic_out,file = "out/HAS3_mod0.rds")


## Fitting models: M1
X <- model.matrix(~dataout$hba1c_categ); Z <- matrix(1,n,1); J <- ncol(X); H <- 1

mu.theta <- c(rep(0,J),rep(0,H)); sigma.2.theta <- c(rep(6.5,J),rep(6.5,H))
datastan <- list(n=n,J=J,H=H,K=1185,lb=0,ub=1185,c=dataout$y_c,h=dataout$y_1.h,X=X,Z=Z,normaliz=dataout$normaliz_factors,
                 mu_beta=mu.theta[1:J],Sigma_beta=sigma.2.theta[1:J]*diag(J),
                 mu_gamma=mu.theta[(J+1):(J+H)],Sigma_gamma=sigma.2.theta[(J+1):(J+H)]*diag(H))
mod_fit <- mod$sample(data = datastan,chains = 12,parallel_chains = 12,iter_sampling = 2.0e3,iter_warmup = 1.5e3)

mod_fit$diagnostic_summary()
mod_fit$print(variables = c("betas","gammas","alpha_h","beta_h"))

log_lik_joint <- mod_fit$draws("log_lik_joint", format = "matrix")
loo_out <- loo::loo(x = log_lik_joint) #PSIS-LOO-CV
waic_out <- loo::waic(log_lik_joint)
summary(loo_out$diagnostics$pareto_k)

save(mod_fit,log_lik_joint,loo_out,waic_out,file = "out/HAS3_mod1.rds")


## Fitting models: M2
bmi <- scale(dataout$bmi); age <- scale(dataout$age)
X <- model.matrix(~dataout$hba1c_categ+dataout$bio_sex+bmi+age); Z <- matrix(1,n,1); J <- ncol(X); H <- 1

mu.theta <- c(rep(0,J),rep(0,H)); sigma.2.theta <- c(rep(6.5,J),rep(6.5,H))
datastan <- list(n=n,J=J,H=H,K=1185,lb=0,ub=1185,c=dataout$y_c,h=dataout$y_1.h,X=X,Z=Z,normaliz=dataout$normaliz_factors,
                 mu_beta=mu.theta[1:J],Sigma_beta=sigma.2.theta[1:J]*diag(J),
                 mu_gamma=mu.theta[(J+1):(J+H)],Sigma_gamma=sigma.2.theta[(J+1):(J+H)]*diag(H))
mod_fit <- mod$sample(data = datastan,chains = 12,parallel_chains = 12,iter_sampling = 2.0e3,iter_warmup = 1.5e3)

mod_fit$diagnostic_summary()
mod_fit$print(variables = c("betas","gammas","alpha_h","beta_h"))

log_lik_joint <- mod_fit$draws("log_lik_joint", format = "matrix")
loo_out <- loo::loo(x = log_lik_joint) #PSIS-LOO-CV
waic_out <- loo::waic(log_lik_joint)
summary(loo_out$diagnostics$pareto_k)

save(mod_fit,log_lik_joint,loo_out,waic_out,file = "out/HAS3_mod2.rds")

# extracting posterior quantities for next analysis
gen_qts <- list("y"=mod_fit$draws(variables = c("y_rep")),"h"=mod_fit$draws(variables = c("h_rep")),"c"=mod_fit$draws(variables = c("c_rep")))
save(gen_qts,file = "out/HAS3_CNAR_ppc.rds")


## Fitting models: M3
bmi <- scale(dataout$bmi); age <- scale(dataout$age)
X <- model.matrix(~dataout$hba1c_categ+dataout$bio_sex+bmi+age+dataout$hba1c_categ:dataout$bio_sex); Z <- matrix(1,n,1); J <- ncol(X); H <- 1

mu.theta <- c(rep(0,J),rep(0,H)); sigma.2.theta <- c(rep(6.5,J),rep(6.5,H))
datastan <- list(n=n,J=J,H=H,K=1185,lb=0,ub=1185,c=dataout$y_c,h=dataout$y_1.h,X=X,Z=Z,normaliz=dataout$normaliz_factors,
                 mu_beta=mu.theta[1:J],Sigma_beta=sigma.2.theta[1:J]*diag(J),
                 mu_gamma=mu.theta[(J+1):(J+H)],Sigma_gamma=sigma.2.theta[(J+1):(J+H)]*diag(H))
mod_fit <- mod$sample(data = datastan,chains = 12,parallel_chains = 12,iter_sampling = 2.0e3,iter_warmup = 1.5e3)

mod_fit$diagnostic_summary()
mod_fit$print(variables = c("betas","gammas","alpha_h","beta_h"),max_rows = 20)

log_lik_joint <- mod_fit$draws("log_lik_joint", format = "matrix")
loo_out <- loo::loo(x = log_lik_joint) #PSIS-LOO-CV
waic_out <- loo::waic(log_lik_joint)
summary(loo_out$diagnostics$pareto_k)

save(mod_fit,log_lik_joint,loo_out,waic_out,file = "out/HAS3_mod3.rds")




# CNAR model comparison ---------------------------------------------------
mc_table <- matrix(NA,4,8); colnames(mc_table) <- c("ELPD_LOO","SE_ELPD_LOO","P_LOO","WAIC","SE_WAIC","Pareto_K_min","Pareto_K_mean","Pareto_K_max"); rownames(mc_table) <- paste0("M",0:3)
pt_table <- matrix(NA,4,6); colnames(pt_table) <- c("min","Q1","med","mean","Q3","max"); rownames(pt_table) <- paste0("M",0:3)

mod_fits <- NULL
for(i in 0:3){
  load(file = paste0("out/HAS3_mod",i,".rds"))
  mod_fits[[i+1]] <- list(log_lik_joint,loo_out,waic_out,mod_fit)
  pt_table[i+1,] <- summary(loo_out$diagnostics$pareto_k)
  mc_table[i+1,] <- c(mod_fits[[i+1]][[2]]$estimates[1,1:2],mod_fits[[i+1]][[2]]$estimates[2,1],mod_fits[[i+1]][[3]]$estimates[3,],summary(mod_fits[[i+1]][[2]]$diagnostics$pareto_k)[c(1,4,6)])
}
print(mc_table)
print(pt_table)

mhat <- 1 #load fitted M1
load(file = paste0("out/HAS3_mod",mhat,".rds")) 


## Figure S1
xsup <- seq(0,1185,by=20)
pdf(file = "out/fig1.pdf",width = 14,height = 6); par(mfrow=c(1,3))

Y_tilde <- dataout[dataout$hba1c_categ=="normal",c(7,8)]
i <- 1; plot(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),col="#ADD8E6",type="b",pch=20,bty="n",lwd=2,xlim=c(0,1185),ylab="",xlab="",cex.lab=2,cex.axis=2)
for(i in 2:nrow(Y_tilde)){points(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),col="#ADD8E6",type="b",pch=20,lwd=2)}
points(xsup,beta_type_FS(xsup,mean(Y_tilde[,1]),mean(Y_tilde[,2]),max(xsup)),col="#4A708B",type="b",pch=15,lwd=2,cex=2)
title("(A) HbA1c: Normal",adj=0,cex.main=2);

Y_tilde <- dataout[dataout$hba1c_categ=="pre_diabetes",c(7,8)]
i <- 1; plot(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),col="#ADD8E6",type="b",pch=20,bty="n",lwd=2,xlim=c(0,1185),ylab="",xlab="",cex.lab=2,cex.axis=2)
for(i in 2:nrow(Y_tilde)){points(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),col="#ADD8E6",type="b",pch=20,lwd=2)}
points(xsup,beta_type_FS(xsup,mean(Y_tilde[,1]),mean(Y_tilde[,2]),max(xsup)),col="#CD7054",type="b",pch=15,lwd=2,cex=2)
title("(B) HbA1c: Pre-diabetes",adj=0,cex.main=2)

Y_tilde <- dataout[dataout$hba1c_categ=="diabetes",c(7,8)]
i <- 1; plot(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),col="#ADD8E6",type="b",pch=20,bty="n",lwd=2,xlim=c(0,1185),ylab="",xlab="",cex.lab=2,cex.axis=2)
for(i in 2:nrow(Y_tilde)){points(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),col="#ADD8E6",type="b",pch=20,lwd=2)}
points(xsup,beta_type_FS(xsup,mean(Y_tilde[,1]),mean(Y_tilde[,2]),max(xsup)),col="#698B22",type="b",pch=15,lwd=2,cex=2)
title("(C) HbA1c: Diabetes",adj=0,cex.main=2)

dev.off()


## Figure S2
pdf(file = "out/fig2.pdf",width = 14,height = 6); par(mfrow=c(1,2))
set.seed(260318)
par(mfrow=c(1,2));
plot(0,0,col="white",xlim=c(0,3.5),ylim=c(0,1190),bty="n",xlab="",ylab="",axes=FALSE);title("(A)",adj=0,cex.main=2)
axis(side = 1,at = seq(0,3.5,length=5),labels = c("","N","PD","D",""),tick = FALSE,cex.axis=2,cex.lab=2); axis(side = 2,at = seq(0,1000,by=100),labels = seq(0,1000,by=100)/1000,cex.axis=2,cex.lab=2)
Y_tilde <- dataout[dataout$hba1c_categ=="normal",c(7,8)]
A_cuts <- matrix(NA,nrow(Y_tilde),2); 
for(i in 1:nrow(Y_tilde)){
  A_cuts[i,] <- range(alpha_cut(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),alpha = 0.01))
  rect(ybottom = A_cuts[i,1],ytop = A_cuts[i,2],xleft = 0+rnorm(1,0.5,0.1),xright = 0.5+rnorm(1,0.5,0.1),lty=2,border="#4A708B",lwd=2)  
}
Y_tilde <- dataout[dataout$hba1c_categ=="pre_diabetes",c(7,8)]
for(i in 1:nrow(Y_tilde)){
  A_cuts[i,] <- range(alpha_cut(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),alpha = 0.01))
  rect(ybottom = A_cuts[i,1],ytop = A_cuts[i,2],xleft = 1.0+rnorm(1,0.5,0.1),xright = 1.5+rnorm(1,0.5,0.1),border="#CD7054",lwd=2,lty=2)  
}
Y_tilde <- dataout[dataout$hba1c_categ=="diabetes",c(7,8)]
for(i in 1:nrow(Y_tilde)){
  A_cuts[i,] <- range(alpha_cut(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),alpha = 0.01))
  rect(ybottom = A_cuts[i,1],ytop = A_cuts[i,2],xleft = 2+rnorm(1,0.5,0.1),xright = 2.5+rnorm(1,0.5,0.1),border="#698B22",lwd=2,lty=2)  
}

plot(0,0,col="white",xlim=c(0,3.5),ylim=c(0,1190),bty="n",xlab="",ylab="",axes=FALSE);title("(B)",adj=0,cex.main=2)
axis(side = 1,at = seq(0,3.5,length=5),labels = c("","N","PD","D",""),tick = FALSE,cex.axis=2,cex.lab=2); axis(side = 2,at = seq(0,1000,by=100),labels = seq(0,1000,by=100)/1000,cex.axis=2,cex.lab=2)
Y_tilde <- dataout[dataout$hba1c_categ=="normal",c(7,8)]
A_cuts <- matrix(NA,nrow(Y_tilde),2); 
for(i in 1:nrow(Y_tilde)){
  A_cuts[i,] <- range(alpha_cut(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),alpha = 0.95))
  rect(ybottom = A_cuts[i,1],ytop = A_cuts[i,2],xleft = 0+rnorm(1,0.5,0.1),xright = 0.5+rnorm(1,0.5,0.1),lty=2,border="#4A708B",lwd=2)  
}
Y_tilde <- dataout[dataout$hba1c_categ=="pre_diabetes",c(7,8)]
for(i in 1:nrow(Y_tilde)){
  A_cuts[i,] <- range(alpha_cut(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),alpha = 0.95))
  rect(ybottom = A_cuts[i,1],ytop = A_cuts[i,2],xleft = 1.0+rnorm(1,0.5,0.1),xright = 1.5+rnorm(1,0.5,0.1),border="#CD7054",lwd=2,lty=2)  
}
Y_tilde <- dataout[dataout$hba1c_categ=="diabetes",c(7,8)]
for(i in 1:nrow(Y_tilde)){
  A_cuts[i,] <- range(alpha_cut(xsup,beta_type_FS(xsup,Y_tilde[i,1],Y_tilde[i,2],max(xsup)),alpha = 0.95))
  rect(ybottom = A_cuts[i,1],ytop = A_cuts[i,2],xleft = 2+rnorm(1,0.5,0.1),xright = 2.5+rnorm(1,0.5,0.1),border="#698B22",lwd=2,lty=2)  
}
dev.off()


## Table S1
Xtab_tex = xtable::xtable(mc_table)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = "tab1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

## Table S2 (CNAR subtable)
draws_array <- as_draws_matrix(mod_fit$draws(variables = c("betas", "gammas", "alpha_h", "beta_h")))
draws_array[,4] <- exp(draws_array[,4])
post_table <- cbind(mod_fit$summary(variables = c("betas","gammas","alpha_h","beta_h"))[,c(2,4,9,10)],
                    coda::HPDinterval(as.mcmc(draws_array)))
post_table <- post_table[,c(1,2,5,6,3,4)]
post_table

Xtab_tex = xtable::xtable(post_table)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = "tab1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})
### !!! ###


## Table S3
Xtab_tex = xtable::xtable(pt_table)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = "tab1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


## Figure S4
posts <- cbind(drop(merge_chains(x = mod_fit$draws(variables = "betas"))),
               exp(drop(merge_chains(x = mod_fit$draws(variables = "gammas")))),
               drop(merge_chains(x = mod_fit$draws(variables = "alpha_h"))),
               drop(merge_chains(x = mod_fit$draws(variables = "beta_h"))))

pdf(file = "out/fig3.pdf",width = 14,height = 6)
par(mfcol=c(2,6))
for(j in 1:3){
  plot(posts[,j],bty="n",xlab="",ylab="",type="l",xlim=c(0,24e3),col="#8B8989");abline(h = mean(posts[,j]),lty=3,col="#CDCD00",lwd=2.5); title(bquote(beta[.(j-1)]),cex.main=2.5)
  hist(posts[,j],col="#8B8989",border = "white",xlab="",ylab="",main="",nclass = 40)
}
j <- 4; plot(posts[,j],bty="n",xlab="",ylab="",type="l",xlim=c(0,24e3),col="#8B8989");abline(h = mean(posts[,j]),lty=3,col="#CDCD00",lwd=2.5); title(bquote(kappa),cex.main=2.5)
hist(posts[,j],col="#8B8989",border = "white",xlab="",ylab="",main="",nclass = 40)
j <- 5; plot(posts[,j],bty="n",xlab="",ylab="",type="l",xlim=c(0,24e3),col="#8B8989");abline(h = mean(posts[,j]),lty=3,col="#CDCD00",lwd=2.5); title(bquote(alpha[s]),cex.main=2.5)
hist(posts[,j],col="#8B8989",border = "white",xlab="",ylab="",main="",nclass = 40)
j <- 6; plot(posts[,j],bty="n",xlab="",ylab="",type="l",xlim=c(0,24e3),col="#8B8989");abline(h = mean(posts[,j]),lty=3,col="#CDCD00",lwd=2.5); title(bquote(beta[s]),cex.main=2.5)
hist(posts[,j],col="#8B8989",border = "white",xlab="",ylab="",main="",nclass = 40)
dev.off()





# CNAR vs CAR-like models -------------------------------------------------

## Fitting CAR-like models using M1 specification

# CAR-like 1 
mod <- cmdstanr::cmdstan_model(stan_file = "CAR_1.stan")

X <- model.matrix(~dataout$hba1c_categ); Z <- matrix(1,n,1); J <- ncol(X); H <- 1
mu.theta <- c(rep(0,J),rep(0,H)); sigma.2.theta <- c(rep(6.5,J),rep(6.5,H))
datastan <- list(n=n,J=J,H=H,K=1185,lb=0,ub=1185,c=dataout$y_c,h=dataout$y_1.h,X=X,Z=Z,normaliz=dataout$normaliz_factors,
                 mu_beta=mu.theta[1:J],Sigma_beta=sigma.2.theta[1:J]*diag(J))
mod_fit <- mod$sample(data = datastan,chains = 12,parallel_chains = 12,iter_sampling = 1.0e3,iter_warmup = 1.5e3)

mod_fit$diagnostic_summary()
gen_qts <- list("h"=mod_fit$draws(variables = c("h_rep")),"c"=mod_fit$draws(variables = c("c_rep")))
save(gen_qts,file = "out/HAS3_CAR1_ppc.rds")

# CAR-like 2
mod <- cmdstanr::cmdstan_model(stan_file = "../CAR_2.stan")

mod_fit <- mod$sample(data = datastan,chains = 12,parallel_chains = 12,iter_sampling = 1.0e3,iter_warmup = 1.5e3)

mod_fit$diagnostic_summary()
gen_qts <- list("h"=mod_fit$draws(variables = c("h_rep")),"c"=mod_fit$draws(variables = c("c_rep")))
save(gen_qts,file = "out/HAS3_CAR2_ppc.rds")


## PPC analysis (basic)

# loading posterior draws
load("HAS3/HAS3_CNAR_ppc.rds"); cnar_posts <- gen_qts 
load("HAS3/HAS3_CAR1_ppc.rds"); car1_posts <- gen_qts
load("HAS3/HAS3_CAR2_ppc.rds"); car2_posts <- gen_qts

C_cnar <- matrix(cnar_posts$c,ncol = n); H_cnar <- matrix(cnar_posts$h,ncol = n)
C_car1 <- matrix(car1_posts$c,ncol = n); H_car1 <- matrix(car1_posts$h,ncol = n)
C_car2 <- matrix(car2_posts$c,ncol = n); H_car2 <- matrix(car2_posts$h,ncol = n)

z1 <- MASS::kde2d(as.vector(C_cnar),as.vector(H_cnar),n = 100)
z2 <- MASS::kde2d(as.vector(C_car1),as.vector(H_car1),n = 100)
z3 <- MASS::kde2d(as.vector(C_car2),as.vector(H_car2),n = 100)


## Figure 1a
pdf(file = "out/fig1a.pdf",width = 16,height = 6)
ymax <- max(dataout$y_1.h)+median(dataout$y_1.h)*2
par(mfrow=c(1,2))
plot(dataout$y_c,dataout$y_1.h,pch=20,cex=1.5,col="white",bty="n",xlab="",ylab="",ylim=c(0,ymax),cex.axis=1.5,cex.lab=2); title(expression(CAR[like]^1 ~vs~ CNAR),adj=0,line=2,cex.main=2)
contour(z1,col = "#CDCDC1",lwd=2.25,add=TRUE,nlevels = 20)
contour(z2,col = hcl.colors(10, "YlGnBu"),lwd=1.5,add=TRUE,nlevels = 20)
points(dataout$y_c,dataout$y_1.h,pch=20,cex=1.5,col="#CD8162")

plot(dataout$y_c,dataout$y_1.h,pch=20,cex=1.5,col="white",bty="n",xlab="",ylab="",ylim=c(0,ymax),cex.axis=1.5,cex.lab=2); title(expression(CAR[like]^2 ~vs~ CNAR),adj=0,line=2,cex.main=2)
contour(z1,col = "#CDCDC1",lwd=2.25,add=TRUE,nlevels = 20)
contour(z3,col = hcl.colors(10, "YlGnBu"),lwd=1.5,add=TRUE,nlevels = 20)
points(dataout$y_c,dataout$y_1.h,pch=20,cex=1.5,col="#CD8162")
dev.off()

## Table S5
# bayesian p-vs (transformed)
c_means <- c(mean(apply(C_cnar,1,mean)>=mean(dataout$y_c)),mean(apply(C_car1,1,mean)>=mean(dataout$y_c)),mean(apply(C_car2,1,mean)>=mean(dataout$y_c)))
c_vars <- c(mean(apply(C_cnar,1,var)>=var(dataout$y_c)),mean(apply(C_car1,1,var)>=var(dataout$y_c)),mean(apply(C_car2,1,var)>=var(dataout$y_c)))
c_means <- sapply(c_means,function(x)abs(0.5-x)); c_vars <- sapply(c_vars,function(x)abs(0.5-x))
# coverages (90%)
c_cvgs <- c(mean(mapply(function(i){q_i<-as.numeric(quantile(C_cnar[,i],c(0.025,0.95)));dataout$y_c[i]>=q_i[1] & dataout$y_c[i]<=q_i[2]},1:n)),
            mean(mapply(function(i){q_i<-as.numeric(quantile(C_car1[,i],c(0.025,0.95)));dataout$y_c[i]>=q_i[1] & dataout$y_c[i]<=q_i[2]},1:n)),
            mean(mapply(function(i){q_i<-as.numeric(quantile(C_car2[,i],c(0.025,0.95)));dataout$y_c[i]>=q_i[1] & dataout$y_c[i]<=q_i[2]},1:n)))
A <- round(cbind(c_means,c_vars,c_cvgs),3); colnames(A) <- c("mean","var","cvg"); rownames(A) <- c("CNAR","CAR1","CAR2")

Xtab_tex = xtable::xtable(A)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = "tab1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


## PPC analysis (energy-like)

# retaing a subset of Bb draws from posteriors
set.seed(200326)
B <- dim(C_cnar)[1]; Bb <- 1000
bbd <- sample(1:B,size = Bb,replace = FALSE)
C_post_0 <- C_cnar[bbd,]; H_post_0 <- 1/H_cnar[bbd,]
C_post_1 <- C_car1[bbd,]; H_post_1 <- 1/H_car1[bbd,]
C_post_2 <- C_car2[bbd,]; H_post_2 <- 1/H_car2[bbd,]
for(i in 1:n){
  H_post_0[H_post_0[,i]<exp(-4.5),i] <- exp(-4.5); # the empirical precision for a beta-type fuzzy set cannot go below exp(-4.5)
  H_post_1[H_post_1[,i]<exp(-4.5),i] <- exp(-4.5);
  H_post_2[H_post_2[,i]<exp(-4.5),i] <- exp(-4.5)
}
save(C_post_0,H_post_0,C_post_1,H_post_1,C_post_2,H_post_2,file = "out/PPC_Bb.rds")

# compute Bertoluzza distance on observed fuzzy counts (u_obs)
DB_obs <- matrix(NA,n,n)
for(i in 1:n){
  cat(i,"\t")
  muy_i <- beta_type_FS(x = 0:1185,dataout$y_c[i],dataout$y_h[i],k = 1185)
  for(j in i:n){
    muy_j <- beta_type_FS(x = 0:1185,dataout$y_c[j],dataout$y_h[j],k = 1185)
    DB_obs[i,j] <- DB_obs[j,i] <- bertoluzza_dist(0:1185,muy_i,muy_j)
  }
}

# compute u_cross and u_reps using a multi-core strategy on a local Linux-based computer
# (it requires GNU Parallel to be properly installed)
ncores <- parallel::detectCores()-1
id_core <- sort(rep(1:ncores, length.out = Bb))
jobs_data <- split(1:Bb, id_core)

jobs <- sapply(1:ncores,function(u){paste0("R CMD BATCH --no-restore --no-save '--args coreId=",u," bbd=c(",paste0(jobs_data[[u]],collapse=","),")' PPC_Bb.R out/out_coreId_",u,".stdout")})
batch_file <- tempfile(tmpdir = "PPC","job_list_"); writeLines(jobs, batch_file)
system(paste0("parallel -j ", ncores, " < ", batch_file," &"));

# collect files and compute desired quantities
DB_reps_0 <- DB_cross_0 <- vector(mode = "list",length = 0); DB_reps_1 <- DB_cross_1 <- vector(mode = "list",length = 0); DB_reps_2 <- DB_cross_2 <- vector(mode = "list",length = 0)
for(u in 1:ncores){
  cat(u,"\t")
  load(paste0("out/PPC_Bb_",u,".rds"))
  
  DB_reps_0 <- c(DB_reps_0,Filter(Negate(is.null),DB_reps_list_0)); DB_cross_0 <- c(DB_cross_0,Filter(Negate(is.null),DB_cross_list_0))
  DB_reps_1 <- c(DB_reps_1,Filter(Negate(is.null),DB_reps_list_1)); DB_cross_1 <- c(DB_cross_1,Filter(Negate(is.null),DB_cross_list_1))
  DB_reps_2 <- c(DB_reps_2,Filter(Negate(is.null),DB_reps_list_2)); DB_cross_2 <- c(DB_cross_2,Filter(Negate(is.null),DB_cross_list_2))
  
  rm(DB_cross_list_0,DB_cross_list_1,DB_cross_list_2,DB_reps_list_0,DB_reps_list_1,DB_reps_list_2)
}


bertDist_B <- 1/(n^2) * sum(DB_obs) #bertoluzza's distance among the observed fuzzy counts (n x n)

bertDist_A_0 <- as.numeric(lapply(DB_cross_0,function(x)1/(n^2)*sum(x))) #bertoluzza's distance among observed and replicated fuzzy counts --fixing the rep (Bb x n x n)
bertDist_C_0 <- as.numeric(lapply(DB_reps_0,function(x)1/(n^2)*sum(x))) #bertoluzza's distance among replicated fuzzy counts --fixing the rep (Bb x n x n)
bertDist_A_1 <- as.numeric(lapply(DB_cross_1,function(x)1/(n^2)*sum(x))); bertDist_C_1 <- as.numeric(lapply(DB_reps_1,function(x)1/(n^2)*sum(x)))
bertDist_A_2 <- as.numeric(lapply(DB_cross_2,function(x)1/(n^2)*sum(x))); bertDist_C_2 <- as.numeric(lapply(DB_reps_2,function(x)1/(n^2)*sum(x)))

## Figure 1b
pdf(file = "out/fig1b.pdf",width = 12,height = 6)
par(mfrow=c(1,2))
boxplot((bertDist_C_0),(bertDist_C_1),(bertDist_C_2),outline=FALSE,frame=FALSE,col=c("#00688B","#548B54","#CD3700"),border = "#404040",cex.axis=1.5,boxwex = 0.5,names=c("CNAR",expression(CAR[like]^1),expression(CAR[like]^2)));abline(h = bertDist_B,col="#FFA500",lwd=2,lty=5); title(expression(u[reps]),adj=0,line=1,cex.main=2); 
boxplot((bertDist_A_0),(bertDist_A_1),(bertDist_A_2),outline=FALSE,frame=FALSE,col=c("#00688B","#548B54","#CD3700"),border = "#404040",boxwex = 0.5,names = c("CNAR",expression(CAR[like]^1),expression(CAR[like]^2)),cex.axis=1.5);abline(h = bertDist_B,col="#FFA500",lwd=2,lty=5); title(expression(u[cross]),adj=0,line=1,cex.main=2)
dev.off()



# Defuzz vs RSEM analysis -------------------------------------------------

mod <- cmdstanr::cmdstan_model(stan_file = "naive_model.stan")
y_hat <- round(mapply(function(i)defuzzify_centroid(dataout$y_c[i],si = dataout$y_1.h[i],what = "mean"),1:n))

X <- model.matrix(~dataout$hba1c_categ); Z <- matrix(1,n,1); J <- ncol(X); H <- 1
mu.theta <- c(rep(0,J),rep(0,H)); sigma.2.theta <- c(rep(6.5,J),rep(6.5,H))
datastan <- list(n=n,J=J,H=H,K=1185,y=y_hat,X=X,Z=Z,normaliz=dataout$normaliz_factors,
                 mu_beta=mu.theta[1:J],Sigma_beta=sigma.2.theta[1:J]*diag(J),
                 mu_gamma=mu.theta[(J+1):(J+H)],Sigma_gamma=sigma.2.theta[(J+1):(J+H)]*diag(H))

mod_fit <- mod$sample(data = datastan,chains = 12,parallel_chains = 12,iter_sampling = 1.0e3,iter_warmup = 2.5e3)
mod_fit$diagnostic_summary()

Y_hat_defuzz <- matrix(mod_fit$draws(variables = c("y_rep")),ncol = n)
draws_array <- as_draws_matrix(mod_fit$draws(variables = c("betas", "gammas")))
draws_array[,4] <- exp(draws_array[,4])
post_table <- cbind(mod_fit$summary(variables = c("betas","phi[1]"))[,c(2,4,9,10)],
                    coda::HPDinterval(as.mcmc(draws_array)))
post_table <- post_table[,c(1,2,5,6,3,4)]
post_table

## Table S2 (defuzz subtable)
Xtab_tex = xtable::xtable(post_table)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = "tab1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


load(file = "data/rsem_exp_counts.rds")
has3_rsem <- expCnts[grep(x = expCnts$gene_id,pattern = "HAS3"),][-1]
has3_rsem <- as.numeric(has3_rsem[iid])

datastan <- list(n=n,J=J,H=H,K=1185,y=has3_rsem,X=X,Z=Z,normaliz=dataout$normaliz_factors,
                 mu_beta=mu.theta[1:J],Sigma_beta=sigma.2.theta[1:J]*diag(J),
                 mu_gamma=mu.theta[(J+1):(J+H)],Sigma_gamma=sigma.2.theta[(J+1):(J+H)]*diag(H))

mod_fit <- mod$sample(data = datastan,chains = 12,parallel_chains = 12,iter_sampling = 1.0e3,iter_warmup = 2.5e3)
mod_fit$diagnostic_summary()

Y_hat_rsem <- matrix(mod_fit$draws(variables = c("y_rep")),ncol = n)

draws_array <- as_draws_matrix(mod_fit$draws(variables = c("betas", "gammas")))
draws_array[,4] <- exp(draws_array[,4])
post_table <- cbind(mod_fit$summary(variables = c("betas","phi[1]"))[,c(2,4,9,10)],
                    coda::HPDinterval(as.mcmc(draws_array)))
post_table <- post_table[,c(1,2,5,6,3,4)]
post_table


## Table S2 (RSEM subtable)
Xtab_tex = xtable::xtable(post_table)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = "tab1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


## Figure S5
pdf(file = "../../SPL/suppmat/fig6.pdf",width = 8,height = 7)
hist(has3_rsem,nclass = 20,border = "gray",col = adjustcolor("#8EE5EE", alpha.f = 0.5),main="",xlab="",ylab="",ylim=c(0,25),xlim=c(0,1e3))
hist(y_hat,nclass = 20,border = "gray",col = adjustcolor("#FFD39B", alpha.f = 0.5),main="",xlab="",ylab="",add=TRUE)
legend(x = 600, y = 20,border = "transparent",legend = c("RSEM","DEFUZZ"),col = c("#8EE5EE","#FFD39B"),ncol = 1,box.col = "white",
       bg = "transparent",x.intersp = 0.5,cex = 1.5,pt.cex=2,pch=15)
dev.off()


set.seed(230326)
B <- min(NROW(Y_hat_defuzz),NROW(Y_hat_rsem))
K <- 1185; bbd <- sample(x = 1:B,size = 2e3,replace = FALSE)
C_hat_defuzz <- t(mapply(function(b){mapply(function(i)rbeta(n = 1,shape1 = (1/K*Y_hat_defuzz[b,i])*dataout$y_1.h[i],shape2 = dataout$y_1.h[i]-(1/K*Y_hat_defuzz[b,i])*dataout$y_1.h[i]),1:n)*K},bbd))
C_hat_rsem <- t(mapply(function(b){mapply(function(i)rbeta(n = 1,shape1 = (1/K*Y_hat_rsem[b,i])*dataout$y_1.h[i],shape2 = dataout$y_1.h[i]-(1/K*Y_hat_rsem[b,i])*dataout$y_1.h[i]),1:n)*K},bbd))  


## Figure S6
pdf(file = "../../SPL/suppmat/fig7.pdf",width = 8,height = 7)  
pt1 <- (mean(dataout$y_c/dataout$y_1.h)); pt2 <- (diff(quantile(dataout$y_c/dataout$y_1.h,c(0.1,0.9))))
xt1 <- (mapply(function(b)mean(C_hat_defuzz[b,]/dataout$y_1.h),1:length(bbd))); yt1 <- (mapply(function(b)diff(quantile(C_hat_defuzz[b,]/dataout$y_1.h,c(0.1,0.9))),1:length(bbd)))
xt2 <- (mapply(function(b)mean(C_hat_rsem[b,]/dataout$y_1.h),1:length(bbd))); yt2 <- (mapply(function(b)diff(quantile(C_hat_rsem[b,]/dataout$y_1.h,c(0.1,0.9))),1:length(bbd)))
xlims <- c(min(min(xt1),min(xt2)),max(max(xt1),max(xt2))); ylims <- c(min(min(yt1),min(yt2)),max(max(yt1),max(yt2)));
plot(xt1,yt1,bty="n",xlab="",ylab="",ylim=ylims,xlim=xlims,col="#CDCDC1",pch=20,cex=3); points(pt1,pt2,pch=20,col=2,cex=4); 
points(xt2,yt2,bty="n",ylim=ylims,xlim=xlims,col=adjustcolor("#548B54", alpha.f = 0.5),pch=20,cex=2); 
legend(x = 10, y = 12,border = "transparent",legend = c("Defuzz","Rsem"),col = c("#CDCDC1","#548B54"),ncol = 1,box.col = "white",bg = "transparent",
       x.intersp = 0.5,cex = 1.5,pt.cex=2,pch=15)
dev.off()



# Case-study-based simulation ---------------------------------------------

#To appropriately run the case-study simulation (see the folder M1_short_simul): 
# run gendata.R to populate simdata/
# run run_sim.sh on a proper HPC system equipped with Slurm and Singularuity 
# the folder fitdata/ is populated with the model estimates 

# Assuming the folder fitdata/ is populated appropriately, run the following code:
fls <- list.files(path = "M1_short_simulation/fitdata/",full.names = TRUE,pattern = ".rds")
Out <- array(NA,c(length(fls),4,2),dimnames = list(NULL,c("mean","err","cvg","leng"),c("CNAR","Naive"))) 

for(i in 1:length(fls)){
  cat(i,"\t")
  load_status <- try(load(fls[i]), silent = TRUE)
  if(!inherits(load_status, "try-error")){
    current_iid <- unlist(strsplit(x = fls[i],split = "fitdata_|.rds"))[2]
    load(paste0("M1_short_simulation/simdata/simdata_",current_iid,".rds")) #it contains theta_m needed to compute mu0 and get kappa0 (theta_m[4])
    
    if(!is.null(naive_draws)){
      
      p_cnar_k <- matrix(posterior::merge_chains(cnar_draws[,,5]))
      p_naive_k <- matrix(posterior::merge_chains(naive_draws[,,5]))
      
      Bb <- min(nrow(p_cnar_k),nrow(p_naive_k))
      p_cnar_k <- p_cnar_k[1:Bb,]; p_naive_k <- p_naive_k[1:Bb,]
      
      Out[i,1,1:2] <- c(mean(p_cnar_k),mean(p_naive_k))
      Out[i,2,1:2] <- apply(cbind(p_cnar_k,p_naive_k)-theta_m[4],2,mean)
      qm1 <- quantile(p_cnar_k,c(0.025,0.975)); qm2 <- quantile(p_naive_k,c(0.025,0.975)); 
      Out[i,3,1:2] <- c(as.numeric(theta_m[4]>qm1[1]&theta_m[4]<qm1[2]),
                        as.numeric(theta_m[4]>qm2[1]&theta_m[4]<qm2[2]))
      Out[i,4,1:2] <- c(as.numeric(diff(qm1)),as.numeric(diff(qm2)))
      
      mu0 <- exp(X%*%theta_m[1:3])
      var0 <- mu0 + mu0^2 / exp(theta_m[4])
      
      p_cnar_b <- matrix(posterior::merge_chains(cnar_draws[,,2:4]),ncol = 3)
      p_naive_b <- matrix(posterior::merge_chains(naive_draws[,,2:4]),ncol = 3)
      p_cnar_b <- p_cnar_b[1:Bb,]; p_naive_b <- p_naive_b[1:Bb,]
    }
  }
  
}
save(Out,file = "M1_short_simulation/data_res.Rdata")

A <- cbind(apply(Out[,2,],2,mean,na.rm=TRUE),
           apply(Out[,3,],2,mean,na.rm=TRUE),
           apply(Out[,4,],2,mean,na.rm=TRUE))
colnames(A) <- c("bias","0.95 coverage","0.95 interval length")
rownames(A) <- c("CNAR","defuzz")

## Table S4
Xtab_tex = xtable::xtable(A)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = "tab1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})



