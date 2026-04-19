#args <- commandArgs(trailingOnly = TRUE); for (arg in args) {eval(parse(text = arg))} #to retrieve input args
args <- commandArgs(trailingOnly = TRUE)
print(args)
currdir <- args[2]
mmd <- as.numeric(unlist(strsplit(gsub("[c()]", "", args[3]) , ",")))


require(cmdstanr)
set_cmdstan_path("../cmdstan-2.38.0"); cmdstan_version()

load(file = paste0(currdir,"data_common.rds"))

mod_cnar <- cmdstanr::cmdstan_model(stan_file = "CNAR_negbin.stan")
mod_naive <- cmdstanr::cmdstan_model(stan_file = "naive_model.stan")

iter_sampling <- 4e3; iter_warmup <- 1.5e3
parallel_chains <- 1
num_chains <- 2
K_chain <- 2 

for(m in mmd){
cat("\n\t Sample:",m)

load(file = paste0(currdir,"simdata/simdata_",m,".rds"))
mod_cnar_fit <- mod_naive_fit <- NULL

#CNAR
cat("\n CNAR \n")
datastan <- list(n=n,J=J,H=H,K=1185,lb=0,ub=1185,c=c,h=h,X=X,Z=Z,normaliz=u,mu_beta=mu.theta[1:J],Sigma_beta=sigma.2.theta[1:J]*diag(J),mu_gamma=mu.theta[(J+1):(J+H)],Sigma_gamma=sigma.2.theta[(J+1):(J+H)]*diag(H))
mod_cnar_fit <- mod_cnar$sample(data = datastan,chains = num_chains,parallel_chains = parallel_chains,iter_sampling = iter_sampling,iter_warmup = iter_warmup)
cnar_draws <- tryCatch(mod_cnar_fit$draws(), error = function(e) NULL); cnar_summary <- tryCatch(mod_cnar_fit$summary(), error = function(e) NULL)

#Naive
cat("\n Naive \n")
datastan <- list(n=n,J=J,H=H,K=1185,y=round(cd),X=X,Z=Z,normaliz=u,mu_beta=mu.theta[1:J],Sigma_beta=sigma.2.theta[1:J]*diag(J),mu_gamma=mu.theta[(J+1):(J+H)],Sigma_gamma=sigma.2.theta[(J+1):(J+H)]*diag(H))
mod_naive_fit <- mod_naive$sample(data = datastan,chains = K_chain*num_chains,parallel_chains = parallel_chains,iter_sampling = iter_sampling,iter_warmup = iter_warmup)
naive_draws <- tryCatch(mod_naive_fit$draws(), error = function(e) NULL); naive_summary <- tryCatch(mod_naive_fit$summary(), error = function(e) NULL)

save(cnar_draws,naive_draws,cnar_summary,naive_summary,file = paste0(currdir,"fitdata/fitdata_",m,".rds"))

}

