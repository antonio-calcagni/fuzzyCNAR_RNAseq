ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

#miss <- as.vector(read.csv("missed_data.dat",header=FALSE,as.is=TRUE))[[1]]
#M <- length(miss)

curr_dir <- "rnaseq_fuzzy/M1_short_simulation/"
M <- length(list.files(paste0(curr_dir,"simdata")))
id_core <- sort(rep(1:ncores, length.out = M))
jobs_data <- split(1:M, id_core)
#jobs_data <- split(miss, id_core)

#jobs <- sapply(1:ncores,function(u){paste0("R CMD BATCH --no-restore --no-save '--args coreId=",u," mmd=c(",paste0(jobs_data[[u]],collapse=","),")' /home/calcagni/rnaseq_fuzzy/M1_short_simulation/M1_run_models.R /home/calcagni/rnaseq_fuzzy/M1_short_simulation/out_coreId_",u,".stdout")})
jobs <- sapply(1:ncores,function(u){paste0("R CMD BATCH --no-restore --no-save '--args coreId=",u," ",curr_dir," c(",paste0(jobs_data[[u]],collapse=","),")' rnaseq_fuzzy/M1_short_simulation/M1_run_models.R rnaseq_fuzzy/M1_short_simulation/out_coreId_",u,".stdout")})

batch_file <- tempfile(tmpdir = "/home/calcagni/rnaseq_fuzzy/M1_short_simulation/","job_list_"); writeLines(jobs, batch_file)
writeLines(jobs, batch_file)

system(paste0("parallel -j ", ncores, " < ", batch_file))
