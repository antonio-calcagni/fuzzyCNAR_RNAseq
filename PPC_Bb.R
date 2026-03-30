args <- commandArgs(trailingOnly = TRUE); for (arg in args) {eval(parse(text = arg))} #to retrieve input args

load("HAS3/PPC_Bb.rds")
dataout <- read.csv(file = "HAS3/HAS3_complete.csv",header = TRUE,dec = ",")
source("utils.R")

n <- dim(C_post_0)[2]
Bb <- length(bbd) #It is computed core-based
DB_cross_list_0 <- DB_cross_list_1 <- DB_cross_list_2 <- vector(mode = "list",length = length(bbd))
DB_reps_list_0 <- DB_reps_list_1 <- DB_reps_list_2 <- vector(mode = "list",length = length(bbd))

for(b in bbd){
  #cat(b,"\t")
  DB0 <- DB1 <- DB2 <- matrix(NA,n,n); DBB0 <- DBB1 <- DBB2 <- matrix(NA,n,n)
  
  for(i in 1:n){
    cat(b,i,"\t")
    muy_i <- beta_type_FS(x = 0:1185,dataout$y_c[i],dataout$y_h[i],k = 1185)
    
    muy_0_i <- beta_type_FS(x = 0:1185,C_post_0[b,i],H_post_0[b,i],k = 1185)
    muy_1_i <- beta_type_FS(x = 0:1185,C_post_1[b,i],H_post_1[b,i],k = 1185)
    muy_2_i <- beta_type_FS(x = 0:1185,C_post_2[b,i],H_post_2[b,i],k = 1185)
    
    for(j in i:n){
      muy_j <- beta_type_FS(x = 0:1185,C_post_0[b,j],H_post_0[b,j],k = 1185)
      DB0[i,j] <- DB0[j,i] <- bertoluzza_dist(0:1185,muy_i,muy_j)
      DBB0[i,j] <- DBB0[j,i] <- bertoluzza_dist(0:1185,muy_0_i,muy_j)
      
      muy_j <- beta_type_FS(x = 0:1185,C_post_1[b,j],H_post_1[b,j],k = 1185)
      DB1[i,j] <- DB1[j,i] <- bertoluzza_dist(0:1185,muy_i,muy_j)
      DBB1[i,j] <- DBB1[j,i] <- bertoluzza_dist(0:1185,muy_1_i,muy_j)
      
      muy_j <- beta_type_FS(x = 0:1185,C_post_2[b,j],H_post_2[b,j],k = 1185)
      DB2[i,j] <- DB2[j,i] <- bertoluzza_dist(0:1185,muy_i,muy_j)
      DBB2[i,j] <- DBB2[j,i] <- bertoluzza_dist(0:1185,muy_2_i,muy_j)
    }
  }
  
  DB_cross_list_0[[b]] <- DB0; DB_reps_list_0[[b]] <- DBB0
  DB_cross_list_1[[b]] <- DB1; DB_reps_list_1[[b]] <- DBB1
  DB_cross_list_2[[b]] <- DB2; DB_reps_list_2[[b]] <- DBB2
}


save(DB_cross_list_0,DB_cross_list_1,DB_cross_list_2,
     DB_reps_list_0,DB_reps_list_1,DB_reps_list_2,file = paste0("HAS3/PPC/PPC_Bb_",coreId,".rds"))

