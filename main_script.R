rm(list=ls())
library(EasyABC)
source('simulate_LM.R')
target <- readRDS('divergence_time.RDS')
prior <- list(c("unif",0,20),c("unif",0,200),c("unif",0,10^5))

try <- ABC_sequential(method="Lenormand",model=simulate_variance,prior=prior,nb_simul=100,
                               summary_stat_target= target, n_cluster=1,
                               use_seed = TRUE, verbose=FALSE,seed_count=42, 
                               alpha=0.5, progress_bar=TRUE, p_acc_min=0.05)
saveRDS(try,'ABCoutput_10.RDS')


7138+3919+1000+100+13606+3003
