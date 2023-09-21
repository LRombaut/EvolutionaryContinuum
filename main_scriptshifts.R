rm(list=ls())
library(EasyABC)
source('simulate_LMshifts.R')
target <- readRDS('divergence_time.RDS')
prior <- list(c("unif",0,20),c("unif",0,200),c("unif",0,10^4),c("unif",0,10))

try <- ABC_sequential(method="Lenormand",model=simulate_variance,prior=prior,nb_simul=100,
                               summary_stat_target= target, n_cluster=1,
                               use_seed = TRUE, verbose=FALSE,seed_count=42, 
                               alpha=0.5, progress_bar=TRUE, p_acc_min=0.05)
saveRDS(try,'ABCoutput10_1.RDS')



