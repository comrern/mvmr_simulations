
args <- as.numeric(commandArgs(T))
set.seed((args[1]*100000))
job_id <- ((args[1]))
message("job number ", job_id)


setwd("/user/work/kb22541/simulations")
output_path <- "./results"
.libPaths("/user/work/kb22541/rlib")



library(dplyr)
library(MASS)
library(TwoSampleMR)
library(MVMR)
library(truncnorm)
library(tidyverse)

source('modes_sims.R')
source('functions_sims.R')

reps = 2
results = data.frame()
results_all = NULL
results_ivw = NULL
results_out = NULL
mvmrres <- NULL

results_all = data.frame()
results_ivw = data.frame()

for (setup_mode in c(1,2,3,4)){
  setup_mode=2
  results_models <- data.frame()
  for (model in c("A","B","C","D"))  {
  results_rep = data.frame()
  run=0
      for(j in 1:reps){  
          
            params <- setup(setup_mode, model)
            snps = params[1]      
            snpsc = params[2]         
            nobs = params[3]
            b1 = params[4]
            b2 = params[5]
            betaC=params[6]
            beta2C=params[7]
            LD_mod=params[9]
            xi=params[8]
            
            
            dat <- data_gen(snps, snpsc, nobs, b1 , b2, betaC, beta2C,  LD_mod)
              #(no of snps, snps for confounding var, samplesize, beta1, beta2, snp-confounder effect)
          #### 
          
            results[1,"model"] <- model
            results[1,"pi"] <- pi
            results[1,"sample.size"] <- nobs

            MR_dat <- GWASres(dat)
          
          #### Univariate two sample MR ####
            
            univariate_results <- univariate_MR(MR_dat)
            
          ### MVMR
            
            mvmr_res <- run_mvmr(MR_dat, dat, setup_mode)
            
        ## format results
          res_run <- rbind(univariate_results, mvmr_res)
          res_run$mode <- model
          
          results_rep <- rbind(results_rep, res_run)  
          
      }
    # results_rep$mod <- mod
  
    results_models <- rbind(results_models, results_rep)  
  }
  results_models$setup_mode <- setup_mode
  results_all <- rbind(results_all, results_models)
  
  

}

## save individual outputs for troubleshooting

save(results_all, file=sprintf(paste0(output_path, "/results_%s.csv"), job_id))


 
 
 
