args <- as.numeric(commandArgs(T))
set.seed((args[1]*100000))
job_id <- ((args[1]))
message("job number ", job_id)


setwd("/user/work/kb22541/simulations/ld_experiments")
output_path <- "./results/"
.libPaths("/user/work/kb22541/rlib")



library(dplyr)
library(MASS)
library(TwoSampleMR)
library(MVMR)
library(truncnorm)
library(tidyverse)
library(metafor)

source('modes_sims.R')
source('functions_sims.R')


reps = 1000
run = 0
results = data.frame()
results_all = NULL
results_out = NULL
mvmrres <- NULL



## LD mod settings;

  # 1: LD varies for outcome, but exposure assoc do not vary by ancestry- replicates use of external (e.g Euro) exposure data
  # 2: LD varies equally for exposure and outcome- but is equal for both- represents two sample MR in an admixed population where exposure and outcome share ancestry
  # 3: LD varies for exposure and outcome, but these samples are from separate strata of the population. 





for (LD_mod in c(1,2,3,4)){
  results_models <- data.frame()

      for (model in c("A","B","C","D","E"))  {

    results_rep = data.frame()
    run=0
    for(j in 1:reps){  
      
      setup_mode=1
      run <- run + 1  
      print(paste("on run", run, model, "and LD mode", LD_mod))
      
      params <- setup(setup_mode, model)
      snps = params[1]      
      snpsc = params[2]         
      nobs = params[3]
      b1 = params[4]
      b2 = params[5]
      betaC=params[6]
      beta2C=params[7]
      LD_mag=params[8]
      
                     
      
      dat <- data_gen(snps, snpsc, nobs, b1 , b2, betaC, beta2C, LD_mod, LD_mag)
      #(no of snps, snps for confounding var, samplesize, beta1, beta2, snp-confounder effect)
      #### 
      
      results[1,"model"] <- model
      results[1,"sample.size"] <- nobs
      
      ####Regression check######
      
      ## allele freq
      
      MR_dat <- GWASres(dat, LD_mod)
      
      #### Univariate two sample MR ####
      
      univariate_results <- univariate_MR(MR_dat, LD_mod)
      
      ### MVMR
      
      # mvmr_res <- run_mvmr(MR_dat, dat, LD_mod)
      
      ## format results
      # res_run <- rbind(univariate_results, mvmr_res)
      res_run <- univariate_results
      
      het <- heterogeneity(MR_dat)
      
      res_run$run <- j
      res_run$mode <- model
      res_run$Q <- het[2]
      res_run$Qpval <- het[3]
      res_run$Q_pct <- het[4]
      res_run$Isq <- het[5]
      
      results_rep <- rbind(results_rep, res_run)  
      
    }
    # results_rep$mod <- mod
    
    results_models <- rbind(results_models, results_rep)  
    
  }
  results_models$setup_mode <- LD_mod
  
  results_all <- rbind(results_all, results_models)
}

save(results_all, file=sprintf(paste0(output_path, "/results_%s.rda"), job_id))


# write.csv(results_averaged, "./results//results_averaged.csv")



