set.seed(1234)

setwd("C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/")


library(dplyr)
library(tidyverse)
library(MASS)
library(TwoSampleMR)
library(MVMR)
library(truncnorm)
library(ggplot2)


source('modes_sims.R')
source('functions_sims.R')

reps = 2
run = 0
results = data.frame()
results_all = NULL
results_ivw = NULL
results_out = NULL
mvmrres <- NULL

results_all = data.frame()
results_ivw = data.frame()

for (setup_mode in c(1,2,3,4)){
  
  for (model in c("A","B","C","D"))  {
  results_rep = data.frame()
  run=0
      for(j in 1:reps){  
          
            run <- run + 1  
            print(paste("on run", run, model))
            
            params <- setup(setup_mode, model)
            snps = params[1]      
            snpsc = params[2]         
            nobs = params[3]
            b1 = params[4]
            b2 = params[5]
            betaC=params[6]
            beta2C=params[7]
            pi = 0.5
            xi=0
            
            
            dat <- data_gen(snps, snpsc, nobs, b1 , b2, betaC, beta2C, pi)
              #(no of snps, snps for confounding var, samplesize, beta1, beta2, snp-confounder effect)
          #### 
          
            results[1,"model"] <- model
            results[1,"pi"] <- pi
            results[1,"sample.size"] <- nobs
            
          ####Regression check######
          
            ols <- summary(lm(Y ~ X1, data = dat))
            
            results[1,"ols_b"] <- ols$coefficients["X1","Estimate"]
            results[1,"ols_se"] <- ols$coefficients["X1","Std. Error"]
            
          ## allele freq
          
            MR_dat <- GWASres(dat)
          
          #### Univariate two sample MR ####
            
            univariate_results <- univariate_MR(MR_dat)
            
          ### MVMR
            
            mvmr_res <- run_mvmr(MR_dat)
            
        ## format results
          res_run <- rbind(univariate_results, mvmr_res)
          res_run$run <- j
          res_run$mode <- model
          
          results_rep <- rbind(results_rep, res_run)  
          
      }
    # results_rep$mod <- mod
  
    results_all <- rbind(results_all, results_rep)  
  }
}
results_ivw <- results_all[results_all$method %in% c("Inverse variance weighted","mvmr"),]
    
results_averaged <- avg_cals(results_ivw,reps)
results_averaged$b <- as.numeric(results_averaged$b)
results_averaged$se <- as.numeric(results_averaged$se)
results_averaged$nsnp <- as.numeric(results_averaged$nsnp)
results_averaged$p <- as.numeric(results_averaged$p)
results_averaged$cov_b <- as.numeric(results_averaged$cov_b)

 # write.csv(results_averaged, "./results/results_averaged.csv")
 
 
 
 
