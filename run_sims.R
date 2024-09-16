set.seed(1234)

setwd("C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/")


library(dplyr)
library(tidyverse)
library(MASS)
library(TwoSampleMR)
library(MVMR)
library(truncnorm)

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


for(j in 1:reps){  
    
      gm <- 0.5
      model <-"A"
      params <- setup(model)
      
      snps = params[1]      
      snpsc = params[2]         
      nobs = params[3]
      b1 = params[4]
      b2 = params[5]
      pi = gm
      
      
      dat <- gendat(snps, snpsc, nobs, b1 , b2, pi)
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
  
    results_all <- rbind(results_all, res_run)  
      
}

results_ivw <- results_all[results_all$method %in% c("IVW","mvmr"),]
    
   
    