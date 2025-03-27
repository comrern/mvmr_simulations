
library(dplyr)
library(tidyverse)

results <- read.csv("./combined_results10k.csv")

# results <- results_models

source("./modes_sims.R")
source("./functions_sims.R")

reps=10
rep_res <- data.frame()
avg_res <- data.frame()

  
  for (model in c("A","B","C","D") ){
    
    params <- setup("1", model)
    b1 = params[4]
    b2 = params[5]
    
    
    row_res <- as.data.frame(matrix(NA, nrow = 3, ncol = 6))
    colnames(row_res) <- c("model","method","exposure","b","se","cov_b")
    row_res$model <- c(model)
    row_res$method <- c("IVW","MVMR","MVMR")
    row_res$exposure <- c(1,1,2)
    
    results_mode <-results[results$mode == model,]
    
    row_res$b     <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1,]$b)
                          , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$b)
                          , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$b))
    
    row_res$se    <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$se)
                          , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$se)
                          , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$se))
    
    row_res$sig   <- list(sum(ifelse(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$p < 0.05,1,0))
                           , sum(ifelse(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$p < 0.05,1,0))
                           , sum(ifelse(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$p < 0.05,1,0)))
    

    ## calculate coverage
    
    results_mode <- cbind.data.frame(results_mode, (results_mode$b - (1.96 * results_mode$se)),((results_mode$b + (1.96 * results_mode$se)))) 
    names(results_mode)[9:10] <- c("lci","uci")
    
    b1_v <- rep(b1, time=reps)
    b2_v <- rep(b2, time=reps)
    
    row_res$cov_b <- list(sum(between(b1_v, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$lci, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$uci))
                          , sum(between(b1_v, results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$lci, results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$uci))
                          ,  sum(between(b2_v, results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$lci, results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$uci)))
    
    row_res$bias   <- list(abs(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$b - b1_v))
                                , abs(mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$b - b1_v))
                                , abs(mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$b - b2_v)))
    
    row_res$MSE   <- list(mean((results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$b - b1_v)^2)
                          , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$b - b1_vmvmr)^2)
                          , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$b - b2_v)^2))
    
    
    rep_res <- rbind(rep_res, row_res)
  }
  avg_res <-rbind(avg_res,rep_res)


  
  avg_res <- avg_res %>%
    mutate(across(c(b, se, cov_b, sig, bias, MSE), as.numeric))  # Replace with actual column names
  
avg_res$cov_p <- (avg_res$cov_b/ reps)*100


  
  # write.table(avg_res, "/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/real_data_10k.csv")
  

