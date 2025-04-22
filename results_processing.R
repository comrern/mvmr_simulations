setwd("C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/")

results <- read.table("./results_full_sims.csv")


source('../modes_sims.R')
source('../functions_sims.R')
reps= 10000
  ## test params:
  # results <- results_all
  avg_res <- data.frame()
  model_res <- data.frame()
  for (setup_mode in c(1,2,3,4)){
    rep_res <- data.frame()
    single_model_res <- results[results$setup_mode == setup_mode,]
    
    for (model in c("A","B","C","D") ){
      
      params <- setup(setup_mode, model)
      b1 = params[4]
      b2 = params[5]
      
      
      row_res <- as.data.frame(matrix(NA, nrow = 3, ncol = 7))
      colnames(row_res) <- c("model","method","exposure","b","se","cov_b","F_stat")
      row_res$model <- c(model, model, model)
      row_res$method <- c("IVW","MVMR","MVMR")
      row_res$exposure <- c(1,1,2)
      
      results_mode <-single_model_res[single_model_res$mode == model,]
      
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
      names(results_mode)[10:11] <- c("lci","uci")
      
      
      b1_v <- rep(b1, time=nrow(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]))
      b2_v <- rep(b2, time=nrow(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]))
      b1_vmvmr <- rep(b1, time=nrow(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]))
      
      row_res$cov_b <- list(sum(between(b1_v, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$lci, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$uci))
                            , sum(between(b1_vmvmr, results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$lci, results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$uci))
                            ,  sum(between(b2_v, results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$lci, results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$uci)))

      
      row_res$bias   <- list(mean((results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$b - b1_v))
                             , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$b - b1_vmvmr))
                             , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$b - b2_v)))
      
      row_res$F_stat <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$F_stat)
                             , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$F_stat))
                             , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$F_stat)))
      
      row_res$nsnp <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$nsnp, na.rm = T)
                           , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$nsnp), na.rm = T)
                           , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$nsnp), na.rm = T))
      
      row_res$mse <- list(mean((results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1,]$b - b1_v)^2),
                          mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 1,]$b - b1_vmvmr)^2),
                          mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 2,]$b - b2_v)^2))
      
      
      rep_res <- rbind(rep_res, row_res)
    }
    rep_res$setup_mode <- setup_mode
    avg_res <-rbind(avg_res,rep_res)
  }
  avg_res <-rbind(avg_res,rep_res)
  
  
  avg_res <- avg_res %>%
    mutate(across(c(b, se, cov_b, sig, bias, F_stat, nsnp, mse), as.numeric))  # Replace with actual column names
  
  avg_res$cov_p <- (avg_res$cov_b/ reps)*100
  
  View(avg_res[avg_res$exposure ==1,])
  
  # write.table(avg_res, "C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/avergaed_results_fullsims.csv")
 