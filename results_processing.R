# setwd("C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/")

results <- read.table("./results_ld_rerun.csv")

reps= 2  
  ## test params:
  # results <- results_all
  avg_res <- data.frame()
  
  mode_res <- data.frame()
  for (LD_mod in c(TRUE,FALSE)){
    rep_res <- data.frame()
    single_model_res <- results[results$setup_mode == LD_mod,]
    
    rep_res <- data.frame()
    setup_mode = 2
    for (model in c("A","B","C","D") ){
      params <- setup(setup_mode, model)
      b1 = params[4]
      b2 = params[5]
      
      
      row_res <- as.data.frame(matrix(NA, nrow = 3, ncol = 6))
      colnames(row_res) <- c("model","method","exposure","b","se","cov_b")
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
                             , mean(abs(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$F_stat))
                             , mean(abs(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$F_stat)))
      
      rep_res <- rbind(rep_res, row_res)
      
    }
    rep_res$LD_mod <- LD_mod
    mode_res <-rbind(mode_res,rep_res)
    
  }
    avg_res <-rbind(avg_res,mode_res)
  
  
  avg_res <- avg_res %>%
    mutate(across(c(b, se, cov_b, sig, bias), as.numeric))  # Replace with actual column names
  
  avg_res$cov_p <- (avg_res$cov_b/ reps)*100
  
  View(avg_res[avg_res$exposure ==1,])
  
  # write.table(avg_res, "C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/avergaed_results_fullsims.csv")
 