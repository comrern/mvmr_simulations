setwd("C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/")

results <- read.table("./results_full_sims.csv")
results <- results[!is.na(results$se),]

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
      
      
      row_res <- as.data.frame(matrix(NA, nrow = 3, ncol = 10))
      colnames(row_res) <- c("model","method","exposure","b","se","cov_b","F_stat", "nsnp", "mc_lci","mc_uci")
      row_res$model <- c(model, model, model)
      row_res$method <- c("IVW","MVMR","MVMR")
      row_res$exposure <- c(1,1,2)
      
      results_mode <-single_model_res[single_model_res$mode == model,]
      
      ## beta vectors for MC CIs
      b_ivw_1  <- results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1, ]$b
      b_mvmr_1 <- results_mode[results_mode$method == "mvmr" & results_mode$exp == 1, ]$b
      b_mvmr_2 <- results_mode[results_mode$method == "mvmr" & results_mode$exp == 2, ]$b
      
      ## calculate MC CIs
      row_res$mc_lci <- c(
        quantile(b_ivw_1,  probs = 0.025, na.rm = TRUE),
        quantile(b_mvmr_1, probs = 0.025, na.rm = TRUE),
        quantile(b_mvmr_2, probs = 0.025, na.rm = TRUE)
      )
      
      row_res$mc_uci <- c(
        quantile(b_ivw_1,  probs = 0.975, na.rm = TRUE),
        quantile(b_mvmr_1, probs = 0.975, na.rm = TRUE),
        quantile(b_mvmr_2, probs = 0.975, na.rm = TRUE)
      )
      
      
      row_res$b     <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1,]$b, na.rm = T)
                            , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$b, na.rm = T)
                            , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$b, na.rm = T))
      
      row_res$se    <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$se, na.rm = T)
                            , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$se, na.rm = T)
                            , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$se, na.rm = T))
      
      row_res$sig   <- list(sum(ifelse(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$p < 0.05,1,0), na.rm = T)
                            , sum(ifelse(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$p < 0.05,1,0), na.rm = T)
                            , sum(ifelse(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$p < 0.05,1,0), na.rm = T))
      
      
      
      ## calculate confidence intervals
      
      results_mode <- cbind.data.frame(results_mode, (results_mode$b - (1.96 * results_mode$se)),((results_mode$b + (1.96 * results_mode$se))))
      names(results_mode)[10:11] <- c("bs_lci","bs_uci")
      
      ## coverage
      
      b1_v <- rep(b1, time=nrow(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]))
      b2_v <- rep(b2, time=nrow(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]))
      b1_vmvmr <- rep(b1, time=nrow(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]))
      
      row_res$cov_b <- list(sum(between(b1_v, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$bs_lci, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$bs_uci), na.rm = T)
                            , sum(between(b1_vmvmr, results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$bs_lci, results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$bs_uci), na.rm = T)
                            ,  sum(between(b2_v, results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$bs_lci, results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$bs_uci), na.rm = T))
      
      
      
      
      row_res$bias   <- list(mean((results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$b - b1_v), na.rm = T)
                             , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$b - b1_vmvmr), na.rm = T)
                             , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$b - b2_v), na.rm = T))
      
      row_res$F_stat <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$F_stat, na.rm = T)
                             , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$F_stat), na.rm = T)
                             , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$F_stat), na.rm = T))
      
      row_res$nsnp <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$nsnp, na.rm = T)
                           , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$nsnp), na.rm = T)
                           , mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$nsnp), na.rm = T))
      
      row_res$mse <- list(mean((results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1,]$b - b1_v)^2 , na.rm = T),
                          mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 1,]$b - b1_vmvmr)^2, na.rm = T),
                          mean((results_mode[results_mode$method == "mvmr" & results_mode$exp == 2,]$b - b2_v)^2, na.rm = T))
      
      
      rep_res <- rbind(rep_res, row_res)
    }
    rep_res$setup_mode <- setup_mode
    avg_res <-rbind(avg_res,rep_res)
  }
  avg_res <-rbind(avg_res,rep_res)
  
  
  avg_res <- avg_res %>%
    mutate(across(c(b, se, cov_b, sig, bias, F_stat, nsnp, mse, mc_lci, mc_uci), as.numeric))
  
  avg_res$cov_p <- (avg_res$cov_b/ reps)*100
  
  View(avg_res[avg_res$exposure ==1,])
  
write.table(avg_res, "C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/mainsims_MCCIs.csv")
 