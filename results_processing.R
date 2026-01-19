setwd("C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/")
source('modes_sims.R')
source('functions_sims.R')


results <- read.table("./results/ld_run5/results_ld_sims.csv")

reps= 10000  
  ## test params:
  # results <- results_all
  avg_res <- data.frame()
  
  mode_res <- data.frame()
  for (LD_mod in c(1,2,3,4)){
    rep_res <- data.frame()
    single_model_res <- results[results$setup_mode == LD_mod,]
    
    rep_res <- data.frame()
    setup_mode = 2
    for (model in c("A","B","C","D") ){
      params <- setup(setup_mode, model)
      b1 = params[4]
      b2 = params[5]
      
      
      row_res <- as.data.frame(matrix(NA, nrow = 1, ncol = 6))
      colnames(row_res) <- c("model","b","se","Q_pct","mean_Qsnps","cov_b")
      row_res$model <- model

      
      results_mode <-single_model_res[single_model_res$mode == model,]
     
      ## beta vectors for MC CIs
      b_ivw_1  <- results_mode$b

      
      ## calculate MC CIs
      row_res$mc_lci <- quantile(b_ivw_1,  probs = 0.025, na.rm = TRUE)

      
      
      row_res$mc_uci <-  quantile(b_ivw_1,  probs = 0.975, na.rm = TRUE)
      
      
         b1 <- rep(b1, time=nrow(results_mode))
      
        row_res$b     <- mean(results_mode$b, na.rm=T)
        row_res$se    <- mean(results_mode$se, na.rm=T)
        row_res$sig   <- sum(ifelse(results_mode$pval < 0.05,1,0), na.rm=T)
        row_res$bias   <- mean((results_mode$b - b1), na.rm=T)
        row_res$F_stat <- mean(results_mode$F_stat, na.rm=T)
        row_res$mse <- mean((results_mode$b - b1)^2, na.rm=T)
        row_res$nsnp <- mean(results_mode$nsnp, na.rm = T)
        row_res$Q_pct <- (sum(as.numeric(results_mode$Qpval) < 0.05) / reps ) * 100
        row_res$mean_Qsnps <- mean(as.numeric(results_mode$Q_pct), na.rm=T)  
        row_res$mean_Isq <- mean(as.numeric(results_mode$Isq))
      ## calculate coverage
      
      results_mode <- cbind.data.frame(results_mode, (results_mode$b - (1.96 * results_mode$se)),((results_mode$b + (1.96 * results_mode$se)))) 
      names(results_mode)[13:14] <- c("lci","uci")
      
      row_res$lci <- mean(results_mode$lci, na.rm=T)
      row_res$uci <- mean(results_mode$uci, na.rm=T)
      
      
      row_res$cov_b <- sum(between(b1, results_mode$lci, results_mode$uci), na.rm=T)

               

      rep_res <- rbind(rep_res, row_res)
      
    }
    rep_res$LD_mod <- LD_mod
    mode_res <-rbind(mode_res,rep_res)
    
  }
    avg_res <-rbind(avg_res,mode_res)
  
  
  # avg_res <- avg_res %>%
  #   mutate(across(c(b, se, cov_b, sig, bias, F_stat, nsnp, mse, lci, uci), as.numeric))  # Replace with actual column names
  # 
  avg_res$cov_p <- (avg_res$cov_b/ reps)*100
  
  # avg_res <- unique(avg_res[avg_res$exposure ==1,])
  
  write.table(avg_res, "C:/Users/kb22541/Desktop/Analyses/simulation/mvmr_simulations/results/LD_sims_MC_CIS.csv")
 