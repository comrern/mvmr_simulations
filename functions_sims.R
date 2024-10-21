


data_gen <- function(nsnps,snpsc,ss,beta1,beta2, betaC, beta2c, pi){
  
  n=2*ss
  
  ## debug code ##      n=20000     nsnps=99     snpsc=99      beta1=0    beta2=0.4   betaC=0.5  beta2C=0.6  pi=0.5
  
  df <- as.data.frame(matrix(nrow=n))
  df$V1 <- seq.int(nrow(df))
  df$X2 <- rtruncnorm(n, a=0.0001, b=0.9999, mean= 0.276, sd= 0.1443219)               ## based on observed data

  G  <- matrix(rbinom(n*snpsc, 2, 0.4), n, nsnps)
  G2 <- matrix(rbinom(n*snpsc, 2, 0.4), n, snpsc)
  
  means <- c(0, 0)                                   
  cov_matrix <- matrix(c(1, 0, 0, 1),
                       ncol = 2)
  
  # create bivariate normal distribution
  errors <- mvrnorm(n = n,
                    mu = means, 
                    Sigma = cov_matrix)
  
  v_x1 <- errors[,1]
  v_y <- errors[,2]
  v_c <- rnorm(n,0,1)
  
  effs_x1 <- abs(rnorm(nsnps,0,0.08))

  
  df <- (cbind(df, G, G2))
  # colnames(df) <- gsub("V","G",colnames(df))

  LD <- 0.5 + (df[,"X2"])
  LD <- ifelse(LD> 1,1 , LD)
  df[,"C"] <-  beta2C*df[,"X2"] + v_c 
  df[,"X1"] <- (G[,]%*%effs_x1)* LD + betaC*df[,"C"] + v_x1
  df[,"Y"] <- beta1*df[,"X1"] + beta2*df[,"X2"] + betaC*df[,"C"] + v_y  
  
  
  data <- df
  return(df)
}


GWASres <- function(dat){
  MR_dat = data.frame()
  
  dat.1 <- dat[1:nobs,]
  dat.2 <- dat[(nobs+1):(2*nobs),]
  
  est.snps <- snps + snpsc
  
  for(i in 1:est.snps){
    a <- summary(lm(dat.1$X1~dat.1[,i]))
    MR_dat[i,"X1_b"] <- a$coefficient[2,1]
    MR_dat[i,"X1_se"] <- a$coefficient[2,2]
    MR_dat[i,"X1_p"] <- a$coefficient[2,4]
    MR_dat[i,"X1_r2"] <- a$r.squared
    b <- summary(lm(dat.1$X2~dat.1[,i]))
    MR_dat[i,"X2_b"] <- b$coefficient[2,1]
    MR_dat[i,"X2_se"] <- b$coefficient[2,2]
    MR_dat[i,"X2_p"] <- b$coefficient[2,4]
    MR_dat[i,"X2_r2"] <- b$r.squared
    c<-summary(lm(dat.2$Y~dat.2[,i]))
    MR_dat[i,"Y_b"] <- c$coefficient[2,1]
    MR_dat[i,"Y_se"] <-c$coefficient[2,2]
    MR_dat[i,"Y_p"] <- c$coefficient[2,4]
    

    
  }
  
  allele_frequencies <- colSums(dat[,3:(snps + snpsc + 2)]) / (2 * nrow(dat))
  MR_dat$af <- allele_frequencies
  MR_dat$id <- seq.int(nrow(MR_dat))
  MR_dat$EA <- "A"
  
  return(MR_dat)
}

univariate_MR <- function(MR_dat){
  
  exp1 <- format_data(MR_dat, type="exposure",
                      snp_col="id",
                      beta_col="X1_b",
                      se_col="X1_se",
                      eaf_col="af",
                      pval_col="X1_p",
                      effect_allele_col = "EA"
  )
  
  exp2 <- format_data(MR_dat, type="exposure",
                      snp_col="id",
                      beta_col="X2_b",
                      se_col="X2_se",
                      eaf_col="af",
                      pval_col="X1_p",
                      effect_allele_col = "EA"
  )
  
  out <- format_data(MR_dat, type="outcome",
                     snp_col="id",
                     beta_col="Y_b",
                     se_col="Y_se",
                     eaf_col="af",
                     pval_col="X1_p",
                     effect_allele_col = "EA"
  )
  
  dat1 <- harmonise_data(exp1, out)
  dat2 <- harmonise_data(exp2, out)
  
  dat1 <- dat1[dat1$pval.exposure <= 5e-8,]
  dat2 <- dat2[dat2$pval.exposure <= 5e-8,]
  
  if (length(dat1$SNP) == 0 | length(dat2$SNP) == 0) {print("No significant SNPs for exposure(s)")}
  if (length(dat1$SNP) < 3 | length(dat2$SNP) < 3) {print("Too few SNPs for MR")}
  
  mr1 <- mr(dat1)
  mr2 <- mr(dat2)
  
  mr1$exp <- 1
  mr2$exp <- 2
  
  univ_results <- rbind(mr1, mr2)
  univ_results <- univ_results[,5:10]
  
  return(univ_results)
}

run_mvmr <- function(MR_dat){
  
  
  dat_formatted <- format_mvmr(
    BXGs = MR_dat[,c("X1_b","X2_b")], 
    BYG = MR_dat$Y_b, 
    seBXGs = MR_dat[,c("X1_se","X2_se")], 
    seBYG=MR_dat$Y_se, 
    RSID=MR_dat$id
  )
  
  dat_IS <- strength_mvmr(r_input= dat_formatted, gencov = 0)
  
  res_mvmr <- as.data.frame(ivw_mvmr(dat_formatted))
  res_mvmr$method <- "mvmr"
  res_mvmr$exposure <- c(1,2)
  res_mvmr$nsnp <- c(50, 50)                                ## fix properly ##
  res_mvmr <- res_mvmr[,c(5,7,1,2,4,6)]
  colnames(res_mvmr) <- c("method","nsnp","b","se","pval","exp")
    
  pres <- pleiotropy_mvmr(r_input = dat_formatted, gencov = 0)
  
  
  return(res_mvmr)
}


avg_cals <- function(results, reps) {
  
  
  
  ## test params:
  results <- results_ivw
  avg_res <- data.frame()

  
  
  for (model in c("A","B","C","D") ){
    
    params <- setup(model)
    b1 = params[4]
    b2 = params[5]
    
    
    row_res <- as.data.frame(matrix(NA, nrow = 4, ncol = 8))
    colnames(row_res) <- c("model","method","exposure","b","se","p","nsnp","cov_b")
    row_res$model <- c(model, model)
    row_res$method <- c("IVW","IVW","MVMR","MVMR")
    row_res$exposure <- c(1,2,1,2)
    
    results_mode <-results[results$mode == model,]
    
    row_res$b     <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1,]$b)
                  , mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 2 ,]$b)
                  , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$b)
                  , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$b))
    
    row_res$se    <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$se)
                  , mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 2 ,]$se)
                  , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$se)
                  , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$se))
    
    row_res$p     <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$pval)
                  , mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 2 ,]$pval)
                  , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$pval)
                  , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$pval))
    
    row_res$nsnp   <- list(mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$nsnp)
                  , mean(results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 2 ,]$nsnp)
                  , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$nsnp)
                  , mean(results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$nsnp))
    
    

    
    
  ## calculate coverage
    
    results_mode <- cbind.data.frame(results_mode, (results_mode$b - (1.96 * results_mode$se)),((results_mode$b + (1.96 * results_mode$se)))) 
    names(results_mode)[9:10] <- c("lci","uci")
    
    b1_v <- rep(b1, time=reps)
    b2_v <- rep(b2, time=reps)
    
    row_res$cov_b <- list(sum(between(b1_v, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$lci, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 1 ,]$uci))
                    , sum(between(b2_v, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 2 ,]$lci, results_mode[results_mode$method == "Inverse variance weighted" & results_mode$exp == 2 ,]$uci))
                    , sum(between(b1_v, results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$lci, results_mode[results_mode$method == "mvmr" & results_mode$exp == 1 ,]$uci))
                    ,  sum(between(b2_v, results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$lci, results_mode[results_mode$method == "mvmr" & results_mode$exp == 2 ,]$uci)))
    
    
    
    avg_res <- rbind(avg_res, row_res)
    }
   
  
 return(avg_res) 
}






