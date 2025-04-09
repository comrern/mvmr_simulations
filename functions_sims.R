


data_gen <- function(nsnps,snpsc,ss,beta1,beta2, betaC, beta2c){
  
  n=2*ss
  
  ## debug code ##      n=25000     nsnps=33     snpsc=33      beta1=0    beta2=0.4   betaC=0.5  beta2C=0.6  pi=0.5
  
  df <- as.data.frame(matrix(nrow=n))
  df$V1 <- seq.int(nrow(df))
  df$X2 <- rtruncnorm(n, a=0.0001, b=0.9999, mean= 0.276, sd= 0.1443219)               ## based on observed data
  

    prob_inc <-  0.3 + 0.1 * df$X2  ## build probability vector based on value of X2 --> 
    ## each observation of G binom distribution has probability dependent on value of X2 meaning higher X2 = higher AF
    
    prob_dec <-  0.4 - 0.1 * df$X2
    
    prob_inc_g <- rep(prob_inc, times = nsnps/3)
    prob_dec_g <- rep(prob_dec, times = nsnps/3)
    
    G_inc <-  matrix(rbinom(n*(nsnps/3), 2, prob_inc), n, (nsnps/3))
    G_dec <-  matrix(rbinom(n*(nsnps/3), 2, prob_dec), n, (nsnps/3))
    G_cont <-  matrix(rbinom(n*(nsnps/3), 2, 0.4), n, (nsnps/3))
    
    G <- cbind(G_inc, G_dec, G_cont)
    
  

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
  df[,"C"] <-  beta2C*df[,"X2"] + v_c 
  
  # model LD
  LD_inc <- 0.5 + 0.5 *(df[,"X2"])
  LD_inc <- ifelse(LD_inc> 1,1 , LD_inc)
  
  LD_dec <- 1 - 0.5 * (df[,"X2"])
  LD_dec <- ifelse(LD_dec> 1,1 , LD_dec)
  
  LD_inc_mat  <- sapply(effs_x1[1:(nsnps/3)], function(y_val) LD_inc * y_val)
  LD_dec_mat  <- sapply(effs_x1[((nsnps/3)+1):(2*(nsnps/3))], function(y_val) LD_dec * y_val)
  
  LD_const_mat <- matrix(rep(effs_x1[(2*(nsnps/3)+1):nsnps], each = n), nrow = n, ncol = nsnps/3, byrow = TRUE)
  
  effs_mat <- cbind(LD_inc_mat, LD_dec_mat, LD_const_mat)
  
  df[,"X1"] <- rowSums(G[,]*effs_mat) + xi*df["X2"] + betaC*df[,"C"] + v_x1
  
  df[,"X1_novar"] <- G[,]%*%effs_x1 + xi*df["X2"] + betaC*df[,"C"] + v_x1  
  
  df[,"Y"] <- beta1*df[,"X1"] + beta2*df[,"X2"] + betaC*df[,"C"] + v_y  
  
  
  data <- df
  return(df)
}



GWASres <- function(dat, LD_mod){
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
    d<-summary(lm(dat.1$X1_novar~dat.1[,i]))
    MR_dat[i,"x_novar_b"] <- d$coefficient[2,1]
    MR_dat[i,"x_novar_se"] <-d$coefficient[2,2]
    MR_dat[i,"x_novar_p"] <- d$coefficient[2,4]
    
    if (LD_mod==T) {
      MR_dat[i,"X1_b"] <- d$coefficient[2,1]
      MR_dat[i,"X1_se"] <- d$coefficient[2,2]
      MR_dat[i,"X1_p"] <- d$coefficient[2,4]
      MR_dat[i,"X1_r2"] <- d$r.squared
      
      
       }
    
    
    }
  
  allele_frequencies <- colSums(dat[,3:(snps + snpsc + 2)]) / (2 * nrow(dat))
  MR_dat$af <- allele_frequencies
  MR_dat$id <- seq.int(nrow(MR_dat))
  MR_dat$EA <- "A"

  return(MR_dat)
}

univariate_MR <- function(MR_dat, LD_mod){
  
  
  

  
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
  
  if (LD_mod==FALSE) {
    
    dat1 <- dat1[dat1$SNP %in% MR_dat[MR_dat$x_novar_p < 5e-8,]$id,]
    
  } else   dat1 <- dat1[dat1$pval.exposure <= 5e-8,]

  dat2 <- dat2[dat2$pval.exposure <= 5e-8,]
  
  f_1 <- mean((dat1$beta.exposure^2)/ (dat1$se.exposure^2))
  f_2 <- mean((dat2$beta.exposure^2)/ (dat2$se.exposure^2))
  
  mr1 <- mr(dat1)
  
  mr2 <- mr1
  
  mr1$exp <- 1
  mr2$exp <- 2
  
  univ_results <- rbind(mr1, mr2)
  univ_results <- univ_results[,5:10]
  univ_results$F_stat <- c(f_1,f_2)
  return(univ_results)
}

run_mvmr <- function(MR_dat, dat, LD_mod){
 
 if (LD_mod==FALSE) {
    
     MR_dat <- MR_dat[(MR_dat$id %in% MR_dat[MR_dat$x_novar_p < 5e-8,]$id) | (MR_dat$X2_p < 5e-8),]
    
  } else    MR_dat <- MR_dat[(MR_dat$X1_p < 5e-8) | (MR_dat$X2_p < 5e-8) ,]
  

  
  if (length(MR_dat) >= 1) {
    
    n_x1 <- sum(MR_dat$X1_p < 5e-8)
    n_x2 <- sum(MR_dat$X2_p < 5e-8)
    
    dat_formatted <- format_mvmr(
      BXGs = MR_dat[,c("X1_b","X2_b")], 
      BYG = MR_dat$Y_b, 
      seBXGs = MR_dat[,c("X1_se","X2_se")], 
      seBYG=MR_dat$Y_se, 
      RSID=MR_dat$id
    )
    
    
    
    if (setup_mode == 4){
      cov <- snpcov_mvmr(dat[,MR_dat$id], dat[,c("X1","X2")])
    } else    cov <- snpcov_mvmr(dat[,MR_dat$id], dat[,c("X1","X2")])
    
    
    
    strength <- strength_mvmr(dat_formatted, gencov = cov)
    res_mvmr <- as.data.frame(ivw_mvmr(dat_formatted, gencov = cov))
    res_mvmr$method <- "mvmr"
    res_mvmr$exposure <- c(1,2)
    res_mvmr$nsnp <- c(n_x1, n_x2)                                ## fix properly ##
    res_mvmr <- res_mvmr[,c(5,7,1,2,4,6)]
    colnames(res_mvmr) <- c("method","nsnp","b","se","pval","exp")
    res_mvmr$F_stat <- c(strength$exposure1, strength$exposure2)
  } else {
    
    res_mvmr <- as.data.frame(matrix(rep(NA, 2 * 3), nrow = 2, ncol = 3))
    res_mvmr$method <- "mvmr"
    res_mvmr$exposure <- c(1,2)
    res_mvmr$nsnp <- c(0, 0)   
    res_mvmr <- res_mvmr[,c(4,6,1,2,3,5)]
    colnames(res_mvmr) <- c("method","nsnp","b","se","pval","exp")
    res_mvmr$F_stat <- c(strength$exposure1, strength$exposure2)
    
    
  }
  
  
  return(res_mvmr)
}








