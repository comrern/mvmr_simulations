


data_gen <- function(nsnps,snpsc,ss,beta1,beta2, betaC, beta2c, LD_mod, LD_mag){
  
  n=2*ss
  
  ## debug code ##      n=50000     nsnps=33     snpsc=33      beta1=0.4    beta2=0.4   betaC=0.5  beta2C=0.6  pi=0.5
  
  df <- as.data.frame(matrix(nrow=n))
  df$V1 <- seq.int(nrow(df))
  df$X2 <- rtruncnorm(n, a=0.0001, b=0.9999, mean= 0.276, sd= 0.1443219)               ## based on observed data
  
  local_anc <- matrix(rbinom(n*nsnps, 2, df[,"X2"]), n, nsnps)
  
  G <- matrix(rbinom(n * nsnps, 2,
                     ifelse(local_anc == 2, 0.4,
                            ifelse(local_anc == 1, 0.35, 0.3))),
              n, nsnps)
    
  
  # G <- matrix(rbinom(n*nsnps, 2, 0.4), n, nsnps)
  
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
  
  
  if(LD_mod ==4){
    snp_dir <- sample(c(-1, 1), size = length(effs_x1), replace = TRUE)
    
    dir_mat <- matrix(snp_dir, nrow = nrow(local_anc), ncol = length(effs_x1), byrow = TRUE)
    
    effs_mat <- matrix(effs_x1, nrow = nrow(local_anc), ncol = length(effs_x1), byrow = TRUE) *
      (1 - ((local_anc / 2) * LD_mag * dir_mat))
    
    
  } else{
    effs_mat <- matrix(effs_x1, nrow = nrow(local_anc), ncol = length(effs_x1), byrow = TRUE) *
      (1 - ((local_anc / 2) * LD_mag))
  }
  df[,"X1"] <- rowSums(G[,]*effs_mat) +  betaC*df[,"C"] + v_x1
  
  
  df[,"X1_novar"] <- G[,]%*%effs_x1 + betaC*df[,"C"] + v_x1  
  
  df[,"Y"] <- beta1*df[,"X1"] + betaC*df[,"C"] + v_y  
  
  
  data <- df
  return(df)
}



GWASres <- function(dat, LD_mod){
  MR_dat = data.frame()
  
  ### Determine sample splitting based on LD_mod
  
  if (LD_mod ==3){
    dat <- dat[order(dat$X2), ]
  } 
  
  dat.1 <- dat[1:nobs,]
  dat.2 <- dat[(nobs+1):(2*nobs),]
  
  est.snps <- snps + snpsc
  
  for(i in 1:est.snps){
    a <- summary(lm(dat.1$X1~dat.1[,i]))
    MR_dat[i,"X1_b"] <- a$coefficient[2,1]
    MR_dat[i,"X1_se"] <- a$coefficient[2,2]
    MR_dat[i,"X1_p"] <- a$coefficient[2,4]
    MR_dat[i,"X1_r2"] <- a$r.squared
    b <- summary(lm(dat.2$X1~dat.2[,i]))
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

    if (LD_mod==1) {
      MR_dat[i,"X1_b"] <- d$coefficient[2,1]
      MR_dat[i,"X1_se"] <- d$coefficient[2,2]
      MR_dat[i,"X1_p"] <- d$coefficient[2,4]


       }
    
    
    }
  
  allele_frequencies <- colSums(dat[,3:(snps + snpsc + 2)]) / (2 * nrow(dat))
  MR_dat$af <- allele_frequencies
  MR_dat$id <- seq.int(nrow(MR_dat))
  MR_dat$EA <- "A"
  MR_dat <- MR_dat[MR_dat$X1_p <= 5e-8,]
  return(MR_dat)
}


heterogeneity <- function(MR_dat){
  
  Q_df <-  data.frame(ID = numeric(),
                           Qsnp = numeric(),
                           Qp = numeric(),
                            X1_b = numeric(),
                          X2_b = numeric(),
                           stringsAsFactors = FALSE)
  
  for (i in 1:nrow(MR_dat)){
    
    betas <- MR_dat[i,c("X1_b","X2_b")] 
    ses <- MR_dat[i,c("X1_se","X2_se")]
    
    w <- 1 / (ses)^2 # get weights
    ivw_b <- sum(betas * w) / sum(w) # ivw betas
    se <- sqrt(1 / sum(w)) # ivw se
    
    Q <- sum(w * (betas - ivw_b)^2)
    
   
    df <- length(betas) -1
  
    Qpval <- stats::pchisq(Q, df, lower.tail=FALSE)
    
    Q_df[i, ] <- list(i, Q, Qpval, betas[1, 1], betas[1, 2])
    
    
  
  }
  
  Qdf <- (length(Q_df$ID)  - 1)
  Q_all <- sum(Q_df$Qsnp)
  
  
  Q_total <-  c("Qsum",
                sum(Q_df$Qsnp), 
                pchisq(Q_all, df =  Qdf, lower.tail = FALSE),
                (mean(Q_df$Qp < 0.05) * 100),
                max(0, ((Q_all - Qdf)/ Q_all) * 100)
                )
  
  return(Q_total)  

  
}


univariate_MR <- function(MR_dat, LD_mod){
  

  f_1 <- mean((MR_dat$X1_b^2)/ (MR_dat$X1_se^2))

  
  mr1 <- as.data.frame(mr_ivw(MR_dat$X1_b,
                MR_dat$Y_b,
                MR_dat$X1_se,
                MR_dat$Y_se))
  
  # mr2 <- mr1
  univ_results <- mr1[,1:4]
  univ_results$exp <- 1
  # mr2$exp <- 2
  
  # univ_results <- rbind(mr1, mr2)

  univ_results$F_stat <- f_1
  return(univ_results)
}

run_mvmr <- function(MR_dat, dat, LD_mod){
 
 if (LD_mod==2) {
    
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








