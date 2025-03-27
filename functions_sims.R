


data_gen <- function(nsnps,snpsc,ss,beta1,beta2, betaC, beta2c, pi, LD_mod){
  
  n=2*ss
  
  ## debug code ##      n=25000     nsnps=28    snpsc=28   xi=0   beta1=0    beta2=0.4   betaC=0.5  beta2C=0.6  pi=0.5
  
  df <- as.data.frame(matrix(nrow=n))
  df$V1 <- seq.int(nrow(df))
  df$X2 <- rtruncnorm(n, a=0.0001, b=0.9999, mean= 0.276, sd= 0.1443219)               ## based on observed data
  df <- df[order(df$X2), ]
  
  obs_af <- read.table("../ld_reports/af_subset", header=T)
  
  # rep allele freqs to feed into binom 
  G <- data.frame(matrix(NA, nrow = nsnps, ncol = n))
  
  for (snp in 1:nrow(obs_af)){
    
    snp_af <- unlist(rep((obs_af[snp,]), each=(n/4)))
    G[snp, ] <- rbinom(n, 2, (snp_af))
    
  }
  
  G <- as.data.frame(t(G))
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
  
  effs_x1 <- abs(rnorm(28,0,0.5))
  
  
  df <- (cbind(df, G, G2))
  df[,"C"] <-  beta2C*df[,"X2"] + v_c 
  
  ### Model LD

  ## read observed LD variance 
  obs_ld_dat <- read.table("../ld_reports/test_ld_mat.txt", header=T)
  
  effs_ld_1 <- effs_x1 * obs_ld_dat$R2_p1
  effs_ld_2 <- effs_x1 * obs_ld_dat$R2_p2
  effs_ld_3 <- effs_x1 * obs_ld_dat$R2_p3
  effs_ld_4 <- effs_x1 * obs_ld_dat$R2_p4
  
  ld1_mat <- matrix(rep(effs_ld_1, each = n/4), nrow = length(effs_ld_1), ncol = n/4)
  ld2_mat <- matrix(rep(effs_ld_2, each = n/4), nrow = length(effs_ld_2), ncol = n/4)
  ld3_mat <- matrix(rep(effs_ld_3, each = n/4), nrow = length(effs_ld_3), ncol = n/4)
  ld4_mat <- matrix(rep(effs_ld_4, each = n/4), nrow = length(effs_ld_4), ncol = n/4)
  
  effs_mat <- t(cbind(ld1_mat, ld2_mat, ld3_mat, ld4_mat))
  
  df[,"X1"] <- rowSums(G[,]*effs_mat) + xi*df["X2"] + betaC*df[,"C"] + v_x1
  
  

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
  
  dat1_temp <- dat1[dat1$pval.exposure <= 5e-8,]
  dat2 <- dat2[dat2$pval.exposure <= 5e-8,]
  
  if (length(dat1_temp$SNP) == 0 | length(dat2$SNP) == 0) {
    print("No significant SNPs for exposure(s), reducing significance threshold")
    dat1_temp  <- dat1[dat1$pval.exposure <= 5e-7,]
    }
  if ((length(dat1_temp$SNP) < 3 | length(dat2$SNP) < 3) & (length(dat1_temp$SNP) == 0 | length(dat2$SNP))) {
    print("Too few SNPs for MR, reducing significance threshold")
    dat1_temp  <- dat1[dat1$pval.exposure <= 5e-7,]
    
    }
  
  mr1 <- mr(dat1_temp)

  mr1$exp <- 1

  univ_results <- mr1[,5:10]
  
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






