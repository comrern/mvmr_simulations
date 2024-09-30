

setup <- function(k){
  
  if(k=='A'){    ## No b1 effect, all obs due to b2
    
    snps = 99       #no of SNPs for X1
    snpsc = 99         #No of SNPs for X2/X3
    nobs = 25000
    b1 = 0
    b2 = 0.4
    betaC = 0.4
    beta2C =0.5
  }

  if(k=='B'){    ## b1 and b2 effect, b2 modifies magnitude of ce
    
    snps = 99       #no of SNPs for X1
    snpsc = 99         #No of SNPs for X2/X3
    nobs = 25000
    b1 = 0.4
    b2 = 0.4
    betaC = 0.4
    beta2C =0.5
  }
  
  if(k=='C'){    ## No b1 or b2 effect
    
    snps = 99       #no of SNPs for X1
    snpsc = 99         #No of SNPs for X2/X3
    nobs = 25000
    b1 = 0
    b2 = 0
    betaC = 0.4
    beta2C =0.5
  }
  
  if(k=='D'){    ## No ancestry effect
    
    snps = 99       #no of SNPs for X1
    snpsc = 99         #No of SNPs for X2/X3
    nobs = 20000
    b1 = 0.4
    b2 = 0
    betaC = 0.4
    beta2C =0.5
  }
  
  
  return(c(snps, snpsc, nobs, b1 , b2, betaC, beta2C, pi))
}