

setup <- function(k){
  
  if(k=='A'){    ## No b1 effect, all obs due to b2
    
    snps = 100       #no of SNPs for X1
    snpsc = 100         #No of SNPs for X2/X3
    nobs = 10000
    b1 = 0
    b2 = 0.8
  }

  if(k=='B'){    ## b1 and b2 effect, b2 modifies magnitude of ce
    
    snps = 100       #no of SNPs for X1
    snpsc = 100         #No of SNPs for X2/X3
    nobs = 10000
    b1 = 0.4
    b2 = 0.8
  }
  
  if(k=='C'){    ## No b1 or b2 effect
    
    snps = 100       #no of SNPs for X1
    snpsc = 100         #No of SNPs for X2/X3
    nobs = 10000
    b1 = 0
    b2 = 0
  }
  
  if(k=='D'){    ## No ancestry effect
    
    snps = 100       #no of SNPs for X1
    snpsc = 100         #No of SNPs for X2/X3
    nobs = 10000
    b1 = 0.4
    b2 = 0
  }
  
  
  return(c(snps, snpsc, nobs, b1, b2))
}