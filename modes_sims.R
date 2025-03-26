

setup <- function(m,k){
  
  if(k=='A'){    ## No b1 effect, all obs due to b2

    b1 = 0
    b2 = 0.8

  }

  if(k=='B'){    ## b1 and b2 effect, b2 modifies magnitude of ce

    b1 = 0.4
    b2 = 0.8

  }
  
  if(k=='C'){    ## No b1 or b2 effect
    
    b1 = 0
    b2 = 0

  }
  
  if(k=='D'){    ## No ancestry effect
    
    b1 = 0.4
    b2 = 0

  }
  
  
  snpsc = ifelse(m==4, 5, 33)
  LD_mod = ifelse(m==2, T, F)
  xi = ifelse(m==3, 1, 0)
  snps=33
  nobs = 25000
  betaC = 0.8
  beta2C = 0.8
  
  return(c(snps, snpsc, nobs, b1 , b2, betaC, beta2C, xi, LD_mod))
}