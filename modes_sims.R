

setup <- function(m,k){
  
  if(k=='A'){    ## No b1 effect, all obs due to b2

    b1 = 0.4
    b2 = 0
    LD_mag = 0

  }

  if(k=='B'){    ## b1 and b2 effect, b2 modifies magnitude of ce

    b1 = 0.4
    b2 = 0
    LD_mag = 0.125
    

  }
  
  if(k=='C'){    ## No b1 or b2 effect
    
    b1 = 0.4
    b2 = 0
    LD_mag = 0.25
    

  }
  
  if(k=='D'){    ## No ancestry effect
    
    b1 = 0.4
    b2 = 0
    LD_mag = 0.5
    
  }
  
  
  snpsc = 33
  snps=33
  nobs = 25000
  betaC <- 0
  beta2C <- 0


  
  return(c(snps, snpsc, nobs, b1 , b2, betaC, beta2C,  LD_mag))
}