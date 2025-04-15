

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
  
  
  if (m == 1) {
    ld_mag = 0
    
  } 
  if (m == 2) {
    ld_mag = 0.25
      } 
  if (m == 3) {
    ld_mag = 0.5
  } 
  if (m == 4) {
    ld_mag = 0.75
      } 
  
  snpsc =  33
  LD_mod = T
  snps= 33
  nobs = 25000
  betaC = 0.8
  beta2C = 0.8
  xi = 0
  
  
  return(c(snps, snpsc, nobs, b1 , b2, betaC, beta2C, xi, LD_mod, ld_mag))
}