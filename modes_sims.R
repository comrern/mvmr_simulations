

setup <- function(k){
  
  if(k=='A'){
    
    snps = 100       #no of SNPs for X1
    snpsc = 100         #No of SNPs for X2/X3
    nobs = 10000
  }

  
  return(c(snps, snpsc, nobs))
}