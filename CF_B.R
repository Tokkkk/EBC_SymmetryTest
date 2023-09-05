################################################
## This script is an auxiliary function for cf statistic from Bahraoui et. al 2018, 
## and this code is tranformed from Matlab code form https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=834
##############################################
### Input-> weight: 'N' Normal, 'DE' Double-exponential
#           sigma: Smoothing parameter > 0
#           a: 1st component
#           b: 2nd component
#Output-> vec: Vector of B and its 1st/2nd order derivatives at (a,b)
##############################################
DiagSym_CopulaCf_B <- function(weight,sigma,a,b) {
  x <- sigma * (a[1] - b[1])
  y <- sigma * (a[2] - b[2])
  z <- sigma * (a[2] - b[1])
  w <- sigma * (a[1] - b[2])
  
  if (weight == 'N') {
    TEMP1 <- exp(-((x^2 + y^2) / 2))
    TEMP2 <- exp(-((z^2 + w^2) / 2))
    
    v <- TEMP1 - TEMP2
    v1 <- sigma * (-x * TEMP1 + w * TEMP2)
    v2 <- sigma * (-y * TEMP1 + z * TEMP2)
    v11 <- sigma^2 * (-(x^2 - 1) * TEMP1 + z * w * TEMP2)
    v22 <- sigma^2 * (-(y^2 - 1) * TEMP1 + z * w * TEMP2)
    v12 <- sigma^2 * (-x * y * TEMP1 + (w^2 - 1) * TEMP2)
    v21 <- sigma^2 * (-x * y * TEMP1 + (z^2 - 1) * TEMP2)
  } else if (weight == 'DE') {
    TEMP1 <- (4 + x^2)^(-1)
    TEMP2 <- (4 + y^2)^(-1)
    TEMP3 <- (4 + z^2)^(-1)
    TEMP4 <- (4 + w^2)^(-1)
    
    v <- TEMP1 * TEMP2 - TEMP3 * TEMP4
    v1 <- sigma * (-2 * x * TEMP1^2 * TEMP2 + 2 * w * TEMP4^2 * TEMP3)
    v2 <- sigma * (-2 * y * TEMP1 * TEMP2^2 + 2 * z * TEMP4 * TEMP3^2)
    v11 <- sigma^2 * (-(6 * x^2 - 8) * TEMP1^3 * TEMP2 + 4 * z * w * TEMP3^2 * TEMP4^2)
    v22 <- sigma^2 * (-(6 * y^2 - 8) * TEMP1 * TEMP2^3 + 4 * z * w * TEMP3^2 * TEMP4^2)
    v12 <- sigma^2 * (-4 * x * y * TEMP1^2 * TEMP2^2 + (6 * w^2 - 8) * TEMP4^3 * TEMP3)
    v21 <- sigma^2 * (-4 * x * y * TEMP1^2 * TEMP2^2 + (6 * z^2 - 8) * TEMP3^3 * TEMP4)
  }
  
  vec <- c(v, v1, v2, v11, v22, v12, v21)
  return(vec)
}
