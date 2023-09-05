#######################################
## This script is the test function for cf statistic from Bahraoui et. al 2018, 
## and this code is tranformed from Matlab code form https://oraprdnt.uqtr.uquebec.ca/pls/public/gscw031?owa_no_site=834
##############################################
# Input  -> X: n x 2 data matrix
#           weight: 'N' Normal, 'DE' Double-exponential
#           sigma: Smoothing parameter > 0
#           M: Number of multiplier bootstrap samples
# Output -> Stat: Test statistic
#           PV : P-value of the test
#           Test: '1' if Ho is rejected, '0' otherwise
# Necessary procedure: DiagSym_CopulaCf_B
####################################################

DiagSym_CopulaCf_Tests <- function(X, weight, sigma, M) {
  n <- nrow(X)
  U <- apply(X, 2, rank) / (n + 1)
  
  # Preliminary computations
  D0 <- matrix(0, n, n)
  D1 <- matrix(0, n, n)
  D2 <- matrix(0, n, n)
  D11 <- matrix(0, n, n)
  D22 <- matrix(0, n, n)
  D12 <- matrix(0, n, n)
  D21 <- matrix(0, n, n)
  I1 <- matrix(0, n, n)
  I2 <- matrix(0, n, n)
  
  for (j in 1:n) {
    for (j1 in 1:n) {
      vec <- DiagSym_CopulaCf_B(weight, sigma, U[j,], U[j1,])
      D0[j, j1] <- vec[1]
      D1[j, j1] <- vec[2]
      D2[j, j1] <- vec[3]
      D11[j, j1] <- vec[4]
      D22[j, j1] <- vec[5]
      D12[j, j1] <- vec[6]
      D21[j, j1] <- vec[7]
      I1[j, j1] <- (U[j, 1] <= U[j1, 1])
      I2[j, j1] <- (U[j, 2] <= U[j1, 2])
    }
  }
  
  # Test statistic
  Stat <- sum(D0) / n
  
  # Computation of the multiplier replicates
  T1 <- I1 %*% D1 + I2 %*% D2 + t(I1 %*% D1) + t(I2 %*% D2)
  T2 <- I1 %*% D11 %*% t(I1) + I1 %*% D12 %*% t(I2) + I2 %*% D21 %*% t(I1) + I2 %*% D22 %*% t(I2)
  D <- D0 + (T1 / n) + (T2 / n^2)
  

  Mult <- numeric(M)
  for (b in 1:M) {
    set.seed(b+1)
    xi <- rexp(n)
    Delta <- (xi / mean(xi)) - 1
    Delta_row <- matrix(Delta, nrow = 1)  # Reshape Delta as a row vector
    Mult[b] <- Delta_row %*% D %*% t(Delta_row) / n
  }
  
  PV <- sum(Mult > Stat) / M
  Test <- (PV < 0.05)
  
  return(list(Stat = Stat, PV = PV, Test = Test))
}
