#########################################
## This script returns the smoothed beta bootstrap from symmetric empirical beta copula.
## Input:
## X: the original data, a matrix with 2 columns;
## K: number of bootstrap replicates;
## n: number of rows in X;
## Output:
## V: an vector with K length. 
#########################################
library(extraDistr)
EmBetaBoot = function(X, K, n, R1, R2){
  I = rdunif(K, 1, n)
  V1 = sapply(R1[I], function(t)rbeta(1, t, n+1-t)); V2 = sapply(R2[I], function(t)rbeta(1, t, n+1-t))
  V = cbind(V1, V2)
  ic = sample(c(0, 1), K, replace=T, c(0.5, 0.5))
  V[ic==1, c("V1", "V2")] = V[ic==1, c("V2", "V1")] 
  return(V)
}
