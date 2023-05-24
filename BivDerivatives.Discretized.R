# Description
#   Returns the estimation of the partial derivatives of a 
#   bivariate copula using the finite difference estimator
# Input
#   U: n x 2 matrix of standardized ranks
#   b: integer acting as a smoothing parameter
#   K: Discretization parameter of the grid in [0,1]^2
#   Var: =1 for the 1st component, =2 for the 2nd component
# Output
#   Chat: K x K matrix of the estimation of the partial derivative
#   w/r to component 'Var' at (i-.5,j-.5), i,j=1,...,K 
         

BivDerivatives.Discretized <- function(U,b,K,Var){
  n <- nrow(U)
  bn <- b/sqrt(n)
  
  Chat <- matrix(0, K, K)
  for(a in 1:K){
    for(b in 1:K){
      u <- (c(a, b)-0.5)/K
      
      if (u[Var]<bn) {
      v1 <- u;v1[Var]=2*bn
      v2 <- u;v2[Var]=0
      } 
      else if (u[Var]>1-bn) {
       v1 <- u;v1[Var] <- 1
       v2 <- u;v2[Var] <- 1-2*bn
      } 
      else {
        v1 <- u;v1[Var] <- u[Var]+bn
        v2 <- u; v2[Var] <- u[Var]-bn
      }
      
      SUM <- 0
      for (i in 1:n) {
        SUM <- SUM + prod(as.double(U[i,] <= v1)) -
          prod(as.double(U[i,] <= v2))
      }
      Chat[a, b] <- SUM/(2*n*bn)
    }
  }
  return(Chat)
}

