BerBivDerivatives.Discretized <- function(X,K,Var,M){
  library(subcopem2D)
  Y <- matrix(c(rep((rep(1:K)-0.5)/K, each=K), rep((rep(1:K)-0.5)/K, K)), ncol =2)
  A <- rep(0, K*K)
  if (Var == 1){
    for (i in 1:(K*K)) {
      A[i]<- Bcopula(X, m=M, TRUE)$du (Y[i, 1], Y[i, 2])  
    }
  }else{
    for (i in 1:(K*K)) {
      A[i]<- Bcopula(X, m=M, TRUE)$dv (Y[i, 1], Y[i, 2])  
    }
  }
  Chat <- matrix(A, byrow = T, ncol = K)
  return(Chat)
}