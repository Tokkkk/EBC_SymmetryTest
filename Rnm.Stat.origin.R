source("BerEmpCopula.R")

Rnm.Stat <- function(X,K,M){
  U <- pobs(X)
  n <- nrow(X)
  w = matrix(c(rep((rep(1:K)-0.5)/K, each=K), rep((rep(1:K)-0.5)/K, K)), ncol =2)
  Rnm <- (BerEmpCopula(w, U, M)-BerEmpCopula(cbind(w[, 2], w[, 1]), U, M))^2
  Rn <- (C.n(cbind(w[, 1], w[, 2]), U)-C.n(cbind(w[, 2], w[, 1]), U))^2
  return(c(n*mean(Rnm), n*mean(Rn)))
}


