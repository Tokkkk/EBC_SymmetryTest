source("BerEmpCopula.R")

Tnm.Stat <- function(X, M){
  U <- pobs(X)
  n <- nrow(X)
  u <- matrix(c(rep((rep(1:n)-0.5)/n, each=n), rep((rep(1:n)-0.5)/n, n)), ncol=2)
  Cnm <- matrix(BerEmpCopula(u,U,M), nrow=n, ncol=n)
  Tnm <- as.vector(max(abs(Cnm-t(Cnm))))
  return(Tnm*sqrt(n))
}
