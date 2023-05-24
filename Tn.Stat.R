Tn.Stat <- function(X){
  n <- nrow(X)
  U <- pobs(X)
  u <- matrix(c(rep((rep(1:n)-0.5)/n, each=n), rep((rep(1:n)-0.5)/n, n)), ncol =2)
  Cn <- matrix(C.n(cbind(u[, 1], u[, 2]),U), n, n)
  Tn <- as.vector(max(abs(Cn-t(Cn))))
  return(Tn*sqrt(n))
}
