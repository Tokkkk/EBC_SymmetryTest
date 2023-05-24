Bn = function(n, N){
  u = seq(0.5/N, (N-0.5)/N, by=1/N)
  rank_index = cbind(rep(1:n, each=n), rep(1:n, times=n)) 
  A = apply(rank_index, 1, function(t){pbeta(u, t[1], n+1-t[1]) * pbeta(u, t[2], n+1-t[2])})
  A1 = colMeans(A)
  B_n = matrix(A1, nrow=n, byrow=T)
}