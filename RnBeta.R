
RnBeta1 = function(n, B, R1, R2){
  B1 = matrix(0, n, n)
  for(i in 1:n){
    for (j in 1:n) {
      B1[j, i] = B[R1[i], R1[j]] * B[R2[i], R2[j]] - B[R1[i], R2[j]] * B[R2[i], R1[j]]
    }
  }
  return(sum(B1)*2/(n^2))
}


