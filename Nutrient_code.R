#############################################################
### This manuscript returns P-values of Table 6.
#############################################################


################################################
### Empirical tests and Bernstein tests
################################################
source("BerEmpCopula.R")
source("Rnm.Stat.origin.R")
source("BerBivDerivatives.Discretized1.R")
source("Tn.Stat.R")
source("BivDerivatives.Discretized.R")
source("Tnm.Stat.R")
library(copula)
library(VineCopula, warn.conflicts = FALSE)
library(Rcpp); sourceCpp("Cppfunctions.cpp")
library(profileR)

X <- cbind(nutrient$a, nutrient$c)
U <- pobs(X)
n=nrow(X)
M=floor(n/90) 
Mat = cbind(rep(0:M, each=M+1),rep(0:M, M+1))
N=20
K=5000


#-------------------------------------------------------------
Rn.Mul <- function(X,H,n,N){
  
  U <- pobs(X)   
  
  Chat1 = BivDerivatives.Discretized(U,b=1,N,1)
  Chat2 = BivDerivatives.Discretized(U,b=1,N,2)
  
  
  
  Q <- array(0, c(n, N, N))
  for (k in 1:N){
    for (ell in 1:N){
      u = k/N; v = ell/N
      for (i in 1:n){
        a1 = (U[i,1]<=u)*(U[i,2]<=v) - (U[i,1]<=v)*(U[i,2]<=u)
        a2 = Chat1[k,ell]*( (U[i,1]<=u) - (U[i,2]<=u) )
        a3 = Chat2[k,ell]*( (U[i,2]<=v) - (U[i,1]<=v) )
        Q[i, k, ell] <-  a1 - a2 - a3
      }
    }
  }
  
  
  Rn_hat <- matrix(0, H, 1)
  set.seed(1)
  A = rexp(n*H, 1)
  Delta = matrix(A, nrow=n, ncol=H)
  for (h in 1:H) {
    xi = Delta[, h]; Delta1 = xi-mean(xi)
    for (k in 1:N)  {
      for (ell in 1:N) {
        Rn_hat[h] = Rn_hat[h] + (Delta1%*%Q[ ,k,ell])^2 # Rn_hat[h] is used to do the summation for each k and ell
      }
    }
  }
  return(Rn_hat/(n*N^2))
}

Rnm.Mul <- function(X,H,n,N,M){
  U = pobs(X)
  grid = matrix(c(rep((rep(1:N)-0.5)/N, each=N), rep((rep(1:N)-0.5)/N, N)), ncol =2)
  grid_index = cbind(rep(1:N, each=N),rep(1:N, N))
  U0 = grid_index/N
  U1 = cbind(U0[, 1], 1)
  U2 = cbind(1, U0[, 2])  
  
  Chat1 = EBC_deriv(grid, U, M)[,1]
  Chat2 = EBC_deriv(grid, U, M)[,2]

  a1 = RA1(n, U0, U, Mat)
  a2 = RA2(n, U1, U, Chat1, Mat)
  a3 = RA3(n, U2, U, Chat2, Mat)
  Q = a1-a2-a3

  set.seed(1)
  A = rexp(n*H, 1)
  Delta = matrix(A, nrow=n, ncol=H)
  Delta1 = sweep(Delta, 2, colMeans(Delta))
  
  B = (Q %*% Delta1)^2
  Rnm_hat = colMeans(B)/n

  return(Rnm_hat)
}

Snm.Mul <- function(X,H,n,M){
  U = pobs(X)
  U1 = cbind(U[,1], 1)
  U2 = cbind(1, U[,2])  
  
  Chat1pseudos <- EBC_deriv(U, U, M)[,1]
  Chat2pseudos <- EBC_deriv(U, U, M)[,2]


  a1 = A1(n, U, U, Mat)
  a2 = A2(n, U1, U, Chat1pseudos, Mat)
  a3 = A3(n, U2, U, Chat2pseudos, Mat)
  Qpseudos = a1-a2-a3
  
  set.seed(1)
  A <- rexp(n*H, 1)
  Delta <- matrix(A, nrow=n, ncol=H)
  Delta1 <- sweep(Delta, 2, colMeans(Delta))
  Snm_rep <- apply((t(Delta1)%*%Qpseudos)^2, 1, sum)
  return(Snm_rep/n^2)
}

Tn.Mul <- function(X, H, n, N){
  U <- pobs(X)   
  
  Chat1 = BivDerivatives.Discretized(U,b=1,N,1)
  Chat2 = BivDerivatives.Discretized(U,b=1,N,2)
  
  Q <- array(0, c(N, N, n))
  for (i in 1:n){
    for (k in 1:N){
      for (ell in 1:N){
        u = k/N; v = ell/N
        a1 = (U[i,1]<=u)*(U[i,2]<=v) - (U[i,1]<=v)*(U[i,2]<=u)
        a2 = Chat1[k,ell]*( (U[i,1]<=u) - (U[i,2]<=u) )
        a3 = Chat2[k,ell]*( (U[i,2]<=v) - (U[i,1]<=v) )
        Q[k, ell, i] <-  a1 - a2 - a3
      }
    }
  }
  
  
  

  Tn_hat <- NULL
  set.seed(1)
  A = rexp(n*H, 1)
  Delta = matrix(A, nrow=n, ncol=H)
  for (h in 1:H){
    xi = Delta[, h]; Delta1 = xi-mean(xi)
    D_hat <- matrix(0, N, N)
    for (k in 1:N) {
      for (ell in 1:N) {
        D_hat[k,ell] = Delta1%*%Q[k,ell, ]
      }
    }
    Tn_hat[h] = max(abs(D_hat))
  }

  return(Tn_hat/sqrt(n))
}

Tnm.Mul <- function(X, H, n, N, M){
  U = pobs(X)
  grid = matrix(c(rep((rep(1:N)-0.5)/N, each=N), rep((rep(1:N)-0.5)/N, N)), ncol =2)
  grid_index = cbind(rep(1:N, each=N),rep(1:N, N))
  U0 = grid_index/N
  U1 = cbind(U0[, 1], 1)
  U2 = cbind(1, U0[, 2])  
  
  Chat1 = EBC_deriv(grid, U, M)[,1]
  Chat2 = EBC_deriv(grid, U, M)[,2]
  
  a1 = RA1(n, U0, U, Mat)
  a2 = RA2(n, U1, U, Chat1, Mat)
  a3 = RA3(n, U2, U, Chat2, Mat)
  Q = a1-a2-a3
  
  set.seed(1)
  A = rexp(n*H, 1)
  Delta = matrix(A, nrow=n, ncol=H)
  Delta1 = sweep(Delta, 2, colMeans(Delta))
  
  B = Q %*% Delta1
  Tnm_hat = apply(abs(B), 2, max)
  
  return(Tnm_hat/sqrt(n))
}

#---------------------------------------------------------
Rn_hat <- Rn.Mul(X, H=K, n, N)
pvRn <- mean(Rn_hat >= Rnm.Stat(X, N, M)[2])

Rnm_hat <- Rnm.Mul(X, H=K, n, N, M)
pvRnm <- mean(Rnm_hat >= Rnm.Stat(X, N, M)[1])

pvSn <- exchTest(X, N=K, m=0)$p.value

Snm_hat <- Snm.Mul(X, H=K, n, M)
pvSnm <- mean(Snm_hat >= Snm_stat(n, U, U, Mat))

Tn_hat <- Tn.Mul(X, H=K, n, N)
pvTn <- mean(Tn_hat >= Tn.Stat(X))

Tnm_hat <- Tnm.Mul(X, H=K, n, N, M)
pvTnm <- mean(Tnm_hat >= Tnm.Stat(X, M))

pv <- cbind(pvRn, pvRnm, pvSn, pvSnm, pvTn, pvTnm)
print(pv)

#####################################################
# empirical beta copula (smoothed beta bootstrap) from Kiriliouk2021
#####################################################
source("EmBetaBoot.R")
source("Matrix_Bn.R")
source("RnBeta.R")

B_n = Bn(n, N)
R1 = rank(X[, 1])
R2 = rank(X[, 2])
res = rep(0, n)
for (i in 1:K) {
  set.seed(i+1)
  A = EmBetaBoot(X1, n, n, R1, R2)
  res[i] = RnBeta1(n, B_n, rank(A[, 1]), rank(A[, 2]))
}

Stat = RnBeta1(n, B_n, rank(X[, 1]), rank(X[, 2]))
pvbeta = mean(res >= Stat)

print(pvbeta)






















