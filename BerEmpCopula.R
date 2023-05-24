########################################
#Empirical Bernstein copula
########################################
library("Rcpp")
empCopGrid <-  function (U,g)  
{  
  grid <- seq(0, 1, by=1/g)
  u =cbind(rep(grid, each=g+1),rep(grid,g+1))
  temp = C.n(u,U)
  C <- matrix(temp,g+1,g+1)
  return(C)
}


BerEmpCopula <- function (u,U,k) #U is pseudo-observations, K is the order
{
  n = nrow(u)
  Ber = rep(0, n)
  emp = empCopGrid(U,k) # (K+1)*(K+1) matrix
  
  # bin1 = lapply(u[,1], function(t){dbinom(0:k, k, t)})
  # bin2 = lapply(u[,2], function(t){dbinom(0:k, k, t)})
  # Ber = mapply(function(r, s){t(r)%*%emp%*%s}, bin1, bin2)

  for(i in 1:n)
  {
    bin1 = t(dbinom(0:k, k, u[i,1]))
    bin2 = dbinom(0:k, k, u[i,2])
    Ber[i] = bin1%*%emp%*%bin2
  }
  
  return(Ber)
}

#######################################################
#Empirical Bernstein copula partial derivative
#######################################################
empCopGrid1 <-  function (U,g)  
{  
  u11 = cbind(rep(seq(1/g, 1, by=1/g), each=g+1),rep(seq(0, 1, by=1/g),g))
  u12 = cbind(rep(seq(0, (g-1)/g, by=1/g), each=g+1),rep(seq(0, 1, by=1/g),g))
  u21 = cbind(rep(seq(0, 1, by=1/g), each=g),rep(seq(1/g, 1, by=1/g),g+1))
  u22 = cbind(rep(seq(0, 1, by=1/g), each=g),rep(seq(0, (g-1)/g, by=1/g),g+1))
  temp11 = C.n(u11,U);temp12 = C.n(u12,U);temp21 = C.n(u21,U);temp22 = C.n(u22,U)
  C11 <- matrix(temp11,g+1,g); C12 <- matrix(temp12,g+1,g) # k is column, ell is row
  C21 <- matrix(temp21,g,g+1); C22 <- matrix(temp22,g,g+1)
  return(list(C11, C12, C21, C22))
}


EBC_deriv <- function(u, U, k){   #u is a matrix
  n = nrow(u)
  Ber = rep(0, n)
  emp = empCopGrid1(U, k)
  bin11 = sapply(u[,1], function(t){dbinom(0:(k-1), k-1, t)})
  bin12 = sapply(u[,2], function(t){dbinom(0:k, k, t)})
  bin21 = sapply(u[,1], function(t){dbinom(0:k, k, t)})
  bin22 = sapply(u[,2], function(t){dbinom(0:(k-1), k-1, t)})
  der1 = k*diag(t(bin12)%*%(emp[[1]]-emp[[2]])%*%bin11)
  der2 = k*diag(t(bin22)%*%(emp[[3]]-emp[[4]])%*%bin21)
  return(cbind(der1, der2))
}



























