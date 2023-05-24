#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]


// [[Rcpp::export]]
double bivQ(NumericVector u, NumericVector U, NumericMatrix M){
  double sumind = 0.0;
  int k = sqrt(M.nrow())-1;
  for (int i = 0; i < M.nrow(); ++i)
    sumind +=((U[0] <= M(i,0)/k) * (U[1] <= M(i,1)/k) - (U[0] <= M(i,1)/k) * (U[1] <= M(i,0)/k))* R::dbinom(M(i, 0), k, u[0], false) * R::dbinom(M(i, 1), k, u[1], false);
  return sumind;
}


// [[Rcpp::export]]

NumericVector A1(int n, NumericMatrix U, NumericMatrix V, NumericMatrix M){
  NumericMatrix v(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      v(j, i) = bivQ(U(i, _ ), V(j, _ ), M);
  return v;
}

// [[Rcpp::export]]

NumericVector A2(int n, NumericMatrix U, NumericMatrix V, NumericVector D, NumericMatrix M){
  NumericMatrix v(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      v(j, i) = bivQ(U(i, _ ), V(j, _ ), M) * D[i];
  return v;
}

// [[Rcpp::export]]

NumericVector A3(int n, NumericMatrix U, NumericMatrix V, NumericVector D, NumericMatrix M){
  NumericMatrix v(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      v(j, i) = bivQ(U(i, _ ), V(j, _ ), M) * D[i];
  return v;
}

// [[Rcpp::export]]

double Snm_stat(int n, NumericMatrix U, NumericMatrix V, NumericMatrix M){
  NumericVector v(n);
  NumericVector res(n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      v[i] += bivQ(U(i, _ ), V(j, _ ), M), res[i] = (v[i]/n) * (v[i]/n);
  return sum(res);    
}

// [[Rcpp::export]]

NumericVector RA1(int n, NumericMatrix U, NumericMatrix V, NumericMatrix M){
  double m = U.nrow();
  NumericMatrix v(m, n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      v(i, j) = bivQ(U(i, _ ), V(j, _ ), M);
  return v;
}

// [[Rcpp::export]]

NumericVector RA2(int n, NumericMatrix U, NumericMatrix V, NumericVector D, NumericMatrix M){
  double m = U.nrow();
  NumericMatrix v(m, n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      v(i, j) = bivQ(U(i, _ ), V(j, _ ), M) * D[i];
  return v;
}

// [[Rcpp::export]]

NumericVector RA3(int n, NumericMatrix U, NumericMatrix V, NumericVector D, NumericMatrix M){
  double m = U.nrow();
  NumericMatrix v(m, n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      v(i, j) = bivQ(U(i, _ ), V(j, _ ), M) * D[i];
  return v;
}











