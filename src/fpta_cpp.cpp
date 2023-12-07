#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat gram_schimdtCpp(arma::mat A) {
  int n = A.n_cols; //number of basis
  int m = A.n_rows; //number of datapoints
  arma::mat Q(m, n);
  Q.fill(0);
  arma::mat R(n, n);
  R.fill(0);
  for (int j = 0; j < n; j++) {
    arma::vec v = arma::vec(A.col(j));
    arma::vec coef_v = arma::vec(n);
    coef_v(j) = 1;
    if (j > 0) {
      for(int i = 0; i < j; i++) {
        double coef = arma::dot(A.col(j), Q.col(i))/m;
        v = v - coef * Q.col(i);
        coef_v -= coef * R.col(i);
      }
    }
    double norm = std::sqrt(arma::dot(v, v)/m);
    Q.col(j) = v / norm;
    R.col(j) = coef_v / norm;
  }
  return R;
}

//[[Rcpp::export]]
Rcpp::List Schur_wrapperCpp(arma::mat projection) {
  arma::mat U;
  arma::mat S;
  arma::schur(U, S, projection);
  return Rcpp::List::create(Rcpp::Named("U") = U,Rcpp::Named("S")=S);
}


//[[Rcpp::export]]
arma::vec predictCpp(arma::mat B1, arma::mat B2, arma::mat Coef){
  arma::mat temp1 = B1 * Coef.t();
  arma::mat temp2 = B2 * Coef.t();
  int p = Coef.n_rows;
  int n = B1.n_rows;
  arma::vec y(n);
  y.fill(0);
  for(int i = 0; i < p /2; i++){
    y += temp1.col(2*i) % temp2.col(2*i + 1) - temp1.col(2*i+1) % temp2.col(2*i);
  }
  return y;
}
