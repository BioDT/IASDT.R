#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
arma::mat Solve1(const arma::mat & x) {
    return arma::inv(x);
}

//' @export
// [[Rcpp::export]]
arma::mat Solve2(const arma::mat& A, const arma::mat& B) {
    return arma::solve(A, B);
}

//' @export
// [[Rcpp::export]]
arma::vec Solve2vect(const arma::mat& A, const arma::vec& B) {
    arma::vec x = arma::solve(A, B);
    // Ensure the result is a vector
    return x;
}

//' @export
// [[Rcpp::export]]
NumericVector fast_pnorm(NumericVector x) {
    int n = x.size();
        NumericVector result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = 0.5 * erfc(-x[i] / sqrt(2.0));
    }
    return result;
}

//' @export
// [[Rcpp::export]]
arma::mat exp_neg_div(const arma::mat& D11, double x) {
    return arma::exp(-D11 / x);
}
