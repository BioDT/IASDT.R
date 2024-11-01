#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Solve1: Invert a matrix
//'
//' Computes the inverse of a square matrix using Armadillo's `inv` function.
//'
//' @param x A square matrix to be inverted.
//' @return The inverted matrix of the same dimensions as `x`.
//' //' @export
// [[Rcpp::export]]
arma::mat Solve1(const arma::mat & x) {
    return arma::inv(x);
}


//' Solve2: Solve a system of linear equations
//'
//' Solves a system of linear equations `A * X = B` for `X`, where `A` is a
//' square matrix and `B` is a matrix of compatible dimensions.
//'
//' @param A A square matrix representing the coefficients.
//' @param B A matrix or vector representing the constants.
//' @return The solution matrix or vector `X`.
//' @export
// [[Rcpp::export]]
arma::mat Solve2(const arma::mat& A, const arma::mat& B) {
    return arma::solve(A, B);
}



//' Solve2vect: Solve a system of linear equations for a vector
//'
//' Solves a system of linear equations `A * x = B` where `A` is a square
//' matrix and `B` is a vector.
//'
//' @param A A square matrix of coefficients.
//' @param B A vector of constants.
//' @return The solution vector `x`.
//' @export
// [[Rcpp::export]]
arma::vec Solve2vect(const arma::mat& A, const arma::vec& B) {
    arma::vec x = arma::solve(A, B);
    // Ensure the result is a vector
    return x;
}



//' fast_pnorm: Fast Normal CDF
//'
//' Computes the cumulative distribution function of the normal distribution
//' for a numeric vector using a fast approximation.
//'
//' @param x A numeric vector for which the CDF should be computed.
//' @return A numeric vector of CDF values for each element in `x`.
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



//' exp_neg_div: Exponential of Negative Division
//'
//' Computes the element-wise exponential of the negative division of each
//' element in a matrix by a scalar value.
//'
//' @param D11 A matrix of values to be divided and exponentiated.
//' @param x A scalar divisor.
//' @return A matrix of the same dimensions as `D11`, containing the computed
//' values.
//'
//' @export
// [[Rcpp::export]]
arma::mat exp_neg_div(const arma::mat& D11, double x) {
    return arma::exp(-D11 / x);
}
