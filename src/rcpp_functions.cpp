#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' @title helper C++ functions for fast matrix computations
//' 
//' @description
//' A collection of efficient C++ functions using RcppArmadillo for common matrix operations, including solving linear systems, computing matrix inverses, approximating normal CDFs, and fast elementwise transformations.
//' 
//' @name cpp_functions
//' @rdname cpp_functions
//' @param A A numeric matrix (for solving, exponential operations, etc.).
//' @param B A numeric matrix or vector (right-hand side for solving linear systems).
//' @param x A numeric matrix, numeric vector, or numeric scalar (depending on the function).
//' 
//' @return 
//' - `solve1()`: A numeric matrix, the inverse of `x`.
//' - `solve2()`: A numeric matrix, the solution to `A * X = B`.
//' - `solve2vect()`: A numeric vector, the solution to `A * x = B`.
//' - `fast_pnorm()`: A numeric vector, CDF values approximated for standard normal distribution.
//' - `exp_neg_div()`: A numeric matrix, elementwise exponential of `-A/x`.
//' 
//' @author Ahmed El-Gabbas
//' @examples
//' 
//' # -----------------------------------------
//' # Example for solve1
//' # -----------------------------------------
//' 
//' N <- 100
//' set.seed(1000)
//' Matrix <- matrix(rnorm(N * N), N, N)
//' all.equal(solve(Matrix), IASDT.R::solve1(Matrix))
//' 
//' # -----------------------------------------
//' # Example for solve2
//' # -----------------------------------------
//' 
//' N <- 100
//' set.seed(1000)
//' A <- matrix(rnorm(N * N), N, N)
//' set.seed(2000)
//' B <- matrix(rnorm(N * N), N, N)
//' identical(solve(A, B), IASDT.R::solve2(A, B))
//' 
//' set.seed(2000)
//' B <- matrix(rnorm(N), N, 1)
//' identical(solve(A, B), IASDT.R::solve2(A, B))
//' 
//' # -----------------------------------------
//' # Example for solve2vect
//' # -----------------------------------------
//' 
//' N <- 100
//' set.seed(1000)
//' A <- matrix(rnorm(N * N), N, N)
//' set.seed(2000)
//' B <- rnorm(N)
//' identical(solve(A, B), as.vector(IASDT.R::solve2vect(A, B)))
//' 
//' # -----------------------------------------
//' # Example for fast_pnorm
//' # -----------------------------------------
//' 
//' set.seed(1000)
//' A <- rnorm(100)
//' all.equal(pnorm(A), IASDT.R::fast_pnorm(A))
//' 
//' # -----------------------------------------
//' # Example for exp_neg_div
//' # -----------------------------------------
//' 
//' N <- 1000
//' set.seed(1000)
//' A <- matrix(rnorm(N * N), N, N)
//' set.seed(2000)
//' x <- rnorm(1)
//' identical(exp(-A / x), IASDT.R::exp_neg_div(A, x))

//' @rdname cpp_functions
//' @export
// [[Rcpp::export]]
arma::mat solve1(const arma::mat& x) {
    return arma::inv(x);
}

//' @rdname cpp_functions
//' @export
// [[Rcpp::export]]
arma::mat solve2(const arma::mat& A, const arma::mat& B) {
    return arma::solve(A, B);
}

//' @rdname cpp_functions
//' @export
// [[Rcpp::export]]
arma::vec solve2vect(const arma::mat& A, const arma::vec& B) {
    arma::vec x = arma::solve(A, B);
    return x;
}

//' @rdname cpp_functions
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

//' @rdname cpp_functions
//' @export
// [[Rcpp::export]]
arma::mat exp_neg_div(const arma::mat& A, double x) {
    return arma::exp(-A / x);
}
