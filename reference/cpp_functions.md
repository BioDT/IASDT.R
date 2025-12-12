# helper C++ functions for fast matrix computations

A collection of efficient C++ functions using RcppArmadillo for common
matrix operations, including solving linear systems, computing matrix
inverses, approximating normal CDFs, and fast elementwise
transformations.

## Usage

``` r
solve1(x)

solve2(A, B)

solve2vect(A, B)

fast_pnorm(x)

exp_neg_div(A, x)
```

## Arguments

- x:

  A numeric matrix, numeric vector, or numeric scalar (depending on the
  function).

- A:

  A numeric matrix (for solving, exponential operations, etc.).

- B:

  A numeric matrix or vector (right-hand side for solving linear
  systems).

## Value

- `solve1()`: A numeric matrix, the inverse of `x`.

- `solve2()`: A numeric matrix, the solution to `A * X = B`.

- `solve2vect()`: A numeric vector, the solution to `A * x = B`.

- `fast_pnorm()`: A numeric vector, CDF values approximated for standard
  normal distribution.

- `exp_neg_div()`: A numeric matrix, elementwise exponential of `-A/x`.

## Author

Ahmed El-Gabbas

## Examples

``` r
# -----------------------------------------
# Example for solve1
# -----------------------------------------

N <- 100
set.seed(1000)
Matrix <- matrix(rnorm(N * N), N, N)
all.equal(solve(Matrix), IASDT.R::solve1(Matrix))
#> [1] TRUE

# -----------------------------------------
# Example for solve2
# -----------------------------------------

N <- 100
set.seed(1000)
A <- matrix(rnorm(N * N), N, N)
set.seed(2000)
B <- matrix(rnorm(N * N), N, N)
all.equal(solve(A, B), IASDT.R::solve2(A, B))
#> [1] TRUE

set.seed(2000)
B <- matrix(rnorm(N), N, 1)
all.equal(solve(A, B), IASDT.R::solve2(A, B))
#> [1] TRUE

# -----------------------------------------
# Example for solve2vect
# -----------------------------------------

N <- 100
set.seed(1000)
A <- matrix(rnorm(N * N), N, N)
set.seed(2000)
B <- rnorm(N)
all.equal(solve(A, B), as.vector(IASDT.R::solve2vect(A, B)))
#> [1] TRUE

# -----------------------------------------
# Example for fast_pnorm
# -----------------------------------------

set.seed(1000)
A <- rnorm(100)
all.equal(pnorm(A), IASDT.R::fast_pnorm(A))
#> [1] TRUE

# -----------------------------------------
# Example for exp_neg_div
# -----------------------------------------

N <- 1000
set.seed(1000)
A <- matrix(rnorm(N * N), N, N)
set.seed(2000)
x <- rnorm(1)
all.equal(exp(-A / x), IASDT.R::exp_neg_div(A, x))
#> [1] TRUE
```
