# MathExtensions

## Advanced Math functions for Jai

The goal of this project is to progressively write a comprehensive
mathematical library for Jai, similar to the Gnu Scientific Library or Java Apache Math.

The focus lies in providing the functionality first, performance
second. 

Feel free to add functionality and performance upgrades.

## Current thoughts

- I have to rethink the focus of the implementations: 
  1. remove all operators and all functions that return a concrete type
  2. remove all automatic casting in the remaining functions
  3. write default implementations (slow, but they work)
  4. specialize to optimize later
  5. generate functions with return parameters maybe with compile-time code generation
  6. use the latter for operator overloading


- Current work focuses on making everything as generic as possible to enable different types of matrices/vector, heap/stack allocated, various number types.
- I'm going through [2] now to improve the algorithms since that book is actually considering special properties of matrices early on and also writes out EVERY algorithm used.
- Implement _views_ for vectors and matrices which themselves behave like vectors and matrices (they implement `VectorType` and `MatrixType`).
- When generalizing the `MatrixType` and `VectorType` such that the dimensions are run-time variables (usefull for lots of applications), we cannot use helper types for the return parameters. Therefore, the operators `+ - / *` etc. have to be individually overloaded while the functions `add sub mul div` etc. can be kept general since the result is stored in one of the arguments.
- I'll probably also use the [BLAS](http://www.netlib.org/blas/) naming convention and functions
(- maybe add `Octonion($T)` etc. later)

There are more thoughts written down in [the document outlining some of my thoughts](https://github.com/shiMusa/MathExtensions/blob/flags-and-other-matrix-types/Thoughts.md).

## Currently Implemented

Quaternion, complex and reals scalars, vectors, matrices interoperate automatically (in most cases).
Many functions specialize during compile-time on the real, complex, or quaternion variant (e.g. all the elementary functions).

- Utils
    - `pool(#code)`; will create and release a Pool of memory for all the allocations in the Code block.
- `Complex($T)` numbers `+`, `-`, `*`, `/`,
    - `str` for pretty printing
- `Quaternion($T)` numbers `+`, `-`, `*`, `/`,
    - `str` for pretty printing
    - [more stuff]
- Elementary (`n` ∈ ℂ, ℍ)
    - `sign(r)`
    - `factorial(n)`
    - `binomial(from, choose)`
    - `conjugate(n)`
    - `abs_sq(n)`
    - `abs(n)`
    - `arg(n)`
    - `exp(n)`
    - `log(n)`
    - `pow(n)`
    - `sqrt(n)`
    - `phase(magnitude ∈ ℂ, angle ∈ ℂ)` [need to update for quaternion]
    - `sin`, `cos`, `tan`, `cot`, `sec`, `csc` on ℂ, ℍ
    - `asin`, `acos`, `atan`, `acot`, `asec`, `acsc` on ℂ, ℍ
    - `sinh`, `cosh`, `tanh`, `coth`, `sech`, `csch` on ℂ, ℍ
    - `asinh`, `acosh`, `atanh`, `acoth`, `asech`, `acsch` on ℂ, ℍ
- polynomials
    - `solve_quadratic` -> ℂ
    - `solve_quadratic_real` -> ℝ
    - `polynom(x, ..a)` = a[0] x^n + a[1] x^{n-1} + ... + a[n] x^0, with a,x ∈ ℝ,ℂ
    - `synthetic_division`
    - `repeated_synthetic_division`
- vector: `VectorType`; real or complex vector with

    `VectorType` is a generic interface/trate that is implemented by any concrete vector struct, e.g. `DenseVector($Type, $Dimensions)`.
    - `str`; for pretty printing
    - operators `[]`, `==`, `+`, `-`, `*`, `/`
    - in-place functions add, sub, neg, mul, div
    - multiple initialization functions (ones, basis, varargs, etc.)
    - `conjugate`
    - `outer_product`
    - `reflect`
    - `norm(vec, n)` with n ∈ ℝ, specialisations `norm_2`, `norm_1`, `norm_inf`
    - `cross`, for 3-dim vectors
    - `angle`
    - `permute`
    - `swap`
- matrix: `MatrixType` ; real or complex

    `MatrixType` is a generic interface/trate that is implemented by any concrete matrix struct, e.g. `DenseMatrix($Type, $Columns, $Rows)`.
    - `pstr`, `str` for pretty printing
    - operators `[]`, `[][]`, `==`, `+`, `-`, `*`, `/`
    - `row`, `column`
    - in-place functions `add`, `sub`, `neg`, `mul`, `div`
    - multiple initialization functions (1, ones, hadamard, varargs, etc.)
    - `reflector`
    - `submatrix`
    - `transpose`
    - `conjugate`
    - `conjugate_transpose` = `dagger`
    - `tensor`
    - `norm_1`, `norm_inf`, `norm_frobenius`
    - `permute_rows`, `permute_columns`, `permute`
    - `swap_columns`, `swap_rows`
- checks
    - `is_diagonal_unit(M)`
    - `is_(left, right)\_triangular(M)`
    - `is_right_quasi_triangular(M)`
    - `is_unitary(M)`
    - `is_(left, right)_trapezoidal(M)`
- linear algebra
    - `solve_linear_2x2` = `sl2`
    - `solve_linear_(left, right)_triangular` = `sl(l, r)ta`
    - `solve_linear_right_quasi_triangular` = `slrqta`
    - `solve_linear_orthogonal_projection` = `solve_linear_unitary` = `slop`
    - `solve_linear_successive_orthogonal_projection` = `slsop`
    - `solve_linear_(left, right)_trapezoidal` = `sl(l, r)tz` for vectors and matrices
    - `decompose_LR`
    - `gaussian_factorization_(no, full, row)_pivot`
    - `solve_linear_gaussian_factorization_(no, full, row)_pivot` = `slgf_(n, f, r)p` = `solve_LR(n, f, r)p`
    - `inverse(M)`; M quadratic matrix ∈ ℂ
    - `determinant` = `det`; ∈ ℂ

### Future Tasks

- linear algebra
    - LAPACK Scaling for Gaussian factorization (Algorithm 3.9.1 [1])
    - Iterative Improvement for Gaussian factorization (Algorithm 3.9.2 [1])

## Important

Give sources for algorithms written so that we can take a look and help debugging.

Write (at least some) tests contained in each file.


## Ressources

[1] Scientific Computing, Vol I: Linear and nonlinear equations, Texts in computational science and engineering 18, Springer

[2] Matrix Computations, 4th Edition; Gene H. Golub & Charles F. Van Loan; Johns Hopkins University Press, Baltimore; 2013
