# MathExtensions

## Advanced Math functions for Jai

The goal of this project is to progressively write a comprehensive
mathematical library for Jai, similar to the Gnu Scientific Library or Java Apache Math.

The focus lies in providing the functionality first, performance
second. 

Feel free to add functionality and performance upgrades.

## Current thoughts

- I just refactored the code to use generic types, e.g. `Matrix($T, $N, $M)`. The next step is to introduce a matrix-interface and make all the algorithms fully generic. `Prototype-Mk-I.jai` shows how to do it.
- Add `Quaternion64` (maybe even `Octonion64`)
- I'm going through [2] now to improve the algorithms since that book is actually considering special properties of matrices early on and also writes out EVERY algorithm used.

There are more thoughts written down in [the document outlining some of my thoughts](https://github.com/shiMusa/MathExtensions/blob/flags-and-other-matrix-types/Thoughts.md).

## Currently Implemented

Complex and reals numbers, vectors, matrices interoperate automatically (in most cases).
Many functions specialize during compile-time on the real or complex variant (e.g. all the elementary functions).

- Utils
    - `calculation(Code)`; will create and release a Pool of memory for all the allocations in the Code block.
- Complex numbers `+`, `-`, `*`, `/`,
    - `str`; for pretty printing
- Elementary (`z` ∈ ℂ)
    - `sign(r)`
    - `factorial(z)`
    - `binomial(from, choose)`
    - `conjugate(z)`
    - `abs_sq(z)`
    - `abs(z)`
    - `arg(z)`
    - `exp(z)`
    - `log(z)`
    - `pow(z)`
    - `sqrt(z)`
    - `phase(magnitude ∈ ℂ, angle ∈ ℂ)`
    - `sin`, `cos`, `tan`, `cot`, `sec`, `csc` on ℂ
    - `asin`, `acos`, `atan`, `acot`, `asec`, `acsc` on ℂ
    - `sinh`, `cosh`, `tanh`, `coth`, `sech`, `csch` on ℂ
    - `asinh`, `acosh`, `atanh`, `acoth`, `asech`, `acsch` on ℂ
- polynomials
    - `solve_quadratic` -> ℂ
    - `solve_quadratic_real` -> ℝ
    - `polynom(x, ..a)` = a[0] x^n + a[1] x^{n-1} + ... + a[n] x^0, with a,x ∈ ℝ,ℂ
    - `synthetic_division`
    - `repeated_synthetic_division`
- vector: `Vector($Type, $Dim)`; real or complex vector with
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
- matrix: `Matrix($Type, $Rows, $Cols)`; real or complex
    - `str`; for pretty printing
    - operators [], `[][]`, `==`, `+`, `-`, `*`, `/`
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

- general
    - due to the constant work and improvment, many functions could be optimized using e.g. in-place functions instead of allocating data.
- linear algebra
    - LAPACK Scaling for Gaussian factorization (Algorithm 3.9.1 [1])
    - Iterative Improvement for Gaussian factorization (Algorithm 3.9.2 [1])


## Structure

Let's keep it as simple as possible:

- Operator overloading where it makes sense (linear algebra, complex numbers, quaternions, etc.), simple function calls otherwise.
- No special type-aliases if not stricly necessary
- A file for each major topic.

## Important

Give sources for algorithms written so that we can take a look and help debugging.

Write (at least some) tests contained in each file.


## Ressources

[1] Scientific Computing, Vol I: Linear and nonlinear equations, Texts in computational science and engineering 18, Springer

[2] Matrix Computations, 4th Edition; Gene H. Golub & Charles F. Van Loan; Johns Hopkins University Press, Baltimore; 2013