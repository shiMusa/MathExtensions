# MathExtensions

## Advanced Math functions for Jai

The goal of this project is to progressively write a comprehensive
mathematical library for Jai, similar to [BLAS](http://www.netlib.org/blas/), the Gnu Scientific Library, or [Hipparchus](https://hipparchus.org/) (former Java Apache Math).

The focus lies in providing the functionality first, performance
second, but no unnecessary performance costs if possible.

Feel free to add functionality and performance upgrades.

## Current thoughts

There was a major rewite. Now, almost all functions do **not** allocate, but rather the results are stored in arguments of the functions themselves. This improves some of the linear algebra functions drastically.
It's basically the old-school way of doing things, and there was a reason to do so.
The new (old-style) system is more general and performant. 

To enable future optimizations, every funtion is basically just a relay to specialized functions, so `foo(...)` will actually call `foo_default(...)`, `foo_dense(...)`, etc. Currently, most functions have a default implementation for `$V/VectorType` or `$M/MatrixType`. It should be ease to extend the system with specialized, high-performance functions for special matrix/vector types (e.g. sparse, triagonal, hermitian, ...).

Almost all operator overloads were removed in that process.
The next step is, to re-write the operator overloads for special matrix and vector types. In these specialized functions, allocation can be handeled efficiently.

Generally:
- Current work focuses on making everything as generic as possible to enable different types of matrices/vector, heap/stack allocated, various number types.
- I'm going through [2] now to improve the algorithms since that book is actually considering special properties of matrices early on and also writes out EVERY algorithm used.
(- I might use some [BLAS](http://www.netlib.org/blas/) naming conventions and functions)
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
- vector: `VectorType`; real or complex vector

    The following types are implemented:
    - `DenseVector(Type, int)`
    - `DenseHeapVector(Type)`
    - `VectorView(Type, int)` (behaves like `DenseVector`)
    - `VectorHeapView(Type)` (behaves like `DenseHeapVector`)
    - `MatrixRowView(Type, int)` (behvaes like `DenseVector`)
    - `MatrixRowHeapView(Type)` (behvaes like `DenseHeapVector`)
    - `MatrixColumnView(Type, int)` (behaves like `DenseVector`)
    - `MatrixColumnHeapView(Type)` (behaves like `DenseHeapVector`)
    

    The following functions are implemented:
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

    The following types are implemented:
    - `DenseMatrix(Type, int, int)`
    - `DenseHeapMatrix(Type)`
    - `MatrixView(Type, int, int)` (behaves like `DenseMatrix`)
    - `MatrixHeapView(Type)` (behaves like `DenseHeapMatrix`)

    The following functions are implemented:
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
