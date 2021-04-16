
# Thoughts

This is a list of plans and thoughts to consider when developing this math library.

We can draw a lot of inspiration from the [blaze C++ library](https://bitbucket.org/blaze-lib/blaze/src/master/).

## Generalization

Currently, there is `Matrix(T)` and `Vector(T)`. Both are _dense_ representations of data, meaning every entry has a place in memory.

Work in Progress is `TriangularMatrix(T)` that only stores the relevant data.

To enable these specializations, algorithms in general should require `$T/Matrix` and just use a common set of basic mathematical operations, such as `get(...)`, `set(...)`, `operator +/-/*//` etc.

Specialized functions can be overloaded to explicitly require e.g. `TriangluarMatrix(T)` for improved performance.

Other planned specializations are:
- `BandMatrix(T)`
- `NeighbourMatrix(T)`; name not final; similar to `BandMatrix` but with periodic boundary conditions
- `SparseMatrix(T)`, `SparseVector(T)`
- and more to come...

### Hybrid Matrices

There should also be a possibility to combine these optimized matric forms to hybrid/mixed matrices.
It's not clear yet how to do that. On first thought, we might do sth like
```
HybridMatrix :: struct(a: A/Matrix, b: B/Matrix) {...}
```
but need also to make sure that `HybridMatrix` fits `T/Matrix` for all the other algorithms;



## Evaluation

Currently, mathematical operation result in the end-type, e.g. 
```
operator * :: (a: Matrix($A), b: Matrix($B)) -> Matrix($C)
```
However, we might instead return a syntax-tree e.g.
```
... -> Operation(.M_M_MULTIPLY, *Matrix(A), *Matrix(B))
```
and using a function `eval(...)`, we'd actually perform the calculations, e.g.
```
expression := A * B * (a + b); // create AST
res := eval(expression); // perform calculations
```
This way, we might be able to optimize the order of operations performed, e.g.
```
(A * B) * (a + b) => Level 3,
A * (B * (a + b)) => Level 2.
```

At the same time, `eval(...)` could automatically use `pool(#code)` to do memory-management.
It might also be possible to choose from a wider array of elementary math functions in the optimization step, such as vectorization `AVX`, multithreading, GPU computation etc. All you'd have to setup is the correct `Engine` for eval.