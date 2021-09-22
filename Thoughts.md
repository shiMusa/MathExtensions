
# Thoughts

This is a list of plans and thoughts to consider when developing this math library.

We can draw a lot of inspiration from the [blaze C++ library](https://bitbucket.org/blaze-lib/blaze/src/master/).

## Generalization

Planned specializations are:

- `TriangularMatrix`
- `BandMatrix`
- `NeighbourMatrix`; name not final; similar to `BandMatrix` but with periodic boundary conditions
- `SparseMatrix`, `SparseVector`
- and more to come...


### Hybrid Matrices

There should also be a possibility to combine these optimized matric forms to hybrid/mixed matrices.
It's not clear yet how to do that. On first thought, we might do sth like
```
HybridMatrix :: struct(a: A/Matrix, b: B/Matrix) {...}
```
but need also to make sure that `HybridMatrix` fits `T/Matrix` for all the other algorithms;

Actually, I'm not sure, if this makes sense at all...



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