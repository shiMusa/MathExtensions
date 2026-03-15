// ? LAPACK *hemv() and *symv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gac23d51cad7f51484371119bb6e3fb1f3.html#gac23d51cad7f51484371119bb6e3fb1f3

// ? original comment
// CHEMV  performs the matrix-vector  operation
//
//    y := alpha*A*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n hermitian matrix.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.


// symv = symmetric matrix * vector
// hemv = hermitian matrix * vector
ssymv :: _hemv; // conjugate() has no effect on real numbers
dsymv :: _hemv; // conjugate() has no effect on real numbers
chemv :: _hemv;
zhemv :: _hemv;


_hemv :: (
    $$uplo: u8,
    $$n: int,
    alpha: $T,
    a: $A/ArrayType2D,
    $$lda: int,
    x: $X/ArrayType1D,
    $$incx: int,
    beta: T,
    y: *$Y/ArrayType1D,
    $$incy: int
) {
    _zero :: #run zero(T);
    _one  :: #run one(T);

    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == X.data_type);
            assert(T == Y.data_type);
        }
        assert(char_eq(uplo, #char "u", #char "l"));
        assert(n >= 0);
        assert(lda >= max(1,n)); // 1 offset?
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible
    if n == 0 || (alpha == _zero && beta == _one) then return;

    // set up the start points in x and y.
    kx := ifx incx > 0 then 1 else 1 - (n-1)*incx;
    ky := ifx incy > 0 then 1 else 1 - (n-1)*incy;

    temp1, temp2 : T = ---;
    ix, iy, jx, jy : int = ---;

    // Start the operations. In this version the elements of A are
    // accessed sequentially with _one pass through the triangular part

    // first form  y := beta*y
    if beta != _one {
        if incy == 1 {
            if beta == _zero {
                for i:1..n {
                    fset(y, i, _zero);
                }
            } else {
                for i:1..n {
                    fset(y, i, beta*fget(y.*, i));
                }
            }
        } else {
            iy = ky;
            if beta == _zero {
                for i:1..n {
                    fset(y, iy, _zero);
                    iy += incy;
                }
            } else {
                for i:1..n {
                    fset(y, iy, beta*fget(y.*, iy));
                    iy += incy;
                }
            }
        }
    }   
    if alpha == _zero then return;

    if char_eq(uplo, #char "u") {
        // Form y when A is stored in the upper triangle.
        if incx == 1 && incy == 1 {
            for j:1..n {
                temp1 = alpha*fget(x, j);
                temp2 = _zero;
                for i:1..j-1 {
                    fset(y, i, fget(y.*, i) + temp1*fget(a, i, j));
                    temp2 += conjugate(fget(a, i, j))*fget(x, i);
                }
                fset(y, j, fget(y.*, j) + temp1*real(fget(a, j, j)) + alpha*temp2);
            }
        } else {
            jx = kx;
            jy = ky;
            for j:1..n {
                temp1 = alpha*fget(x, jx);
                temp2 = _zero;
                ix = kx;
                iy = ky;
                for i:1..j-1 {
                    fset(y, iy, fget(y.*, iy) + temp1*fget(a, i, j));
                    temp2 += conjugate(fget(a, i, j))*fget(x, ix);
                    ix += incx;
                    iy += incy;
                }
                fset(y, jy, fget(y.*, jy) + temp1*real(fget(a, j, j)) + alpha*temp2);
                jx += incx;
                jy += incy;
            }
        }
    } else {
        // Form y when A is stored in the lower triangle.
        if incx == 1 && incy == 1 {
            for j:1..n {
                temp1 = alpha*fget(x, j);
                temp2 = _zero;
                fset(y, j, fget(y.*, j) + temp1*real(fget(a, j, j)));
                for i:j+1..n {
                    fset(y, i, fget(y.*, i) + temp1*fget(a, i, j));
                    temp2 += conjugate(fget(a, i, j))*fget(x, i);
                }
                fset(y, j, fget(y.*, j) + alpha*temp2);
            }
        } else {
            jx = kx;
            jy = ky;
            for j:1..n {
                temp1 = alpha*fget(x, jx);
                temp2 = _zero;
                fset(y, jy, fget(y.*, jy) + temp1*real(fget(a, j, j)));
                ix = jx;
                iy = jy;
                for i:j+1..n {
                    ix += incx;
                    iy += incy;
                    fset(y, iy, fget(y.*, iy) + temp1*fget(a, i, j));
                    temp2 += conjugate(fget(a, i, j))*fget(x, ix);
                }
                fset(y, jy, fget(y.*, jy) + alpha*temp2);
                jx += incx;
                jy += incy;
            }
        }
    }
}
