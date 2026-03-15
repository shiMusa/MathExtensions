// ? LAPACK *gemv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga0983da08821bec7701e90fb1e65c8cd7.html#ga0983da08821bec7701e90fb1e65c8cd7


// ? original comment, cgemv()
// CGEMV performs one of the matrix-vector operations
//
//    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
//
//    y := alpha*A**H*x + beta*y,
//
// where alpha and beta are scalars, x and y are vectors and A is an
// m by n matrix.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
//  
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

sgemv :: _gemv;
dgemv :: _gemv;
cgemv :: _gemv;
zgemv :: _gemv;

_gemv :: (
    $$trans: u8,         // normal (N), transposed (T), or hermitian conjugated (H)
    $$m: int,            // rows of a
    $$n: int,            // columns of a
    alpha: $T,           // scalar
    a: $A/ArrayType2D,   // matrix A 
    $$lda: int,          // first dimension of A  // unused
    x: $X/ArrayType1D,   // vector
    $$incx: int,         // increment of elements of x
    beta: T,             // scalar
    y: *$Y/ArrayType1D,  // vector
    $$incy: int          // increment of elements of y
) {
    DEBUG :: false;
    #if DEBUG println_push("_gemv");

    cone  := one(T);
    czero := zero(T);

    temp : T = ---;
    ix, iy, jx, jy, lenx, leny: int = ---;

    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == X.data_type);
            assert(T == Y.data_type);
        }
        assert(char_eq(trans, #char "n", #char "t", #char "c"));
        assert(m >= 0);
        assert(n >= 0);
        assert(lda >= max(1,m)); // 1 ?
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible

    if m == 0 || n == 0 || (alpha == czero && beta == cone) then return;

    noconj := char_eq(trans, #char "t");

    // set lenx and leny, the lengths of the vectors x and y, and set
    // up the start points in x and y.

    if char_eq(trans, #char "n") {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }
    kx := ifx incx > 0 then 1 else 1 - (lenx-1)*incx;
    ky := ifx incy > 0 then 1 else 1 - (leny-1)*incy;

    // start the operations. In this version the elements of a are
    // accessed sequantially with cone pass through a.

    // first form y := beta*y

    #if DEBUG println("beta == czero ? %", beta == czero, color=.FG_CYAN);
    if beta != cone {
        if incy == 1 {
            if beta == czero {
                for i:1..leny {
                    fset(y, i, czero);
                }
            } else {
                for i:1..leny {
                    fset(y, i, beta*fget(y.*, i));
                }
            }
        } else {
            iy = ky;
            if beta == czero {
                for i:1..leny {
                    fset(y, iy, czero);
                    iy += incy;
                }
            } else {
                for i:1..leny {
                    fset(y, iy, beta*fget(y.*, iy));
                    iy += incy;
                }
            }
        }
    }
    
    #if DEBUG println("alpha == czero ? %", alpha == czero, color=.FG_CYAN);
    if alpha == czero then return;

    #if DEBUG println("trans == 'n' ? %", char_eq(trans, #char "n"), color=.FG_CYAN);
    if char_eq(trans, #char "n") {
        // form y := alpha*A**x + y
        jx = kx;
        if incy == 1 {
            for j:1..n {
                temp = alpha * fget(x, jx);
                for i:1..m {
                    fset(y, i, fget(y.*, i) + temp * fget(a, i, j));
                }
                jx += incx;
            }
        } else {
            for j:1..n {
                temp = alpha * fget(x, jx);
                iy = ky;
                for i:1..m {
                    fset(y, iy, fget(y.*, iy) + temp * fget(a, i, j));
                    iy += incy;
                }
                jx += incx;
            }
        }
    } else {
        // form y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
        jy = ky;
        if incx == 1 {
            for j:1..n {
                temp = czero;
                if noconj {
                    for i:1..m {
                        temp += fget(a, i, j)*fget(x, i);
                    }
                } else {
                    for i:1..m {
                        temp += conjugate(fget(a, i, j))*fget(x, i);
                    }
                }
                fset(y, jy, fget(y.*, jy) + alpha*temp);
                jy += incy;
            }
        } else {
            for j:1..n {
                temp = czero;
                ix = kx;
                if noconj {
                    for i:1..m {
                        temp += fget(a, i, j)*fget(x, ix);
                        ix += incx;
                    }
                } else {
                    for i:1..m {
                        temp += conjugate(fget(a, i, j))*fget(x, ix);
                        ix += incx;
                    }
                }
                fset(y, jy, fget(y.*, jy) + alpha*temp);
                jy += incy;
            }
        }
    }
}
