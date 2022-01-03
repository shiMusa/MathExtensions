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
    $$alpha: $T,         // scalar
    $$a: $A/ArrayType2D, // matrix A
    $$lda: int,          // first dimension of A
    $$x: $X/ArrayType1D, // vector
    $$incx: int,         // increment of elements of x
    $$beta: T,           // scalar
    $y: *$Y/ArrayType1D, // vector
    $$incy: int          // increment of elements of y
) {
    one := one(T);
    zero := zero(T);

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

    if m == 0 || n == 0 || (alpha == zero && beta == one) then return;

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
    kx := ifx incx > 0 then 0 else (lenx-1)*(-incx);
    ky := ifx incy > 0 then 0 else (leny-1)*(-incy);

    // start the operations. In this version the elements of a are
    // accessed sequantially with one pass through a.

    // first form y := beta*y

    if beta != one {
        if incy == 1 {
            if beta == zero {
                for i:0..leny-1 {
                    set(y,i,zero);
                }
            } else {
                for i:0..leny-1 {
                    set(y,i,beta*get(y,i));
                }
            }
        } else {
            iy = ky;
            if beta == zero {
                for i:0..leny-1 {
                    set(y,iy, zero);
                    iy += incy;
                }
            } else {
                for i:0..leny-1 {
                    set(y,iy, beta*get(y,iy));
                    iy += incy;
                }
            }
        }
    }
    
    if alpha == zero then return;

    if char_eq(trans, #char "n") {
        // form y := alpha*A**x + y

        jx = kx;
        if incy == 1 {
            for j:0..n-1 {
                temp = alpha * get(x,jx);
                for i:0..m-1 {
                    set(y,i, get(y,i) + temp * get(a,i,j));
                }
                jx += incx;
            }
        } else {
            for j:0..n-1 {
                temp = alpha * get(x,jx);
                iy = ky;
                for i:0..m-1 {
                    set(y,iy, get(y,iy) + temp * get(a,i,j));
                    iy += incy;
                }
                jx += incx;
            }
        }
    } else {
        // form y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.

        jy = ky;
        if incx == 1 {
            for j:0..n-1 {
                temp = zero;
                if noconj {
                    for i:0..m-1 {
                        temp += get(a,i,j)*get(x,i);
                    }
                } else {
                    for i:0..m-1 {
                        temp += conjugate(get(a,i,j))*get(x,i);
                    }
                }
                set(y,iy, get(y,iy) + alpha*temp);
                jy += incy;
            }
        } else {
            for j:0..n-1 {
                temp = zero;
                ix = kx;
                if noconj {
                    for i:0..m-1 {
                        temp += get(a,i,j)*get(x,ix);
                        ix += incx;
                    }
                } else {
                    for i:0..m-1 {
                        temp += conjugate(get(a,i,j))*get(x,ix);
                        ix += incx;
                    }
                }
                set(y,iy, get(y,iy) + alpha*temp);
                jy += incy;
            }
        }
    }
}