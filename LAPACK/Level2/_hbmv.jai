// ? LAPACK *hbmv() or *sbmv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga426c198e1c8863d9ecbdbe16efa0f6cf.html#ga426c198e1c8863d9ecbdbe16efa0f6cf

// ? original comment
// CHBMV  performs the matrix-vector  operation
//
//    y := alpha*A*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n hermitian band matrix, with k super-diagonals.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

ssbmv :: _hbmv;
dsbmv :: _hbmv;
chbmv :: _hbmv;
zhbmv :: _hbmv;

_hbmv :: (
    $$uplo: u8,
    $$n: int,
    $$k: int,
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
        assert(k >= 0);
        assert(lda >= k+1);
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible
    if n == 0 || (alpha == _zero && beta == _one) then return;

    // set up the start points in x and y
    kx := ifx incx > 0 then 1 else 1 - (n-1)*incx;
    ky := ifx incy > 0 then 1 else 1 - (n-1)*incy;

    temp1, temp2 : T = ---;
    ix, iy, jx, jy, kplus1, l : int = ---;

    // start the operations. in this version the elements of the array a
    // are accessed sequentially with _one pass through a.

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
        // Form  y  when upper triangle of A is stored.
        kplus1 = k + 1;
        if incx == 1 && incy == 1 {
            for j:1..n {
                temp1 = alpha*fget(x, j);
                temp2 = _zero;
                l = kplus1 - j;
                for i:max(1,j-k)..j-1 {
                    fset(y, i, fget(y.*, i) + temp1*fget(a, l+i, j));
                    temp2 += conjugate(fget(a, l+i, j))*fget(x, i);
                }
                fset(y, j, fget(y.*, j) + temp1*real(fget(a, kplus1, j)) + alpha*temp2);
            }
        } else {
            jx = kx;
            jy = ky;
            for j:1..n {
                temp1 = alpha * fget(x, jx);
                temp2 = _zero;
                ix = kx;
                iy = ky;
                l = kplus1 - j;
                for i:max(1,j-k)..j-1 {
                    fset(y, iy, fget(y.*, iy) + temp1*fget(a, l+i, j));
                    temp2 += conjugate(fget(a, l+i, j))*fget(x, ix);
                    ix += incx;
                    iy += incy;
                }
                fset(y, jy, fget(y.*, jy) + temp1*real(fget(a, kplus1, j)) + alpha*temp2);
                jx += incx;
                jy += incy;
                if j > k {
                    kx += incx;
                    ky += incy;
                }
            }
        }
    } else {
        // Form  y  when lower triangle of A is stored.
        if incx == 1 && incy == 1 {
            for j:1..n {
                temp1 = alpha*fget(x, j);
                temp2 = _zero;
                fset(y, j, fget(y.*, j) + temp1*real(fget(a, 1, j)));
                l = 1 - j;
                for i:j+1..min(n,j+k) {
                    fset(y, i, fget(y.*, i) + temp1*fget(a, l+i, j));
                    temp2 += conjugate(fget(a, l+i, j))*fget(x, i);
                }
                fset(y, j, fget(y.*, j) + alpha*temp2);
            }
        } else {
            jx = kx;
            jy = ky;
            for j:1..n {
                temp1 = alpha*fget(x, jx);
                temp2 = _zero;
                fset(y, jy, fget(y.*, jy) + temp1*real(fget(a, 1, j)));
                l = 1 - j;
                ix = jx;
                iy = jy;
                for i:j+1..min(n,j+k) {
                    ix += incx;
                    iy += incy;
                    fset(y, iy, fget(y.*, iy) + temp1*fget(a, l+i, j));
                    temp2 += conjugate(fget(a, l+i, j))*fget(x, ix);
                }
                fset(y, jy, fget(y.*, jy) + alpha*temp2);
                jx += incx;
                jy += incy;
            }
        }
    }
}
