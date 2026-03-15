// ? LAPACK *hpmv() and *spmv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gad1dc14308bf0c1c8d9ba6ee068ac4d60.html#gad1dc14308bf0c1c8d9ba6ee068ac4d60

// ? hpmv = hermitian packed matrix * vector

// ? original comment
// CHPMV  performs the matrix-vector operation
//
//    y := alpha*A*x + beta*y,
//
// where alpha and beta are scalars, x and y are n element vectors and
// A is an n by n hermitian matrix, supplied in packed form.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

sspmv :: _hpmv; // conjugate has no effect on real numbers
dspmv :: _hpmv; // conjugate has no effect on real numbers
chpmv :: _hpmv;
zhpmv :: _hpmv;

_hpmv :: (
    $$uplo: u8,
    $$n: int,
    alpha: $T,
    ap: $A/ArrayType1D,
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
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible
    if n == 0 || (alpha == _zero && beta == _one) then return;

    // Set up the start points in x and y.
    kx := ifx incx > 0 then 1 else 1 - (n-1)*incx;
    ky := ifx incy > 0 then 1 else 1 - (n-1)*incy;

    temp1, temp2 : T = ---;
    ix, iy, jx, jy, k, kk : int = ---;

    // Start the operations. In this version the elements of the array AP
    // are accessed sequentially with _one pass through AP.

    // First form  y := beta*y

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
    
    kk = 1;
    if char_eq(uplo, #char "u") {
        // Form  y  when AP contains the upper triangle.
        if incx == 1 && incy == 1 {
            for j:1..n {
                temp1 = alpha*fget(x, j);
                temp2 = _zero;
                k = kk;
                for i:1..j-1 {
                    fset(y, i, fget(y.*, i) + temp1*fget(ap, k));
                    temp2 += conjugate(fget(ap, k))*fget(x, i);
                    k += 1;
                }
                fset(y, j, fget(y.*, j) + temp1*real_nc(fget(ap, kk+j-1)) + alpha*temp2);
                kk += j;
            }
        } else {
            jx = kx;
            jy = ky;
            for j:1..n {
                temp1 = alpha*fget(x, jx);
                temp2 = _zero;
                ix = kx;
                iy = ky;
                for k:kk..kk+j-2 {
                    fset(y, iy, fget(y.*, iy) + temp1*fget(ap, k));
                    temp2 += conjugate(fget(ap, k))*fget(x, ix);
                    ix += incx;
                    iy += incy;
                }
                fset(y, jy, fget(y.*, jy) + temp1*real_nc(fget(ap, kk+j-1)) + alpha*temp2);
                jx += incx;
                jy += incy;
                kk += j;
            }
        }
    } else {
        // Form  y  when AP containse the lower triangle.
        if incx == 1 && incy == 1 {
            for j:1..n {
                temp1 = alpha*fget(x, j);
                temp2 = _zero;
                fset(y, j, fget(y.*, j) + temp1*real_nc(fget(ap, kk)));
                k = kk + 1;
                for i:j+1..n {
                    fset(y, i, fget(y.*, i) + temp1*fget(ap, k));
                    temp2 += conjugate(fget(ap, k))*fget(x, i);
                    k += 1;
                }
                fset(y, j, fget(y.*, j) + alpha * temp2);
                kk += n-j+1;
            }
        } else {
            jx = kx;
            jy = ky;
            for j:1..n {
                temp1 = alpha*fget(x, jx);
                temp2 = _zero;
                fset(y, jy, fget(y.*, jy) + temp1*real_nc(fget(ap, kk)));
                ix = jx;
                iy = jy;
                for k:kk+1..kk+n-j {
                    ix += incx;
                    iy += incy;
                    fset(y, iy, fget(y.*, iy) + temp1*fget(ap, k));
                    temp2 += conjugate(fget(ap, k))*fget(x, ix);
                }
                fset(y, jy, fget(y.*, jy) + alpha*temp2);
                jx += incx;
                jy += incy;
                kk += n-j+1;
            }
        }
    }
}
