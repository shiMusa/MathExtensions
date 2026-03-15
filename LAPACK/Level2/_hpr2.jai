// ? LAPACK *hpr2() and *spr2()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga537b8bcd05cbe626ebce0fbc8d66a2d6.html#ga537b8bcd05cbe626ebce0fbc8d66a2d6

// ? hpr2 = hermitian packed matrix rank 2

// ? original comment
// CHPR2  performs the hermitian rank 2 operation
//
//    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
//
// where alpha is a scalar, x and y are n element vectors and A is an
// n by n hermitian matrix, supplied in packed form.

// Level 2 Blas routine.
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.


sspr2 :: _hpr2;
dspr2 :: _hpr2;
chpr2 :: _hpr2;
zhpr2 :: _hpr2;


_hpr2 :: (
    $$uplo: u8,
    $$n: int,
    alpha: $T,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    ap: *$A/ArrayType1D
) {
    zero :: #run zero(T);
    one  :: #run one(T);
    nonreal :: #run is_non_real_number(T);
    
    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == X.data_type);
        }
        assert(char_eq(uplo, #char "u", #char "l"));
        assert(n >= 0);
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible
    if n == 0 || alpha == zero then return;

    temp1, temp2 : T = ---;
    ix, iy, jx, jy, k, kk, kx, ky : int = ---;

    // Set up the start points in x and y if the increments are not both
    // unity.

    if incx != 1 || incy != 1 {
        kx = ifx incx > 0 then 1 else 1 - (n-1)*incx;
        ky = ifx incy > 0 then 1 else 1 - (n-1)*incy;
        jx = kx;
        jy = ky;
    }

    // Start the operations. In this version the elements of the array AP
    // are accessed sequentially with one pass through AP.

    kk = 1;
    if char_eq(uplo, #char "u") {
        // Form  A  when upper triangle is stored in AP.
        if incx == 1 && incy == 1 {
            for j:1..n {
                if fget(x, j) != zero || fget(y, j) != zero {
                    temp1 = alpha*conjugate(fget(y, j));
                    temp2 = conjugate(alpha*fget(x, j));
                    k = kk;
                    for i:1..j-1 {
                        fset(ap, k, fget(ap.*, k) + fget(x, i)*temp1 + fget(y, i)*temp2);
                        k += 1;
                    }
                    fset(ap, kk+j-1, real_nc(fget(ap.*, kk+j-1)) + real_nc(fget(x, j)*temp1 + fget(y, j)*temp2));
                } else {
                    #if nonreal fset(ap, kk+j-1, real_nc(fget(ap.*, kk+j-1)));
                }
                kk += j;
            }
        } else {
            for j:1..n {
                if fget(x, jx) != zero || fget(y, jy) != zero {
                    temp1 = alpha*conjugate(fget(y, jy));
                    temp2 = conjugate(alpha*fget(x, jx));
                    ix = kx;
                    iy = ky;
                    for k:kk..kk+j-2 {
                        fset(ap, k, fget(ap.*, k) + fget(x, ix)*temp1 + fget(y, iy)*temp2);
                        ix += incx;
                        iy += incy;
                    }
                    fset(ap, kk+j-1, real_nc(fget(ap.*, kk+j-1)) + real_nc(fget(x, jx)*temp1 + fget(y, jy)*temp2));
                } else {
                    #if nonreal fset(ap, kk+j-1, real_nc(fget(ap.*, kk+j-1)));
                }
                jx += incx;
                jy += incy;
                kk += j;
            }
        }
    } else {
        // Form  A  when lower triangle is stored in AP.
        if incx == 1 && incy == 1 {
            for j:1..n {
                if fget(x, j) != zero || fget(y, j) != zero {
                    temp1 = alpha*conjugate(fget(y, j));
                    temp2 = conjugate(alpha*fget(x, j));
                    fset(ap, kk, real_nc(fget(ap.*, kk)) + real_nc(fget(x, j)*temp1 + fget(y, j)*temp2));
                    k = kk + 1;
                    for i:j+1..n {
                        fset(ap, k, fget(ap.*, k) + fget(x, i)*temp1 + fget(y, i)*temp2);
                        k += 1;
                    }
                } else {
                    #if nonreal fset(ap, kk, real_nc(fget(ap.*, kk)));
                }
                kk += n-j+1;
            }
        } else {
            for j:1..n {
                if fget(x, jx) != zero || fget(y, jy) != zero {
                    temp1 = alpha*conjugate(fget(y, jy));
                    temp2 = conjugate(alpha*fget(x, jx));
                    fset(ap, kk, real_nc(fget(ap.*, kk)) + real_nc(fget(x, jx)*temp1 + fget(y, jy)*temp2));
                    ix = jx;
                    iy = jy;
                    for k:kk+1..kk+n-j {
                        ix += incx;
                        iy += incy;
                        fset(ap, k, fget(ap.*, k) + fget(x, ix)*temp1 + fget(y, iy)*temp2);
                    }
                } else {
                    #if nonreal fset(ap, kk, real_nc(fget(ap.*, kk)));
                }
                jx += incx;
                jy += incy;
                kk += n-j+1;
            }
        }
    }
}
