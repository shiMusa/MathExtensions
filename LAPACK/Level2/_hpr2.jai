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
    $$alpha: $T,
    $$x: $X/ArrayType1D,
    $$incx: int,
    $$y: $Y/ArrayType1D,
    $$incy: int,
    $$ap: *$A/ArrayType1D
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

    temp1, temp1 : T = ---;
    ix, iy, jx, jy, k, kk, kx, ky : int = ---;

    // Set up the start points in x and y if the increments are not both
    // unity.

    if incx != 1 || incy != 1 {
        kx = ifx incx > 0 then 0 else (n-1)*(-incx);
        ky = ifx incy > 0 then 0 else (n-1)*(-incy);
        jx = kx;
        jy = ky;
    }

    // Start the operations. In this version the elements of the array AP
    // are accessed sequentially with one pass through AP.

    kk = 1;
    if char_eq(uplo, #char "u") {
        // Form  A  when upper triangle is stored in AP.
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                if get(x,j) != zero || get(y,j) != zero {
                    temp1 = alpha*conjugate(get(y,j));
                    temp2 = conjugate(alpha*get(x,j));
                    k = kk;
                    for i:0..j-1 {
                        set(ap, k, get(ap,k) + get(x,i)*temp1 + get(y,i)*temp2);
                        k += 1;
                    }
                } else {
                    #if nonreal set(ap, kk+j-1, real_nc(get(ap,kk+j-1)));
                }
                kk += j;
            }
        } else {
            for j:0..n-1 {
                if get(x,jx) != zero || get(y,jy) != zero {
                    temp1 = alpha*conjugate(get(y,jy));
                    temp2 = conjugate(alpha*get(x,jx));
                    ix = kx;
                    iy = ky;
                    for k:kk..kk+j-2 {
                        set(ap,k, get(ap,k) + get(x,ix)*temp1 + get(y,iy)*temp2);
                        ix += incx;
                        iy += incy;
                    }
                    set(ap,kk+j-1, real_nc(get(ap,kk+j-1)) + 
                        real_nc(get(x,jx)*temp1 + get(y,jy)*temp2)
                    );
                } else {
                    #if nonreal set(ap,kk+j-1, real_nc(get(ap,kk+j-1)));
                }
                jx += incx;
                jy += incy;
                kk += j;
            }
        }
    } else {
        // Form  A  when lower triangle is stored in AP.
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                if get(x,j) != zero || get(y,j) != zero {
                    temp1 = alpha*conjugate(get(y,j));
                    temp2 = conjugate(alpha*get(x,j));
                    set(ap,kk, real_nc(get(ap,kk)) +
                        real_nc(get(x,j)*temp1 + get(y,j)*temp2)
                    );
                    k = kk + 1;
                    for i:j+1..n-1 {
                        set(ap,k, get(ap,k) + get(x,i)*temp1 + get(y,i)*temp2);
                        k += 1;
                    }
                } else {
                    #if nonreal set(ap,kk, real_nc(get(ap,kk)));
                }
                kk += n-j+1;
            }
        } else {
            for j:0..n-1 {
                if get(x,jx) != zero || get(y,jy) != zero {
                    temp1 = alpha*conjugate(get(y,jy));
                    temp2 = conjugate(alpha*get(x,jx));
                    set(ap,kk, real_nc(get(ap,kk)) + 
                        real_nc(get(x,jx)*temp1 + get(y,jy)*temp2)
                    );
                    ix = jx;
                    iy = jy;
                    for k:kk+1..kk+n-j {
                        ix += incx;
                        iy += incy;
                        set(ap,k, get(ap,k) + get(x,ix)*temp1 + get(y,iy)*temp2);
                    }
                } else {
                    #if nonreal set(ap,kk, real_nc(get(ap,kk)));
                }
                jx += incx;
                jy += incy;
                kk += n-j+1;
            }
        }
    }
}