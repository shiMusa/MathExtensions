// ? LAPACK *her2() & *syr2()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gaa59d93fbbd8d0b1be4a51634cb437cc1.html#gaa59d93fbbd8d0b1be4a51634cb437cc1

// ? original comment
// CHER2  performs the hermitian rank 2 operation
//
//    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
//
// where alpha is a scalar, x and y are n element vectors and A is an n
// by n hermitian matrix.

// Level 2 Blas routine.
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

dsyr2 :: _her2; // conjugate has no effect on real numbers
ssyr2 :: _her2; // conjugate has no effect on real numbers
cher2 :: _her2;
zher2 :: _her2;

_her2 :: (
    $$uplo: u8,
    $$n: int,
    alpha: $T,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    a: *$A/ArrayType2D,
    $$lda: int
) {
    zero :: #run zero(T);
    one  :: #run one(T);
    nonreal :: #run is_non_real_number(T);

    #if CHECKS {
        #run {
            assert(T == X.data_type);
            assert(T == Y.data_type);
            assert(T == A.data_type);
        }
        assert(char_eq(uplo, #char "u", #char "l"));
        assert(n >= 0);
        assert(lda >= max(1,n)); // 1 offset?
        assert(incx != 0);
        assert(incy != 0);
    }

    // quick return if possible
    if n == 0 || alpha == zero then return;

    temp1, temp2 : T = ---;
    ix, iy, jx, jy, kx, ky : int = ---;

    // Set up the start points in x and y if the increments are not both unity.
    if incx != 1 || incy != 1 {
        kx = ifx incx > 0 then 1 else 1 - (n-1)*incx;
        ky = ifx incy > 0 then 1 else 1 - (n-1)*incy;
        jx = kx;
        jy = ky;
    }

    // Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through the triangular part
    // of a.

    if char_eq(uplo, #char "u") {
        // Form  A  when A is stored in the upper triangle.
        if incx == 1 && incy == 1 {
            for j:1..n {
                if fget(x, j) != zero || fget(y, j) != zero {
                    temp1 = alpha*conjugate(fget(y, j));
                    temp2 = conjugate(alpha*fget(x, j));
                    for i:1..j-1 {
                        fset(a, i, j, fget(a.*, i, j) + fget(x, i)*temp1 + fget(y, i)*temp2);
                    }
                    fset(a, j, j, real_nc(fget(a.*, j, j)) + real_nc(fget(x, j) * temp1 + fget(y, j) * temp2));
                } else {
                    #if nonreal fset(a, j, j, real_nc(fget(a.*, j, j)));
                }
            }
        } else {
            for j:1..n {
                if fget(x, jx) != zero || fget(y, jy) != zero {
                    temp1 = alpha*conjugate(fget(y, jy));
                    temp2 = conjugate(alpha*fget(x, jx));
                    ix = kx;
                    iy = ky;
                    for i:1..j-1 {
                        fset(a, i, j, fget(a.*, i, j) + fget(x, ix)*temp1 + fget(y, iy)*temp2);
                        ix += incx;
                        iy += incy;
                    }
                    fset(a, j, j, real_nc(fget(a.*, j, j)) + real_nc(fget(x, jx)*temp1 + fget(y, jy)*temp2));
                } else {
                    #if nonreal fset(a, j, j, real_nc(fget(a.*, j, j)));
                }
                jx += incx;
                jy += incy;
            }
        }
    } else {
        // Form  A  when A is stored in the lower triangle.
        if incx == 1 && incy == 1 {
            for j:1..n {
                if fget(x, j) != zero || fget(y, j) != zero {
                    temp1 = alpha*conjugate(fget(y, j));
                    temp2 = conjugate(alpha*fget(x, j));
                    fset(a, j, j, real_nc(fget(a.*, j, j)) + real_nc(fget(x, j)*temp1 + fget(y, j)*temp2));
                    for i:j+1..n {
                        fset(a, i, j, fget(a.*, i, j) + fget(x, i)*temp1 + fget(y, i)*temp2);
                    }
                } else {
                    #if nonreal fset(a, j, j, real_nc(fget(a.*, j, j)));
                }
            }
        } else {
            for j:1..n {
                if fget(x, jx) != zero || fget(y, jy) != zero {
                    temp1 = alpha*conjugate(fget(y, jy));
                    temp2 = conjugate(alpha*fget(x, jx));
                    fset(a, j, j, real_nc(fget(a.*, j, j)) + real_nc(fget(x, jx)*temp1 + fget(y, jy)*temp2));
                    ix = jx;
                    iy = jy;
                    for i:j+1..n {
                        ix += incx;
                        iy += incy;
                        fset(a, i, j, fget(a.*, i, j) + fget(x, ix)*temp1 + fget(y, iy)*temp2);
                    }
                } else {
                    #if nonreal fset(a, j, j, real_nc(fget(a.*, j, j)));
                }
                jx += incx;
                jy += incy;
            }
        }
    }
}
