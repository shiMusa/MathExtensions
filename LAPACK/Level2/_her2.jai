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
        kx = ifx incx > 0 then 0 else (n-1)*(-incx);
        ky = ifx incy > 0 then 0 else (n-1)*(-incy);
        jx = kx;
        jy = ky;
    }

    // Start the operations. In this version the elements of A are
    // accessed sequentially with one pass through the triangular part
    // of a.

    if char_eq(uplo, #char "u") {
        // Form  A  when A is stored in the upper triangle.
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                if get(x,j) != zero || get(y,j) != zero {
                    temp1 = alpha*conjugate(get(y,j));
                    temp2 = conjugate(alpha*get(x,j));
                    for i:0..j-1 {
                        set(a,i,j, get(a,i,j) + get(x,i)*temp1 + get(y,i)*temp2);
                    }
                    set(a,j,j, real_nc(get(a,j,j)) + real_nc( get(x,j) * temp1 + get(y,j) * temp2 ) );
                } else {
                    #if nonreal set(a,j,j, real_nc(get(a,j,j)));
                }
            }
        } else {
            for j:0..n-1 {
                if get(x,jx) != zero || get(y,iy) != zero {
                    temp1 = alpha*conjugate(get(y,iy));
                    temp2 = conjugate(alpha*get(x,jx));
                    ix = kx;
                    iy = ky;
                    for i:0..j-1 {
                        set(a,i,j, get(a,i,j) + get(x,ix)*temp1 + get(y,iy)*temp2);
                        ix += incx;
                        iy += incy;
                    }
                    set(a,j,j, real_nc(get(a,j,j)) + 
                        real_nc(get(x,jx)*temp1 + get(y,iy)*temp2)
                    );
                } else {
                    #if nonreal set(a,j,j, real_nc(get(a,j,j)));
                }
            }
        }
    } else {
        // Form  A  when A is stored in the lower triangle.
        if incx == 1 && incy == 1 {
            for j:0..n-1 {
                if get(x,j) != zero || get(y,iy) != zero {
                    temp1 = alpha*conjugate(get(y,j));
                    temp2 = conjugate(alpha*get(x,j));
                    set(a,j,j, real_nc(get(a,j,j)) +
                        real_nc(get(x,j)*temp1 + get(y,j)*temp2)
                    );
                    for i:j..n-1 {
                        set(a,i,j, get(a,i,j) + get(x,i)*temp1 + get(y,i)*temp2);
                    }
                } else {
                    #if nonreal set(a,j,j, real_nc(get(a,j,j)));
                }
            }
        } else {
            for j:0..n-1 {
                if get(x,jx) != zero || get(y,iy) != zero {
                    temp1 = alpha*conjugate(get(y,iy));
                    temp2 = conjugate(alpha*get(x,jx));
                    set(a,j,j, real_nc(get(a,j,j)) + 
                        real_nc(get(x,jx)*temp1 + get(y,iy)*temp2)
                    );
                    ix = jx;
                    iy = jy;
                    for i:j..n-1 {
                        ix += incx;
                        iy += incy;
                        set(a,i,j, get(a,i,j) + get(x,ix)*temp1 + get(y,iy)*temp2);
                    }
                } else {
                    #if nonreal set(a,j,j, real_nc(get(a,j,j)));
                }
                jx += incx;
                jy += incy;
            }
        }
    }
}