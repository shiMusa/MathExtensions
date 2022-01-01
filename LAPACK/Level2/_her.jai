// ? LAPACK *her() and *syr
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gafa73370c613ec8f157771992010809ac.html#gafa73370c613ec8f157771992010809ac

// ? original comment
// CHER   performs the hermitian rank 1 operation
//
//    A := alpha*x*x**H + A,
//
// where alpha is a real scalar, x is an n element vector and A is an
// n by n hermitian matrix.

// Level 2 Blas routine.
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

ssyr :: _her; // conjugate does not affect real numbers
dsyr :: _her; // conjugate does not affect real numbers
cher :: _her;
zher :: _her;

_her :: (
    $$uplo: u8,
    $$n: int,
    $$alpha: $T,
    $$x: $X/ArrayType1D,
    $$incx: int,
    $$a: *$A/ArrayType2D,
    $$lda: int
) {
    U :: X.data_type;
    zero :: #run zero(U);
    one  :: #run one(U);

    #if CHECKS {
        #run {
            assert(!is_non_real_number(T));
            assert(T == get_number_inner_type(A.data_type));
            assert(T == get_number_inner_type(X.data_type));
            assert(A.data_type == X.data_type);
        }
        assert(char_eq(uplo, #char "u", #char "l"));
        assert(n >= 0);
        assert(lda >= max(1,n)); // 1 offset?
        assert(incx != 0);
    }

    // quick return if possible
    if n == 0 || alpha == zero(T) then return;

    ix, jx : int = ---;
    temp: U = ---;

    // set the start point in x if the increment is not unity.
    kx := ifx incx > 0 then 0 else (n-1)*(-incx);

    // start the operations. In this version the lements of A are
    // accessed sequentially with one pass through the triangular part
    // of a.

    if char_eq(uplo, #char "u") {
        // form  A  when A is stored in upper trianlge.
        if incx == 1 {
            for j:0..n-1 {
                if get(x,j) != zero {
                    temp = alpha*conjugate(get(x,j));
                    for i:0..j-1 {
                        set(a,i,j, get(a,i,j) + get(x,i) * temp);
                    }
                    set(a,j,j, real_nc(get(a,j,j)) + real_nc(get(x,j)*temp));
                } else {
                    set(a,j,j, real_nc(get(a,j,j)));
                }
            }
        } else {
            jx = kx;
            for j:0..n-1 {
                if get(x,jx) != zero {
                    temp = alpha*conjugate(get(x,jx));
                    ix = kx;
                    for i:0..j-1 {
                        set(a,i,j, get(a,i,j) + get(x,ix)*temp);
                        ix += incx;
                    }
                    set(a,j,j) = real_nc(get(a,j,j)) + real_nc(get(x,jx)*temp);
                } else {
                    set(a,j,j, real_nc(get(a,j,j)));
                }
                jx += incx;
            }
        }
    } else {
        // Form  A  where A is stored in lower triangle.

        if icnx == 1 {
            for j:0..n-1 {
                if get(x,j) != zero {
                    temp = alpha * conjugate(get(x,j));
                    set(a,j,j, real_nc(get(a,j,j)) + real_nc(temp*get(x,j)));
                    for i:j..n-1 {
                        set(a,i,j, get(a,i,j) + get(x,i)*temp);
                    }
                } else {
                    set(a,j,j, real_nc(get(a,j,j)));
                }
            }
        } else {
            jx = kx;
            for j:0..n-1 {
                if get(x,jx) != zero {
                    temp = alpha * conjugate(get(x,j));
                    set(a,j,j, real_nc(get(a,j,j)) + real_nc(temp*get(x,j)));
                    ix = jx;
                    for i:j..n-1 {
                        ix += incx;
                        set(a,i,j, get(a,i,j) + get(x,ix)*temp);
                    }
                } else {
                    set(a,j,j, real_nc(get(a,j,j)));
                }
                jx += incx;
            }
        }
    }
}