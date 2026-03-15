// ? LAPACK *hpr() and *spr()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gac98100eee08124fb7a5d55effbb85a65.html#gac98100eee08124fb7a5d55effbb85a65

// ? hpr = hermitian packed matrix * vector
// ? spr = symmetric packed matrix * vector

// ? original comment
// CHPR    performs the hermitian rank 1 operation
//
//    A := alpha*x*x**H + A,
//
// where alpha is a real scalar, x is an n element vector and A is an
// n by n hermitian matrix, supplied in packed form.

//  Level 2 Blas routine.
//  
//  -- Written on 22-October-1986.
//     Jack Dongarra, Argonne National Lab.
//     Jeremy Du Croz, Nag Central Office.
//     Sven Hammarling, Nag Central Office.
//     Richard Hanson, Sandia National Labs.

sspr :: _hpr;
dspr :: _hpr;
chpr :: _hpr;
zhpr :: _hpr;

_hpr :: (
    $$uplo: u8,
    $$n: int,
    alpha: $T,
    x: $X/ArrayType1D,
    $$incx: int,
    ap: $A/ArrayType1D
) {
    U :: X.data_type;
    zero :: #run zero(U);
    one  :: #run one(U);
    nonreal :: #run is_non_real_number(U);
    
    #if CHECKS {
        #run {
            assert(!is_non_real_number(T));
            assert(T == get_number_inner_type(A.data_type));
            assert(T == get_number_inner_type(X.data_type));
            assert(A.data_type == X.data_type);
        }
        assert(char_eq(uplo, #char "u", #char "l"));
        assert(n >= 0);
        assert(incx != 0);
    }

    // Quick return if possible
    if n == 0 || alpha == zero(T) then return;

    temp : U = ---;
    ix, jx, k, kk, kx : int = ---;

    // Set the start point in x if the increment is not unity
    kx := ifx incx > 0 then 1 else 1 - (n-1)*incx;

    // Start the operations. In this version the elements of the array AP
    // are accessed sequentially with one pass through AP.

    kk = 1;
    if char_eq(uplo, #char "u") {
        // Form  A  when upper triangle is stored in AP.
        if incx == 1 {
            for j:1..n {
                if fget(x, j) != zero {
                    temp = alpha * conjugate(fget(x, j));
                    k = kk;
                    for i:1..j-1 {
                        fset(ap, k, fget(ap, k) + fget(x, i)*temp);
                        k += 1;
                    }
                    fset(ap, kk+j-1, real_nc(fget(ap, kk+j-1)) + real_nc(fget(x, j))*temp);
                } else {
                    #if nonreal fset(ap, kk+j-1, real_nc(fget(ap, kk+j-1)));
                }
                kk += j;
            }
        } else {
            jx = kx;
            for j:1..n {
                if fget(x, jx) != zero {
                    temp = alpha*conjugate(fget(x, jx));
                    ix = kx;
                    for k:kk..kk+j-2 {
                        fset(ap, k, fget(ap, k) + fget(x, ix)*temp);
                        ix += incx;
                    }
                    fset(ap, kk+j-1, real_nc(fget(ap, kk+j-1)) + real_nc(fget(x, jx))*temp);
                } else {
                    #if nonreal fset(ap, kk+j-1, real_nc(fget(ap, kk+j-1)));
                }
                jx += incx;
                kk += j;
            }
        }
    } else {
        // Form  A  when lower triangle is stored in AP.
        if incx == 1 {
            for j:1..n {
                if fget(x, j) != zero {
                    temp = alpha*conjugate(fget(x, j));
                    fset(ap, kk, real_nc(fget(ap, kk)) + real_nc(temp*fget(x, j)));
                    k = kk + 1;
                    for i:j+1..n {
                        fset(ap, k, fget(ap, k) + fget(x, i)*temp);
                        k += 1;
                    }
                } else {
                    #if nonreal fset(ap, kk, real_nc(fget(ap, kk)));
                }
                kk += n-j+1;
            }
        } else {
            jx = kx;
            for j:1..n {
                if fget(x, jx) != zero {
                    temp = alpha*conjugate(fget(x, jx));
                    fset(ap, kk, real_nc(fget(ap, kk)) + real_nc(temp*fget(x, jx)));
                    ix = jx;
                    for k:kk+1..kk+n-j {
                        ix += incx;
                        fset(ap, k, fget(ap, k) + fget(x, ix)*temp);
                    }
                } else {
                    #if nonreal fset(ap, kk, real_nc(fget(ap, kk)));
                }
                jx += incx;
                kk += n-j+1;
            }
        }
    }
}
