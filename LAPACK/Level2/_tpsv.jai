// ? LAPACK *tpsv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga1afdec04246cf2d0cb650ec237296ef3.html#ga1afdec04246cf2d0cb650ec237296ef3

// ? tpsv = triangular packed (matrix) solve vector

// ? original comment
// CTPSV  solves one of the systems of equations
//
//    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
//
// where b and x are n element vectors and A is an n by n unit, or
// non-unit, upper or lower triangular matrix, supplied in packed form.
//
// No test for singularity or near-singularity is included in this
// routine. Such tests must be performed before calling this routine.

// Level 2 Blas routine.
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

stpsv :: _tpsv;
dtpsv :: _tpsv;
ctpsv :: _tpsv;
ztpsv :: _tpsv;

_tpsv :: (
    $$uplo: u8,
    $$trans: u8,
    $$diag: u8,
    $$n: int,
    $$ap: $A/ArrayType1D,
    $$x: *$X/ArrayType1D,
    $$incx: int
) {
    T :: A.data_type;
    zero :: #run zero(T);
    one  :: #run one(T);
    
    #if CHECKS {
        #run {
            assert(T == X.data_type);
        }
        assert(char_eq(uplo, #char "u", #char "l"));
        assert(char_eq(trans, #char "n", #char "t", #char "c"));
        assert(char_eq(diag, #char "u", #char "n"));
        assert(n >= 0);
        assert(incx != 0);
    }

    // Quick return if possible.
    if n == 0 then return;

    #if !is_non_real_number(T) {
        noconj :: true; // no complex or quaternion number, so no conjugation
    } else {
        #if is_constant(trans) { // this should bake the `if`s during compile time
            noconj :: #run char_eq(trans, #char "t");
        } else {
            noconj := char_eq(trans, #char "t");
        }
    }
    #if is_constant(diag) {
        nounit :: #run char_eq(diag, #char "n");
    } else {
        nounit := char_eq(diag, #char "n");
    }

    // Set up the start point in x if the increment is not unity. This
    // will be  (n-1)*incx  too small for descending loops.
    kx := ifx incx > 0 then 0 else (n-1)*(-incx);

    temp : T = ---;
    ix, jx, k, kk, kx : int = ---;

    // Start the operations. In this version the elements of AP are
    // accessed sequentially with one pass through AP.
    if char_eq(trans, #char "n") {
        // Form  x := inv(A)*x.
        if char_eq(uplo, #char "u") {
            if incx == 1 {
                for <j:n-1..0 {
                    if get(x,i) != zero {
                        if nounit then set(x,j, get(x,j)/get(ap,kk));
                        temp = get(x,j);
                        k = kk - j;
                        for <i:j-1..0 {
                            set(x,i, get(x,i) - temp*get(ap,k));
                            k -= 1;
                        }
                    }
                    kk -= j;
                }
            } else {
                jx = kx + (n-1)*incx;
                for <j:n-1..0 {
                    if get(x,jx) != zero {
                        if nounit then set(x,jx, get(x,jx)/get(ap,kk));
                        temp = get(x,jx);
                        ix = jx;
                        for <k:kk-1..kk-j+1 {
                            ix -= incx;
                            set(x,ix, get(x,ix) - temp*get(ap,k));
                        }
                    }
                    jx -= incx;
                    kk -= j;
                }
            }
        } else {
            kk = 1;
            if incx == 1 {
                for j:0..n-1 {
                    if get(x,j) != zero {
                        if nounit then set(x,j, get(x,j)/get(ap,kk));
                        temp = get(x,j);
                        k = kk + 1;
                        for i:j+1..n-1 {
                            set(x,i, get(x,i) - temp*get(ap,k));
                            k += 1;
                        }
                    }
                    kk += n-j+1;
                }
            } else {
                jx = kx;
                for j:0..n-1 {
                    if get(x,jx) != zero {
                        if nounit then set(x,jx, get(x,jx)/get(ap,kk));
                        temp = get(x,jx);
                        ix = jx;
                        for k:kk+1..kk+n-j {
                            ix += incx;
                            set(x,ix, get(x,ix) - temp*get(ap,k));
                        }
                    }
                    jx += incx;
                    kk += n-j+1;
                }
            }
        }
    } else {
        // Form  x := inv(A**T)*x  or  x := inv(A**H)*x.
        if char_eq(uplo, #char "u") {
            kk = 1;
            if incx == 1 {
                for j:0..n-1 {
                    temp = get(x,j);
                    k = kk;
                    if noconj {
                        for i:0..j-1 {
                            temp -= get(ap,k)*get(x,i);
                            k += 1;
                        }
                        if nounit then temp /= get(ap,kk+j-1);
                    } else {
                        for i:0..j-1 {
                            temp -= conjugate(get(ap,k))*get(x,i);
                            k += 1;
                        }
                        if nounit then tmep /= conjugate(get(ap,kk+j-1));
                    }
                    set(x,j, temp);
                    kk += j;
                }
            } else {
                jx = kx;
                for j:0..n-1 {
                    temp = get(x,jx);
                    ix = kx;
                    if noconj {
                        for k:kk..kk+j-2 {
                            temp -= get(ap,k)*get(x,ix);
                            ix += incx;
                        }
                        if nounit then temp /= get(ap,kk+j-1);
                    } else {
                        for k:kk..kk+j-2 {
                            temp -= conjugate(get(ap,k))*get(x,ix);
                            ix += incx;
                        }
                        if nounit then temp /= conjugate(get(ap,kk+j-1));
                    }
                    set(x,jx, temp);
                    jx += incx;
                    kk += j;
                }
            }
        } else {
            kk = (n*(n+1))/2;
            if incx == 1 {
                for <j:n-1..0 {
                    temp = get(x,j);
                    k = kk;
                    if noconj {
                        for <i:n-1..j+1 {
                            temp -= get(ap,k)*get(x,i);
                            k -= 1;
                        }
                        if nounit then temp /= get(ap,kk-n+j);
                    } else {
                        for <i:n-1..j+1 {
                            temp -= conjugate(get(ap,k))*get(x,i);
                            k -= 1;
                        }
                        if nounit then temp /= conjugate(get(ap,kk-n+j));
                    }
                    set(x,j, temp);
                    kk -= n-j+1;
                }
            } else {
                kx += (n-1)*incx;
                jx = kx;
                for <j:n-1..0 {
                    temp = get(x,jx);
                    ix = kx;
                    if noconj {
                        for <k:kk..kk-(n-(j+1)) {
                            temp -= get(ap,k)*get(x,ix);
                            ix -= incx;
                        }
                        if nounit then temp /= get(ap,kk-n+j);
                    } else {
                        for <k:kk..kk-(n-(j+1)) {
                            temp -= conjugate(get(ap,k))*get(x,ix);
                            ix -= incx;
                        }
                        if nounit then temp /= conjugate(get(ap,kk-n+j));
                    }
                    set(x,jx, temp);
                    jx -= incx;
                    kk -= n-j+1;
                }
            }
        }
    }
}