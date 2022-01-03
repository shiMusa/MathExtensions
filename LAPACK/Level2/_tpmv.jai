// ? LAPACK *tpmv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gafb7324ce48931e58b392dd6eef9a286c.html#gafb7324ce48931e58b392dd6eef9a286c

// ? tpmv = triangular packed matrix * vector

// ? original comment
// CTPMV  performs one of the matrix-vector operations
//
//    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit,
// upper or lower triangular matrix, supplied in packed form.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

stpmv :: _tpmv;
dtpmv :: _tpmv;
ctpmv :: _tpmv;
ztpmv :: _tpmv;

_tpmv :: (
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
    ix, jx, k, kk : int = ---;

    // Start the operations. In this version the elements of AP are
    // accessed sequentially with one pass through AP.

    if char_eq(trans, #char "n") {
        // Form  x := A*x.
        if char_eq(uplo, #char "u") {
            kk = 1; // ! 1 offset?
            if incx == 1 {
                for j:0..n-1 {
                    if get(x,j) != zero {
                        temp = get(x,j);
                        k = kk;
                        for i:0..j-1 {
                            set(x,i, get(x,i) + temp * get(ap,k));
                            k += 1;
                        }
                        if nounit then set(x,j, get(x,j) * get(ap,kk+j-1));
                    }
                    kk += j;
                }
            } else {
                jx = kx;
                for j:0..n-1 {
                    if get(x,jx) != zero {
                        temp = get(x,jx);
                        ix = kx;
                        for k:kk..kk+j-2 {
                            set(x,ix, get(x,ix) + temp * get(ap,k));
                            ix += incx;
                        }
                        if nounit then set(x,jx, get(x,jx)*get(ap,kk+j-1));
                    }
                    jx += incx;
                    kk += j;
                }
            }
        } else {
            kk = (n*(n+1))/2;
            if incx == 1 {
                for <j:n-1..0 {
                    if get(x,j) != zero {
                        temp = get(x,j);
                        k = kk;
                        for <i:n-1..j+1 {
                            set(x,i, get(x,i) + temp*get(ap,k));
                            k -= 1;
                        }
                        if nounit then set(x,j, get(x,j) * get(ap,kk-n+j));
                    }
                    kk -= (n-1)*incx;
                }
            } else {
                kx += (n-1)*incx;
                jx = kx;
                for <j:n-1..0 {
                    if get(x,jx) != zero {
                        temp = get(x,jx);
                        ix = kx;
                        for k:kk..kk-(n-(j+1)) {
                            set(x,ix, get(x,ix) + temp*get(ap,k));
                            ix -= incx;
                        }
                        if nounit then set(x,jx, get(x,jx)*get(ap,kk-n+j));
                    }
                    jx -= incx;
                    kk -= n-j+1;
                }
            }
        }
    } else {
        // Form  x := A**T*x  or  x := A**H*x.
        if char_eq(uplo, #char "u") {
            kk = (n*(n+1))/2;
            if incx == 1 {
                for <j:n-1..0 {
                    temp = get(x,j);
                    k = kk - 1;
                    if noconj {
                        if nounit then temp *= get(ap,kk);
                        for <i:j-1..0 {
                            temp += get(ap,k)*get(x,i);
                            k -= 1;
                        }
                    } else {
                        if nounit then temp *= conjugate(get(ap,kk));
                        for <i:j-1..0 {
                            temp += conjugate(get(ap,k))*get(x,i);
                            k -= 1;
                        }
                    }
                    set(x,j, temp);
                    kk = kk - j;
                }
            } else {
                jx = kx + (n-1)*incx;
                for <j:n-1..0 {
                    temp = get(x,jx);
                    ix = jx;
                    if noconj {
                        if nounit then temp *= get(ap,kk);
                        for <k:kk-1..kk-j+1 {
                            ix -= incx;
                            temp += get(ap,k)*get(x,ix);
                        }
                    } else {
                        if nounit then temp *= conjugate(get(ap,kk));
                        for <k:kk-1..kk-j+1 {
                            ix -= incx;
                            temp += conjugate(get(ap,k))*get(x,ix);
                        }
                    }
                    set(x,jx, temp);
                    jx -= incx;
                    kk -= j;
                }
            }
        } else {
            kk = 1;
            if incx == 1 {
                for j:0..n-1 {
                    temp = get(x,j);
                    k = kk + 1;
                    if noconj {
                        if nounit then temp *= get(ap,kk);
                        for i:j+1..n-1 {
                            temp += get(ap,k)*get(x,i);
                            k += 1;
                        }
                    } else {
                        if nounit then temp *= conjugate(get(ap,kk));
                        for i:j+1..n-1 {
                            temp += conjugate(get(ap,k))*get(x,i);
                            k += 1;
                        }
                    }
                    set(x,j, temp);
                    kk += n-j+1;
                }
            } else {
                jx = kx;
                for j:0..n-1 {
                    temp = get(x,jx);
                    ix = jx;
                    if noconj {
                        if nounit then temp *= get(ap,kk);
                        for k:kk+1..kk+n-j {
                            ix += incx;
                            temp += get(ap,k)*get(x,ix);
                        }
                    } else {
                        if nounit then temp *= conjugate(get(ap,kk));
                        for k:kk+1..kk+n-j {
                            ix += incx;
                            temp += conjugate(get(ap,k))*get(x,ix);
                        }
                    }
                    set(x,jx, temp);
                    jx += incx;
                    kk += n-j+1;
                }
            }
        }
    }
}