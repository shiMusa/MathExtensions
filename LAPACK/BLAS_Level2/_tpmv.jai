// ? LAPACK *tpmv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_gafb7324ce48931e58b392dd6eef9a286c.html#gafb7324ce48931e58b392dd6eef9a286c

// ? tpmv = triangular packed matrix * vector

// ? original comment
// CTPMV  performs _one of the matrix-vector operations
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
    ap: $A/ArrayType1D,
    x: *$X/ArrayType1D,
    $$incx: int
) {
    T :: A.data_type;
    _zero :: #run zero(T);
    _one  :: #run one(T);
    
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

    #if #run !is_non_real_number(T) {
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
    kx := ifx incx > 0 then 1 else 1 - (n-1)*incx;

    temp : T = ---;
    ix, jx, k, kk : int = ---;

    // Start the operations. In this version the elements of AP are
    // accessed sequentially with _one pass through AP.

    if char_eq(trans, #char "n") {
        // Form  x := A*x.
        if char_eq(uplo, #char "u") {
            kk = 1;
            if incx == 1 {
                for j:1..n {
                    if fget(x, j) != _zero {
                        temp = fget(x, j);
                        k = kk;
                        for i:1..j-1 {
                            fset(x, i, fget(x.*, i) + temp * fget(ap, k));
                            k += 1;
                        }
                        if nounit then fset(x, j, fget(x.*, j) * fget(ap, kk+j-1));
                    }
                    kk += j;
                }
            } else {
                jx = kx;
                for j:1..n {
                    if fget(x, jx) != _zero {
                        temp = fget(x, jx);
                        ix = kx;
                        for k:kk..kk+j-2 {
                            fset(x, ix, fget(x.*, ix) + temp * fget(ap, k));
                            ix += incx;
                        }
                        if nounit then fset(x, jx, fget(x.*, jx)*fget(ap, kk+j-1));
                    }
                    jx += incx;
                    kk += j;
                }
            }
        } else {
            kk = (n*(n+1))/2;
            if incx == 1 {
                for <j:1..n {
                    if fget(x, j) != _zero {
                        temp = fget(x, j);
                        k = kk;
                        for <i:j+1..n {
                            fset(x, i, fget(x.*, i) + temp*fget(ap, k));
                            k -= 1;
                        }
                        if nounit then fset(x, j, fget(x.*, j) * fget(ap, kk-n+j));
                    }
                    kk -= n-j+1;
                }
            } else {
                kx += (n-1)*incx;
                jx = kx;
                for <j:1..n {
                    if fget(x, jx) != _zero {
                        temp = fget(x, jx);
                        ix = kx;
                        k = kk;
                        for <i:j+1..n {
                            fset(x, ix, fget(x.*, ix) + temp*fget(ap, k));
                            ix -= incx;
                            k -= 1;
                        }
                        if nounit then fset(x, jx, fget(x.*, jx)*fget(ap, kk-n+j));
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
                for <j:1..n {
                    temp = fget(x, j);
                    k = kk - 1;
                    if noconj {
                        if nounit then temp *= fget(ap, kk);
                        if j > 1 for <i:1..j-1 {
                            temp += fget(ap, k)*fget(x, i);
                            k -= 1;
                        }
                    } else {
                        if nounit then temp *= conjugate(fget(ap, kk));
                        if j > 1 for <i:1..j-1 {
                            temp += conjugate(fget(ap, k))*fget(x, i);
                            k -= 1;
                        }
                    }
                    fset(x, j, temp);
                    kk -= j;
                }
            } else {
                jx = kx + (n-1)*incx;
                for <j:1..n {
                    temp = fget(x, jx);
                    ix = jx;
                    k = kk - 1;
                    if noconj {
                        if nounit then temp *= fget(ap, kk);
                        if j > 1 for <i:1..j-1 {
                            ix -= incx;
                            temp += fget(ap, k)*fget(x, ix);
                            k -= 1;
                        }
                    } else {
                        if nounit then temp *= conjugate(fget(ap, kk));
                        if j > 1 for <i:1..j-1 {
                            ix -= incx;
                            temp += conjugate(fget(ap, k))*fget(x, ix);
                            k -= 1;
                        }
                    }
                    fset(x, jx, temp);
                    jx -= incx;
                    kk -= j;
                }
            }
        } else {
            kk = 1;
            if incx == 1 {
                for j:1..n {
                    temp = fget(x, j);
                    k = kk + 1;
                    if noconj {
                        if nounit then temp *= fget(ap, kk);
                        for i:j+1..n {
                            temp += fget(ap, k)*fget(x, i);
                            k += 1;
                        }
                    } else {
                        if nounit then temp *= conjugate(fget(ap, kk));
                        for i:j+1..n {
                            temp += conjugate(fget(ap, k))*fget(x, i);
                            k += 1;
                        }
                    }
                    fset(x, j, temp);
                    kk += n-j+1;
                }
            } else {
                jx = kx;
                for j:1..n {
                    temp = fget(x, jx);
                    ix = jx;
                    k = kk + 1;
                    if noconj {
                        if nounit then temp *= fget(ap, kk);
                        for i:j+1..n {
                            ix += incx;
                            temp += fget(ap, k)*fget(x, ix);
                            k += 1;
                        }
                    } else {
                        if nounit then temp *= conjugate(fget(ap, kk));
                        for i:j+1..n {
                            ix += incx;
                            temp += conjugate(fget(ap, k))*fget(x, ix);
                            k += 1;
                        }
                    }
                    fset(x, jx, temp);
                    jx += incx;
                    kk += n-j+1;
                }
            }
        }
    }
}
