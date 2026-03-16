// ? LAPACK *tbsv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga47f0128bb830b67b37b285af34eaad5f.html#ga47f0128bb830b67b37b285af34eaad5f

// ? tbsv = triangular band (matrix) solve vector

// ? original comment
// CTBSV  solves _one of the systems of equations
//
//    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
//
// where b and x are n element vectors and A is an n by n unit, or
// non-unit, upper or lower triangular band matrix, with ( k + 1 )
// diagonals.
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


stbsv :: _tbsv;
dtbsv :: _tbsv;
ctbsv :: _tbsv;
ztbsv :: _tbsv;


_tbsv :: (
    $$uplo: u8,
    $$trans: u8,
    $$diag: u8,
    $$n: int,
    $$k: int,
    a: $A/ArrayType2D,
    $$lda: int,
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
        assert(k >= 0);
        assert(lda >= k+1);
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
    ix, jx, kplus1, l : int = ---;

    // Start the operations. In this version the elements of A are
    // accessed by sequentially with _one pass through A.

    if char_eq(trans, #char "n") {
        // Form x := inv(A)*x.
        if char_eq(uplo, #char "u") {
            kplus1 = k + 1;
            if incx == 1 {
                for <j:1..n {
                    if fget(x, j) != _zero {
                        l = kplus1 - j;
                        if nounit then fset(x, j, fget(x.*, j)/fget(a, kplus1, j));
                        temp = fget(x, j);
                        if j > 1 for <i:max(1,j-k)..j-1 {
                            fset(x, i, fget(x.*, i) - temp*fget(a, l+i, j));
                        }
                    }
                }
            } else {
                kx += (n-1)*incx;
                jx = kx;
                for <j:1..n {
                    kx -= incx;
                    if fget(x, jx) != _zero {
                        ix = kx;
                        l = kplus1 - j;
                        if nounit then fset(x, jx, fget(x.*, jx)/fget(a, kplus1, j));
                        temp = fget(x, jx);
                        if j > 1 for <i:max(1,j-k)..j-1 {
                            fset(x, ix, fget(x.*, ix) - temp*fget(a, l+i, j));
                            ix -= incx;
                        }
                    }
                    jx -= incx;
                }
            }
        } else {
            if incx == 1 {
                for j:1..n {
                    if fget(x, j) != _zero {
                        l = 1 - j;
                        if nounit then fset(x, j, fget(x.*, j)/fget(a, 1, j));
                        temp = fget(x, j);
                        for i:j+1..min(n,j+k) {
                            fset(x, i, fget(x.*, i) - temp*fget(a, l+i, j));
                        }
                    }
                }
            } else {
                jx = kx;
                for j:1..n {
                    kx += incx;
                    if fget(x, jx) != _zero {
                        ix = kx;
                        l = 1 - j;
                        if nounit then fset(x, jx, fget(x.*, jx)/fget(a, 1, j));
                        temp = fget(x, jx);
                        for i:j+1..min(n,j+k) {
                            fset(x, ix, fget(x.*, ix) - temp*fget(a, l+i, j));
                            ix += incx;
                        }
                    }
                    jx += incx;
                }
            }
        }
    } else {
        // Form   x := inv(A**T)*x  or  x := inv(A**H)*x.
        if char_eq(uplo, #char "u") {
            kplus1 = k + 1;
            if incx == 1 {
                for j:1..n {
                    temp = fget(x, j);
                    l = kplus1 - j;
                    if noconj {
                        for i:max(1,j-k)..j-1 {
                            temp -= fget(a, l+i, j)*fget(x, i);
                        }
                        if nounit then temp /= fget(a, kplus1, j);
                    } else {
                        for i:max(1,j-k)..j-1 {
                            temp -= conjugate(fget(a, l+i, j))*fget(x, i);
                        }
                        if nounit then temp /= conjugate(fget(a, kplus1, j));
                    }
                    fset(x, j, temp);
                }
            } else {
                jx = kx;
                for j:1..n {
                    temp = fget(x, jx);
                    ix = kx;
                    l = kplus1 - j;
                    if noconj {
                        for i:max(1,j-k)..j-1 {
                            temp -= fget(a, l+i, j)*fget(x, ix);
                            ix += incx;
                        }
                        if nounit then temp /= fget(a, kplus1, j);
                    } else {
                        for i:max(1,j-k)..j-1 {
                            temp -= conjugate(fget(a, l+i, j))*fget(x, ix);
                            ix += incx;
                        }
                        if nounit then temp /= conjugate(fget(a, kplus1, j));
                    }
                    fset(x, jx, temp);
                    jx += incx;
                    if j > k then kx += incx;
                }
            }
        } else {
            if incx == 1 {
                for <j:1..n {
                    temp = fget(x, j);
                    l = 1 - j;
                    if noconj {
                        for <i:j+1..min(n,j+k) {
                            temp -= fget(a, l+i, j)*fget(x, i);
                        }
                        if nounit then temp /= fget(a, 1, j);
                    } else {
                        for <i:j+1..min(n,j+k) {
                            temp -= conjugate(fget(a, l+i, j))*fget(x, i);
                        }
                        if nounit then temp /= conjugate(fget(a, 1, j));
                    }
                    fset(x, j, temp);
                }
            } else {
                kx += (n-1)*incx;
                jx = kx;
                for <j:1..n {
                    temp = fget(x, jx);
                    ix = kx;
                    l = 1 - j;
                    if noconj {
                        for <i:j+1..min(n,j+k) {
                            temp -= fget(a, l+i, j)*fget(x, ix);
                            ix -= incx;
                        }
                        if nounit then temp /= fget(a, 1, j);
                    } else {
                        for <i:j+1..min(n,j+k) {
                            temp -= conjugate(fget(a, l+i, j))*fget(x, ix);
                            ix -= incx;
                        }
                        if nounit then temp /= conjugate(fget(a, 1, j));
                    }
                    fset(x, jx, temp);
                    jx -= incx;
                    if (n-j) >= k then kx -= incx;
                }
            }
        }
    }
}
