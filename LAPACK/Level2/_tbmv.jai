// ? LAPACK *tbmv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga2483c6533c8a24ff715137629f939074.html#ga2483c6533c8a24ff715137629f939074

// ? tbmv = triangular band matrix * vector

// ? original comment
// CTBMV  performs _one of the matrix-vector operations
//
//    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit,
// upper or lower triangular band matrix, with ( k + 1 ) diagonals.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

stbmv :: _tbmv;
dtbmv :: _tbmv;
ctbmv :: _tbmv;
ztbmv :: _tbmv;

_tbmv :: (
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

    // quick return if possible.
    if n == 0 then return;

    noconj := char_eq(trans, #char "t");
    nounit := char_eq(diag, #char "n");

    temp : T = ---;
    ix, jx, kplus1, l : int = ---;

    // Set up the start point in x if the increment is not unity. This
    // will be  (n-1)*incx  too small for descending loops.
    kx := ifx incx > 0 then 1 else 1 - (n-1)*incx;

    // Start the operations. In this version the elements of A are
    // accessed sequentially with _one pass through A.

    if char_eq(trans, #char "n") {
        // Form  x := A*x
        if char_eq(uplo, #char "u") {
            kplus1 = k + 1;
            if incx == 1 {
                for j:1..n {
                    if fget(x, j) != _zero {
                        temp = fget(x, j);
                        l = kplus1 - j;
                        for i:max(1,j-k)..j-1 {
                            fset(x, i, fget(x.*, i) + temp*fget(a, l+i, j));
                        }
                        if nounit then fset(x, j, fget(x.*, j)*fget(a, kplus1, j));
                    }
                }
            } else {
                jx = kx;
                for j:1..n {
                    if fget(x, jx) != _zero {
                        temp = fget(x, jx);
                        ix = kx;
                        l = kplus1 - j;
                        for i:max(1,j-k)..j-1 {
                            fset(x, ix, fget(x.*, ix) + temp*fget(a, l+i, j));
                            ix += incx;
                        }
                        if nounit then fset(x, jx, fget(x.*, jx)*fget(a, kplus1, j));
                    }
                    jx += incx;
                    if j > k then kx += incx;
                }
            }
        } else {
            if incx == 1 {
                for <j:1..n {
                    if fget(x, j) != _zero {
                        temp = fget(x, j);
                        l = 1 - j;
                        for <i:j+1..min(n,j+k) {
                            fset(x, i, fget(x.*, i) + temp*fget(a, l+i, j));
                        }
                        if nounit then fset(x, j, fget(x.*, j)*fget(a, 1, j));
                    }
                }
            } else {
                kx += (n-1)*incx;
                jx = kx;
                for <j:1..n {
                    if fget(x, jx) != _zero {
                        temp = fget(x, jx);
                        ix = kx;
                        l = 1 - j;
                        for <i:j+1..min(n,j+k) {
                            fset(x, ix, fget(x.*, ix) + temp*fget(a, l+i, j));
                            ix -= incx;
                        }
                        if nounit then fset(x, jx, fget(x.*, jx)*fget(a, 1, j));
                    }
                    jx -= incx;
                    if (n-j) >= k then kx -= incx;
                }
            }
        }
    } else {
        // Form  x := A**T*x  or  x := A**H*x.
        if char_eq(uplo, #char "u") {
            kplus1 = k + 1;
            if incx == 1 {
                for <j:1..n {
                    temp = fget(x, j);
                    l = kplus1 - j;
                    if noconj {
                        if nounit then temp *= fget(a, kplus1, j);
                        if j > 1 for <i:max(1,j-k)..j-1 {
                            temp += fget(a, l+i, j)*fget(x, i);
                        }
                    } else {
                        if nounit then temp *= conjugate(fget(a, kplus1, j));
                        if j > 1 for <i:max(1,j-k)..j-1 {
                            temp += conjugate(fget(a, l+i, j))*fget(x, i);
                        }
                    }
                    fset(x, j, temp);
                }
            } else {
                kx += (n-1)*incx;
                jx = kx;
                for <j:1..n {
                    temp = fget(x, jx);
                    kx -= incx;
                    ix = kx;
                    l = kplus1 - j;
                    if noconj {
                        if nounit then temp *= fget(a, kplus1, j);
                        if j > 1 for <i:max(1,j-k)..j-1 {
                            temp += fget(a, l+i, j)*fget(x, ix);
                            ix -= incx;
                        }
                    } else {
                        if nounit then temp *= conjugate(fget(a, kplus1, j));
                        if j > 1 for <i:max(1,j-k)..j-1 {
                            temp += conjugate(fget(a, l+i, j))*fget(x, ix);
                            ix -= incx;
                        }
                    }
                    fset(x, jx, temp);
                    jx -= incx;
                }
            }
        } else {
            if incx == 1 {
                for j:1..n {
                    temp = fget(x, j);
                    l = 1 - j;
                    if noconj {
                        if nounit then temp *= fget(a, 1, j);
                        for i:j+1..min(n,j+k) {
                            temp += fget(a, l+i, j)*fget(x, i);
                        }
                    } else {
                        if nounit then temp *= conjugate(fget(a, 1, j));
                        for i:j+1..min(n,j+k) {
                            temp += conjugate(fget(a, l+i, j))*fget(x, i);
                        }
                    }
                    fset(x, j, temp);
                }
            } else {
                jx = kx;
                for j:1..n {
                    temp = fget(x, jx);
                    kx += incx;
                    ix = kx;
                    l = 1 - j;
                    if noconj {
                        if nounit then temp *= fget(a, 1, j);
                        for i:j+1..min(n,j+k) {
                            temp += fget(a, l+i, j)*fget(x, ix);
                            ix += incx;
                        }
                    } else {
                        if nounit then temp *= conjugate(fget(a, 1, j));
                        for i:j+1..min(n,j+k) {
                            temp += conjugate(fget(a, l+i, j))*fget(x, ix);
                            ix += incx;
                        }
                    }
                    fset(x, jx, temp);
                    jx += incx;
                }
            }
        }
    }
}
