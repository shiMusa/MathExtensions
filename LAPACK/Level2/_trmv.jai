// ? LAPACK *trmv()
// ? http://www.netlib.org/lapack/explore-html/d6/dee/group__complex__blas__level2_ga62930e76219d094e424db1712a5356cb.html#ga62930e76219d094e424db1712a5356cb

// ? trmv = triangular matrix * vector

// ? original comment
// CTRMV  performs one of the matrix-vector operations
//
//    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit,
// upper or lower triangular matrix.

// Level 2 Blas routine.
// The vector and matrix arguments are not referenced when N = 0, or M = 0
// 
// -- Written on 22-October-1986.
//    Jack Dongarra, Argonne National Lab.
//    Jeremy Du Croz, Nag Central Office.
//    Sven Hammarling, Nag Central Office.
//    Richard Hanson, Sandia National Labs.

strvm :: _trmv;
dtrmv :: _trmv;
ctrmv :: _trmv;
ztrmv :: _trmv;

_trmv :: (
    $$uplo: u8,
    $$trans: u8,
    $$diag: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    x: *$X/ArrayType1D,
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
        assert(lda >= max(1,n));
        assert(incx != 0);
    }

    // Quick return if possible.
    if n == 0 then return;

    #if !is_non_real_number(T) {
        noconj :: true; // no complex or quaternion number, so no conjugation
    } else {
        #if is_constant(trans) {
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
    kx := ifx incx > 0 then 0 else (n-1)*incx;

    temp : T = ---;
    ix, jx, kx : int = ---;

    // Start the operations. In this version the lements of A are
    // accessed sequentially with one pass through A.

    if char_eq(trans, #char "n") {
        // From  x := A*x.
        if char_eq(uplo, #char "u") {
            if incx == 1 {
                for j:0..n-1 {
                    if get(x,j) != zero {
                        temp = get(x,j);
                        for i:0..j-1 {
                            set(x,i, get(x,i) + temp*get(a,i,j));
                        }
                        if nounit then set(x,j, get(x,j) * get(a,j,j));
                    }
                }
            } else {
                jx = kx;
                for j:0..n-1 {
                    if get(x,jx) != zero {
                        temp = get(x,jx);
                        ix = kx;
                        for i:0..j-1 {
                            set(x,ix, get(x,ix) + temp*get(a,i,j));
                            ix += incx;
                        }
                        if nounit then set(x,jx, get(x,jx)*get(a,j,j));
                    }
                    jx += incx;
                }
            }
        } else {
            if incx == 1 {
                for <j:n-1..0 {
                    if get(x,j) != zero {
                        temp = get(x,j);
                        for <i:n-1..j+1 {
                            set(x,i, get(x,i) + temp*get(a,i,j));
                        }
                        if nounit then set(x,j, get(x,j)*get(a,j,j));
                    }
                }
            } else {
                kx += (n-1)*incx;
                jx = kx;
                for <j:n-1..0 {
                    if get(x,jx) != zero {
                        temp = get(x,jx);
                        ix = kx;
                        for <i:n-1..j+1 {
                            set(x,ix, get(x,ix) + temp*get(a,i,j));
                            ix -= incx;
                        }
                        if nounit then set(x,jx, get(x,jx)*get(a,j,j));
                    }
                    jx -= incx;
                }
            }
        }
    } else {
        // Form  x := A**T*x  or  x := A**H*x.
        if char_eq(uplo, #char "u") {
            if incx == 1 {
                for <j:n-1..0 {
                    temp = get(x,j);
                    if noconj {
                        if nounit then tmp *= get(a,j,j);
                        for <i:j-1..0 {
                            temp += get(a,i,j)*get(x,i);
                        }
                    } else {
                        if nounit then temp *= conjugate(get(a,j,j));
                        for <i:j-1..0 {
                            temp += conjugate(get(a,i,j))*get(x,i);
                        }
                    }
                    set(x,j, temp);
                }
            } else {
                jx = kx + (n-1)*incx;
                for <j:n-1..0 {
                    temp = get(x,jx);
                    ix = jx;
                    if noconj {
                        if nounit then temp *= get(a,j,j);
                        for <i:j-1..0 {
                            ix -= incx;
                            temp += get(a,i,j)*get(x,ix);
                        }
                    } else {
                        if nounit then temp *= conjugate(get(a,j,j));
                        for <i:j-1..0 {
                            ix -= incx;
                            temp += conjugate(get(a,i,j))*get(x,ix);
                        }
                    }
                    set(x,jx, temp);
                    jx -= incx;
                }
            }
        } else {
            if incx == 1 {
                for j:0..n-1 {
                    temp = get(x,j);
                    if noconj {
                        if nounit then temp *= get(a,j,j);
                        for i:j+1..n-1 {
                            temp += get(a,i,j)*get(x,i);
                        }
                    } else {
                        if nounit then temp *= conjugate(get(a,j,j));
                        for i:j+1..n-1 {
                            temp += conjugate(get(a,i,j))*get(x,i);
                        }
                    }
                    set(x,j, temp);
                }
            } else {
                jx = kx;
                for j:0..n-1 {
                    temp = get(x,jx);
                    ix = jx;
                    if noconj {
                        if nounit then temp *= get(a,j,j);
                        for i:j+1..n-1 {
                            ix += incx;
                            temp += get(a,i,j)*get(x,ix);
                        }
                    } else {
                        if nounit then temp *= conjugate(get(a,j,j));
                        for i:j+1..n-1 {
                            ix += incx;
                            temp += conjugate(get(a,i,j))*get(x,ix);
                        }
                    }
                    set(x,jx, temp);
                    jx += incx;
                }
            }
        }
    }
}