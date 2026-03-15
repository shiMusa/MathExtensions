// ? LAPACK *trsm()
// ? http://www.netlib.org/lapack/explore-html/db/def/group__complex__blas__level3_gaf33844c7fd27e5434496d2ce0c1fc9d4.html#gaf33844c7fd27e5434496d2ce0c1fc9d4

// ? trsm = triangular solve matrix

// ? original comment
// CTRSM  solves one of the matrix equations
//
//    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
//
// where alpha is a scalar, X and B are m by n matrices, A is a unit, or
// non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
//
//    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
//
// The matrix X is overwritten on B.

// Level 3 Blas routine.
// 
// -- Written on 8-February-1989.
//    Jack Dongarra, Argonne National Laboratory.
//    Iain Duff, AERE Harwell.
//    Jeremy Du Croz, Numerical Algorithms Group Ltd.
//    Sven Hammarling, Numerical Algorithms Group Ltd.

strsm :: _trsm;
dtrsm :: _trsm;
ctrsm :: _trsm;
ztrsm :: _trsm;

_trsm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    alpha: $T,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int
) {
    zero :: #run zero(T);
    one  :: #run one(T);

    // Set nrowa as the number of rows of a.
    lside := char_eq(side, "l");
    nrowa := ifx lside then m else n;
    noconj := char_eq(transa, #char "t");
    nounit := char_eq(diag, #char "n");
    upper := char_eq(uplo, #char "u");
    
    // Test the input parameters.
    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == B.data_type);
        }
        assert(lside || char_eq(side, #char "r"));
        assert(upper || char_eq(uplo, #char "l"));
        assert(char_eq(transa, #char "n", #char "t", #char "c"));
        assert(char_eq(diag, #char "u", #char "n"));
        assert(m >= 0);
        assert(n >= 0);
        assert(lda >= max(1,nrowa));
        assert(ldb >= max(1,m));
    }

    // Quick return if possible.
    if m == 0 || n == 0 then return;

    // And when  alpha == zero
    if alpha == zero {
        for j:1..n {
            for i:1..m {
                fset(b, i, j, zero);
            }
        }
        return;
    }

    temp : T = ---;

    // Start the operations.
    if lside {
        if char_eq(transa, #char "n") {
            // Form  B := alpha*inv(A)*B.
            if upper {
                for j:1..n {
                    if alpha != one {
                        for i:1..m {
                            fset(b, i, j, alpha*fget(b, i, j));
                        }
                    }
                    for <k:1..m {
                        if fget(b, k, j) != zero {
                            if nounit then fset(b, k, j, fget(b, k, j)/fget(a, k, k));
                            for i:1..k-1 {
                                fset(b, i, j, fget(b, i, j) - fget(b, k, j)*fget(a, i, k));
                            }
                        }
                    }
                }
            } else {
                for j:1..n {
                    if alpha != one {
                        for i:1..m {
                            fset(b, i, j, alpha*fget(b, i, j));
                        }
                    }
                    for k:1..m {
                        if fget(b, k, j) != zero {
                            if nounit then fset(b, k, j, fget(b, k, j)/fget(a, k, k));
                            for i:k+1..m {
                                fset(b, i, j, fget(b, i, j) - fget(b, k, j)*fget(a, i, k));
                            }
                        }
                    }
                }
            }
        } else {
            // Form  B := alpha*inv(A**T)*B
            // or    B := alpha*inv(A**H)*B.
            if upper {
                for j:1..n {
                    for i:1..m {
                        temp = alpha*fget(b, i, j);
                        if noconj {
                            for k:1..i-1 {
                                temp -= fget(a, k, i)*fget(b, k, j);
                            }
                            if nounit then temp /= fget(a, i, i);
                        } else {
                            for k:1..i-1 {
                                temp -= conjugate(fget(a, k, i))*fget(b, k, j);
                            }
                            if nounit then temp /= conjugate(fget(a, i, i));
                        }
                        fset(b, i, j, temp);
                    }
                }
            } else {
                for j:1..n {
                    for <i:1..m {
                        temp = alpha*fget(b, i, j);
                        if noconj {
                            for k:i+1..m {
                                temp -= fget(a, k, i)*fget(b, k, j);
                            }
                            if nounit then temp /= fget(a, i, i);
                        } else {
                            for k:i+1..m {
                                temp -= conjugate(fget(a, k, i))*fget(b, k, j);
                            }
                            if nounit then temp /= conjugate(fget(a, i, i));
                        }
                        fset(b, i, j, temp);
                    }
                }
            }
        }
    } else {
        if char_eq(transa, #char "n") {
            // Form  B := alpha*B*inv(A).
            if upper {
                for j:1..n {
                    if alpha != one {
                        for i:1..m {
                            fset(b, i, j, alpha*fget(b, i, j));
                        }
                    }
                    for k:1..j-1 {
                        if fget(a, k, j) != zero {
                            for i:1..m {
                                fset(b, i, j, fget(b, i, j) - fget(a, k, j)*fget(b, i, k));
                            }
                        }
                    }
                    if nounit {
                        temp = one/fget(a, j, j);
                        for i:1..m {
                            fset(b, i, j, temp*fget(b, i, j));
                        }
                    }
                }
            } else {
                for <j:1..n {
                    if alpha != one {
                        for i:1..m {
                            fset(b, i, j, alpha*fget(b, i, j));
                        }
                    }
                    for k:j+1..n {
                        if fget(a, k, j) != zero {
                            for i:1..m {
                                fset(b, i, j, fget(b, i, j) - fget(a, k, j)*fget(b, i, k));
                            }
                        }
                    }
                    if nounit {
                        temp = one/fget(a, j, j);
                        for i:1..m {
                            fset(b, i, j, temp*fget(b, i, j));
                        }
                    }
                }
            }
        } else {
            // Form  B := alpha*B*inv(A**T)
            // or    B := alpha*B*inv(A**H).
            if upper {
                for <k:1..n {
                    if nounit {
                        temp = ifx noconj then one/fget(a, k, k) else one/conjugate(fget(a, k, k));
                        for i:1..m {
                            fset(b, i, k, temp*fget(b, i, k));
                        }
                    }
                    for j:1..k-1 {
                        if fget(a, j, k) != zero {
                            temp = ifx noconj then fget(a, j, k) else conjugate(fget(a, j, k));
                            for i:1..m {
                                fset(b, i, j, fget(b, i, j) - temp*fget(b, i, k));
                            }
                        }
                    }
                    if alpha != one {
                        for i:1..m {
                            fset(b, i, k, alpha*fget(b, i, k));
                        }
                    }
                }
            } else {
                for k:1..n {
                    if nounit {
                        temp = ifx noconj then one/fget(a, k, k) else one/conjugate(fget(a, k, k));
                        for i:1..m {
                            fset(b, i, k, temp*fget(b, i, k));
                        }
                    }
                    for j:k+1..n {
                        if fget(a, j, k) != zero {
                            temp = ifx noconj then fget(a, j, k) else conjugate(fget(a, j, k));
                            for i:1..m {
                                fset(b, i, j, fget(b, i, j) - temp*fget(b, i, k));
                            }
                        }
                    }
                    if alpha != one {
                        for i:1..m {
                            fset(b, i, k, alpha*fget(b, i, k));
                        }
                    }
                }
            }
        }
    }
}
