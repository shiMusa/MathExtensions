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
        for j:0..n-1 {
            for i:0..m-1 {
                set(b,i,j, zero);
            }
        }
    }

    temp : T = ---;

    // Start the operations.
    if lside {
        if char_eq(transa, #char "n") {
            // Form  B := alpha*inv(A)*B.
            if upper {
                for j:0..n-1 {
                    if alpha != one {
                        for i:0..m-1 {
                            set(b,i,j, alpha*get(b,i,j));
                        }
                    }
                    for <k:m-1..0 {
                        if get(b,k,j) != zero {
                            if nounit then set(b,k,j, get(b,k,j)/get(a,k,k));
                            for i:0..k-1 {
                                set(b,i,j, get(b,i,j) - get(b,k,j)*get(a,i,k));
                            }
                        }
                    }
                }
            } else {
                for j:0..n-1 {
                    if alpha != one {
                        for i:0..m-1 {
                            set(b,i,j, alpha*get(b,i,j));
                        }
                    }
                    for k:0..m-1 {
                        if get(b,k,j) != zero {
                            if nounit then set(b,k,j, get(b,k,j)/get(a,k,k));
                            for i:k+1..m-1 {
                                set(b,i,j, get(b,i,j) - get(b,k,j)*get(a,i,k));
                            }
                        }
                    }
                }
            }
        } else {
            // Form  B := alpha*inv(A**T)*B
            // or    B := alpha*inv(A**H)*B.
            if upper {
                for j:0..n-1 {
                    for i:0..m-1 {
                        temp = alpha*get(b,i,j);
                        if noconj {
                            for k:0..i-1 {
                                temp -= get(a,k,i)*get(b,k,j);
                            }
                            if nounit then temp /= get(a,i,i);
                        } else {
                            for k:0..i-1 {
                                temp -= conjugate(get(a,k,i))*get(b,k,j);
                            }
                            if nounit then temp /= conjugate(get(a,i,i));
                        }
                        set(b,i,j, temp);
                    }
                }
            } else {
                for j:0..n-1 {
                    for <i:m-1..0 {
                        temp = alpha*get(b,i,j);
                        if noconj {
                            for k:i+1..m-1 {
                                temp -= get(a,k,i)*get(b,k,j);
                            }
                            if nounit then temp /= get(a,i,i);
                        } else {
                            for k:i+1..m-1 {
                                temp -= conjugate(get(a,k,i))*get(b,k,j);
                            }
                            if nounit then temp /= conjugate(get(a,i,i));
                        }
                        set(b,i,j, temp);
                    }
                }
            }
        }
    } else {
        if char_eq(transa, #char "n") {
            // Form  B := alpha*B*inv(A).
            if upper {
                for j:0..n-1 {
                    if alpha != one {
                        for i:0..m-1 {
                            set(b,i,j, alpha*get(b,i,j));
                        }
                    }
                    for k:0..j-1 {
                        if get(a,k,j) != zero {
                            for i:0..m-1 {
                                set(b,i,j, get(b,i,j) - get(a,k,j)*get(b,i,k));
                            }
                        }
                    }
                    if nounit {
                        temp = one/get(a,j,j);
                        for i:0..m-1 {
                            set(b,i,j, temp*get(b,i,j));
                        }
                    }
                }
            } else {
                for <j:n-1..0 {
                    if alpha != one {
                        for i:0..m-1 {
                            set(b,i,j, alpha*get(b,i,j));
                        }
                    }
                    for k:j+1..n-1 {
                        if get(a,k,j) != zero {
                            for i:0..m-1 {
                                set(b,i,j, get(b,i,j) - get(a,k,j)*get(b,i,k));
                            }
                        }
                    }
                    if nounit {
                        temp = one/get(a,j,j);
                        for i:0..m-1 {
                            set(b,i,j, temp*get(b,i,j));
                        }
                    }
                }
            }
        } else {
            // Form  B := alpha*B*inv(A**T)
            // or    B := alpha*B*inv(A**H).
            if upper {
                for <k:n-1..0 {
                    if nounit {
                        temp = ifx noconj then one/get(a,k,k) else one/conjugate(get(a,k,k));
                        for i:0..m-1 {
                            set(b,i,k, temp*get(b,i,k));
                        }
                    }
                    for j:0..k-1 {
                        if get(a,j,k) != zero {
                            temp = ifx noconj then get(a,j,k) else conjugate(get(a,j,k));
                            for i:0..m-1 {
                                set(b,i,j, get(b,i,j) - temp*get(b,i,k));
                            }
                        }
                    }
                    if alpha != one {
                        for i:0..m-1 {
                            set(b,i,k, alpha*get(b,i,k));
                        }
                    }
                }
            } else {
                for k:0..n-1 {
                    if nounit {
                        temp = ifx noconj then one/get(a,k,k) else one/conjugate(get(a,k,k));
                        for i:0..m-1 {
                            set(b,i,k, temp*get(b,i,k));
                        }
                    }
                    for j:k+1..n-1 {
                        if get(a,j,k) != zero {
                            temp = ifx noconj then get(a,j,k) else conjugate(get(a,j,k));
                            for i:0..m-1 {
                                set(b,i,j, get(b,i,j) - temp*get(b,i,k));
                            }
                        }
                    }
                    if alpha != one {
                        for i:0..m-1 {
                            set(b,i,k, alpha*get(b,i,k));
                        }
                    }
                }
            }
        }
    }
}