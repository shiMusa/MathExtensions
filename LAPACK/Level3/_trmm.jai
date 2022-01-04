// ? LAPACK *trmm()
// ? http://www.netlib.org/lapack/explore-html/db/def/group__complex__blas__level3_gad7c297c05b482699b6d60a29c8d4a165.html#gad7c297c05b482699b6d60a29c8d4a165

// ? trmm = triangular matrix * matrix

// ? original comment
// CTRMM  performs one of the matrix-matrix operations
//
//    B := alpha*op( A )*B,   or   B := alpha*B*op( A )
//
// where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
// non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
//
//    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.

// Level 3 Blas routine.
// 
// -- Written on 8-February-1989.
//    Jack Dongarra, Argonne National Laboratory.
//    Iain Duff, AERE Harwell.
//    Jeremy Du Croz, Numerical Algorithms Group Ltd.
//    Sven Hammarling, Numerical Algorithms Group Ltd.

strmm :: _trmm;
dtrmm :: _trmm;
ctrmm :: _trmm;
ztrmm :: _trmm;

_trmm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    $$alpha: $T,
    $$a: $A/ArrayType2D,
    $$lda: int,
    $$b: *$B/ArrayType2D,
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

    // And when alpha == zero
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
            // Form  B := alpha*A*B
            if upper {
                for j:0..n-1 {
                    for k:0..m-1 {
                        if get(b,k,j) != zero {
                            temp = alpha*get(b,k,j);
                            for i:0..k-1 {
                                set(b,i,j, get(b,i,j) + temp*get(a,i,k));
                            }
                            if nounit then temp *= get(a,k,k);
                            set(b,k,j, temp);
                        }
                    }
                }
            } else {
                for j:0..n-1 {
                    for <k:m-1..0 {
                        if get(b,k,j) != zero {
                            temp = alpha*get(b,k,j);
                            set(b,k,j, temp);
                            if nounit then set(b,k,j, get(b,k,j)*get(a,k,k));
                            for i:k+1..m-1 {
                                set(b,i,j, get(b,i,j) + temp*get(a,i,k));
                            }
                        }
                    }
                }
            }
        } else {
            // Form  B := alpha*A**T*B  or  B := alpha*A**H*B.
            if upper {
                for j:0..n-1 {
                    for <i:m-1..0 {
                        temp = get(b,i,j);
                        if noconj {
                            if nounit then temp *= get(a,i,i);
                            for k:0..i-1 {
                                temp += get(a,k,i)*get(b,k,j);
                            }
                        } else {
                            if nounit then temp *= conjugate(get(a,i,i));
                            for k:0..i-1 {
                                temp += conjugate(get(a,k,i))*get(b,k,j);
                            }
                        }
                        set(b,i,j, alpha*temp);
                    }
                }
            } else {
                for j:0..n-1 {
                    for i:0..m-1 {
                        temp = get(b,i,j);
                        if noconj {
                            if nounit then temp *= get(a,i,i);
                            for k:i+1..m-1 {
                                temp += get(a,k,i)*get(b,k,j);
                            }
                        } else {
                            if nounit then temp *= conjugate(get(a,i,i));
                            for k:i+1..m-1 {
                                temp += conjugate(get(a,k,i))*get(b,k,j);
                            }
                        }
                        set(b,i,j, alpha*temp);
                    }
                }
            }
        }
    } else {
        if char_eq(transa, "n") {
            // Form  B := alpha*B*A.
            if upper {
                for <j:n-1..0 {
                    temp = alpha;
                    if nounit then temp *= get(a,j,j);
                    for i:0..m-1 {
                        set(b,i,j, temp*get(b,i,j));
                    }
                    for k:0..j-1 {
                        if get(a,k,j) != zero {
                            temp = alpha*get(a,k,j);
                            for i:0..m-1 {
                                set(b,i,j, get(b,i,j) + temp*get(b,i,k));
                            }
                        }
                    }
                }
            } else {
                for j:0..n-1 {
                    temp = alpha;
                    if nounit then temp *= get(a,j,j);
                    for i:0..m-1 {
                        set(b,i,j, temp*get(b,i,j));
                    }
                    for k:j+1..n-1 {
                        if get(a,k,j) != zero {
                            temp = alpha*get(a,k,j);
                            for i:0..m-1 {
                                set(b,i,j, get(b,i,j) + temp*get(b,i,k));
                            }
                        }
                    }
                }
            }
        } else {
            // Form  B := alpha*B*A**T  or  B := alpha*B*A**H.
            if upper {
                for k:0..n-1 {
                    for j:0..k-1 {
                        if get(a,j,k) != zero {
                            if noconj {
                                temp = alpha*get(a,j,k);
                            } else {
                                temp = alpha*conjugate(get(a,j,k));
                            }
                            for i:0..m-1 {
                                set(b,i,j, get(b,i,j) + temp*get(b,i,k));
                            }
                        }
                    }
                    temp = alpha;
                    if nounit {
                        if noconj {
                            temp *= get(a,k,k);
                        } else {
                            temp *= conjugate(get(a,k,k));
                        }
                    }
                    if temp != one {
                        for i:0..m-1 {
                            set(b,i,k, temp*get(b,i,k));
                        }
                    }
                }
            } else {
                for <k:n-1..0 {
                    for j:k+1..n-1 {
                        if get(a,j,k) != zero {
                            temp = ifx noconj then alpha*get(a,j,k) else alpha*conjugate(get(a,j,k));
                            for i:0..m-1 {
                                set(b,i,j, get(b,i,j) + temp*get(b,i,k));
                            }
                        }
                    }
                    temp = alpha;
                    if nounit {
                        temp = ifx noconj then temp*get(a,k,k) else temp*conjugate(get(a,k,k));
                    }
                    if temp != one {
                        for i:0..m-1 {
                            set(b,i,k, temp*get(b,i,k));
                        }
                    }
                }
            }
        }
    }
}