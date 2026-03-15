// ? LAPACK *her2k() and *syr2k()
// ? http://www.netlib.org/lapack/explore-html/db/def/group__complex__blas__level3_gaf5266b622e0fbbd972cfc2df3061984f.html#gaf5266b622e0fbbd972cfc2df3061984f

// ? her2k = hermitian (rank) 2k

// ? original comment
// CHER2K  performs one of the hermitian rank 2k operations
//
//    C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
//
// or
//
//    C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
//
// where  alpha and beta  are scalars with  beta  real,  C is an  n by n
// hermitian matrix and  A and B  are  n by k matrices in the first case
// and  k by n  matrices in the second case.

// Level 3 Blas routine.
// 
// -- Written on 8-February-1989.
//    Jack Dongarra, Argonne National Laboratory.
//    Iain Duff, AERE Harwell.
//    Jeremy Du Croz, Numerical Algorithms Group Ltd.
//    Sven Hammarling, Numerical Algorithms Group Ltd.
// 
// -- Modified 8-Nov-93 to set C(J,J) to REAL( C(J,J) ) when BETA = 1.
//    Ed Anderson, Cray Research Inc.

ssyr2k :: #bake_arguments _her2k(hermitian=true);
dsyr2k :: #bake_arguments _her2k(hermitian=true);
cher2k :: #bake_arguments _her2k(hermitian=true);
zher2k :: #bake_arguments _her2k(hermitian=true);
csyr2k :: #bake_arguments _her2k(hermitian=false);
zsyr2k :: #bake_arguments _her2k(hermitian=false);

_her2k :: (
    $$uplo: u8,
    $$trans: u8,
    $$n: int,
    $$k: int,
    alpha: $T,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: $U,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $hermitian: bool = true
) {
    zero :: #run zero(T);
    one  :: #run one(T); 
    zeroU :: #run zero(U); 
    oneU  :: #run one(U);  

    conj :: inline (val: T) -> T #expand {
        #if hermitian {
            return conjugate(val);
        } else {
            return val;
        }
    }
    re :: inline (val: T) -> T #expand {
        #if hermitian {
            return real_nc(val);
        } else {
            return val;
        }
    }

    // Set nrowa as the number of rows of a.
    nrowa := ifx char_eq(trans, #char "n") then n else k;
    upper := char_eq(uplo, #char "u");
    
    // Test the input parameters.
    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == B.data_type);
            assert(T == C.data_type);
            assert((U == get_number_inner_type(T) && hermitian) || (U == T && !hermitian));
        }
        assert(upper || char_eq(uplo, #char "l"));
        assert(char_eq(trans, #char "n", #char "c"));
        assert(n >= 0);
        assert(k >= 0);
        assert(lda >= max(1,nrowa));
        assert(ldb >= max(1,nrowa));
        assert(ldc >= max(1,n));
    }

    // Quick return if possible.
    if n == 0 || ((alpha == zero || k == 0) && beta == oneU) then return;

    // And when  alpha == zero
    if alpha == zero {
        if upper {
            if beta == zeroU {
                for j:1..n {
                    for i:1..j {
                        fset(c, i, j, zero);
                    }
                }
            } else {
                for j:1..n {
                    for i:1..j-1 {
                        fset(c, i, j, beta*fget(c.*, i, j));
                    }
                    fset(c, j, j, beta*re(fget(c.*, j, j)));
                }
            }
        } else {
            if beta == zeroU {
                for j:1..n {
                    for i:j..n {
                        fset(c, i, j, zero);
                    }
                }
            } else {
                for j:1..n {
                    fset(c, j, j, beta*re(fget(c.*, j, j)));
                    for i:j+1..n {
                        fset(c, i, j, beta*fget(c.*, i, j));
                    }
                }
            }
        }
        return;
    }

    // Start the operations.
    if char_eq(trans, #char "n") {
        // From  C := alpha*A*B**H + conjugate(alpha)*B*A**H + C.
        if upper {
            for j:1..n {
                if beta == zeroU {
                    for i:1..j {
                        fset(c, i, j, zero);
                    }
                } else if beta != oneU {
                    for i:1..j-1 {
                        fset(c, i, j, beta*fget(c.*, i, j));
                    }
                    fset(c, j, j, beta*re(fget(c.*, j, j)));
                } else {
                    fset(c, j, j, re(fget(c.*, j, j)));
                }
                for l:1..k {
                    if fget(a, j, l) != zero || fget(b, j, l) != zero {
                        temp1 = alpha*conj(fget(b, j, l));
                        temp2 = conj(alpha*fget(a, j, l));
                        for i:1..j-1 {
                            fset(c, i, j, fget(c.*, i, j) + fget(a, i, l)*temp1 + fget(b, i, l)*temp2);
                        }
                        fset(c, j, j, re(fget(c.*, j, j)) + re(fget(a, j, l)*temp1 + fget(b, j, l)*temp2));
                    }
                }
            }
        } else {
            for j:1..n {
                if beta == zeroU {
                    for i:j..n {
                        fset(c, i, j, zero);
                    }
                } else if beta != oneU {
                    for i:j+1..n {
                        fset(c, i, j, beta*fget(c.*, i, j));
                    }
                } else {
                    fset(c, j, j, beta*re(fget(c.*, j, j)));
                }
                for l:1..k {
                    if fget(a, j, l) != zero || fget(b, j, l) != zero {
                        temp1 = alpha*conj(fget(b, j, l));
                        temp2 = conj(alpha*fget(a, j, l));
                        for i:j+1..n {
                            fset(c, i, j, fget(c.*, i, j) + fget(a, i, l)*temp1 + fget(b, i, l)*temp2);
                        }
                        fset(c, j, j, re(fget(c.*, j, j)) + re(fget(a, j, l)*temp1 + fget(b, j, l)*temp2));
                    }
                }
            }
        }
    } else {
        // Form  C := alpha*A**H*B + conj(alpha)*B**H*A + C.
        if upper {
            for j:1..n {
                for i:1..j {
                    temp1 = zero;
                    temp2 = zero;
                    for l:1..k {
                        temp1 += conj(fget(a, l, i))*fget(b, l, j);
                        temp2 += conj(fget(b, l, i))*fget(a, l, j);
                    }
                    if i == j {
                        if beta == zeroU {
                            fset(c, j, j, re(alpha*temp1 + conj(alpha)*temp2));
                        } else {
                            fset(c, j, j, beta*re(fget(c.*, j, j)) + re(alpha*temp1 + conj(alpha)*temp2));
                        }
                    } else {
                        if beta == zeroU {
                            fset(c, i, j, alpha*temp1 + conj(alpha)*temp2);
                        } else {
                            fset(c, i, j, beta*fget(c.*, i, j) + alpha*temp1 + conj(alpha)*temp2);
                        }
                    }
                }
            }
        } else {
            for j:1..n {
                for i:j..n {
                    temp1 = zero;
                    temp2 = zero;
                    for l:1..k {
                        temp1 += conj(fget(a, l, i))*fget(b, l, j);
                        temp2 += conj(fget(b, l, i))*fget(a, l, j);
                    }
                    if i == j {
                        if beta == zeroU {
                            fset(c, j, j, re(alpha*temp1 + conj(alpha)*temp2));
                        } else {
                            fset(c, j, j, beta*re(fget(c.*, j, j)) + re(alpha*temp1 + conj(alpha)*temp2));
                        }
                    } else {
                        if beta == zeroU {
                            fset(c, i, j, alpha*temp1 + conj(alpha)*temp2);
                        } else {
                            fset(c, i, j, beta*fget(c.*, i, j) + alpha*temp1 + conj(alpha)*temp2);
                        }
                    }
                }
            }
        }
    }
}
