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
                for j:0..n-1 {
                    for i:0..j {
                        set(c,i,j, zero);
                    }
                }
            } else {
                for j:0..n-1 {
                    for i:0..j-1 {
                        set(c,i,j, beta*get(c,i,j));
                    }
                    set(c,j,j, beta*re(get(c,j,j)));
                }
            }
        } else {
            if beta == zeroU {
                for j:0..n-1 {
                    for i:j..n-1 {
                        set(c,i,j, zero);
                    }
                }
            } else {
                for j:0..n-1 {
                    set(c,j,j, beta*re(get(c,j,j)));
                    for i:j+1..n-1 {
                        set(c,i,j, beta*get(c,i,j));
                    }
                }
            }
        }
    }

    // Start the operations.
    if char_eq(trans, #char "n") {
        // From  C := alpha*A*B**H + conjugate(alpha)*B*A**H + C.
        if upper {
            for j:0..n-1 {
                if beta == zeroU {
                    for i:0..j {
                        set(c,i,j, zero);
                    }
                } else if beta == oneU {
                    for i:0..j-1 {
                        set(c,i,j, beta*get(c,i,j));
                    }
                    set(c,j,j, beta*re(get(c,j,j)));
                } else {
                    set(c,j,j, re(get(c,j,j)));
                }
                for l:0..k {
                    if get(a,j,l) != zero || get(b,l) != zero {
                        temp1 = alpha*conj(get(b,j,l));
                        temp2 = conj(alpha*get(a,j,l));
                        for i:0..j-1 {
                            set(c,i,j, get(c,i,j) + get(a,i,l)*temp1 + get(b,i,l)*temp2);
                        }
                        set(c,j,j, re(get(c,j,j)) + 
                            re(get(a,j,l)*temp1 + get(b,j,l)*temp2));
                    }
                }
            }
        } else {
            for j:0..n-1 {
                if beta == zeroU {
                    for i:j..n-1 {
                        set(c,i,j, zero);
                    }
                } else if beta != oneU {
                    for i:j+1..n-1 {
                        set(c,i,j, beta*get(c,i,j));
                    }
                } else {
                    set(c,j,j, beta*re(get(c,j,j)));
                }
                for l:0..k {
                    if get(a,j,l) != zero || get(b,j,l) != zero {
                        temp1 = alpha*conj(get(b,j,l));
                        temp2 = conj(alpha*get(a,j,l));
                        for i:j+1..n-1 {
                            set(c,i,j, get(c,i,j) + get(a,i,l)*temp1 + get(b,i,l)*temp2);
                        }
                        set(c,j,j, re(get(c,j,j)) + 
                            re(get(a,j,l)*temp1 + get(b,j,l)*temp2));
                    }
                }
            }
        }
    } else {
        // Form  C := alpha*A**H*B + conj(alpha)*B**H*A + C.
        if upper {
            for j:0..n-1 {
                for i:0..j {
                    temp1 = zero;
                    temp2 = zero;
                    for l:0..k {
                        temp1 += conj(get(a,l,i))*get(b,l,j);
                        temp2 += conj(get(b,l,i))*get(a,l,j);
                    }
                    if i == j {
                        if beta == zeroU {
                            set(c,j,j, re( alpha*temp1 + conj(alpha)*temp2 ));
                        } else {
                            set(c,j,j, beta*re(get(c,j,j)) + 
                                re( alpha*temp1 + conj(alpha)*temp2 ));
                        }
                    } else {
                        if beta == zeroU {
                            set(c,i,j, alpha*temp1 + conj(alpha)*temp2);
                        } else {
                            set(c,i,j, beta*get(c,i,j) + alpha*temp1 + conj(alpha)*temp2);
                        }
                    }
                }
            }
        } else {
            for j:0..n-1 {
                for i:j..n-1 {
                    temp1 = zero;
                    temp2 = zero;
                    for l:0..k {
                        temp1 += conj(get(a,l,i))*get(b,l,j);
                        temp2 += conj(get(b,l,i))*get(a,l,j);
                    }
                    if i == j {
                        if beta == zeroU {
                            set(c,j,j, re( alpha*temp1 + conj(alpha*temp2) ));
                        } else {
                            set(c,j,j, beta*re(get(c,j,j)) + 
                                re( alpha*temp1 + conj(alpha)*temp2 ));
                        }
                    } else {
                        if beta == zeroU {
                            set(c,i,j, alpha*temp1 + conj(alpha)*temp2);
                        } else {
                            set(c,i,j, beta*get(c,i,j) + alpha*temp1 + conj(alpha)*temp2);
                        }
                    }
                }
            }
        }
    }
}