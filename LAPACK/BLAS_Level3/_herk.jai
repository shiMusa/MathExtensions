// ? LAPACK *herk() and *syrk()
// ? http://www.netlib.org/lapack/explore-html/db/def/group__complex__blas__level3_gade9f14cf41f0cefea7918d716f3e1c20.html#gade9f14cf41f0cefea7918d716f3e1c20

// ? herk = hermitian rank k

// ? original comment
// CHERK  performs one of the hermitian rank k operations
//
//    C := alpha*A*A**H + beta*C,
//
// or
//
//    C := alpha*A**H*A + beta*C,
//
// where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
// matrix and  A  is an  n by k  matrix in the  first case and a  k by n
// matrix in the second case.

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

ssyrk :: #bake_arguments _herk(hermitian = false);
dsyrk :: #bake_arguments _herk(hermitian = false);
cherk :: #bake_arguments _herk(hermitian = true);
zherk :: #bake_arguments _herk(hermitian = true);
csyrk :: #bake_arguments _herk(hermitian = false);
zsyrk :: #bake_arguments _herk(hermitian = false);

_herk :: (
    $$uplo: u8,
    $$trans: u8,
    $$n: int,
    $$k: int,
    alpha: $U,
    a: $A/ArrayType2D,
    $$lda: int,
    beta: U,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $hermitian: bool = true
) {
    T :: A.data_type;
    zero_t :: #run zero(T);
    one_t  :: #run one(T); 
    zeroU :: #run zero(U);
    oneU  :: #run one(U); 

    conj :: inline (val: $V) -> V #expand {
        #if hermitian {
            return conjugate(val);
        } else {
            return val;
        }
    }
    re :: inline (val: $V) -> V #expand {
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
            assert(T == C.data_type);
            assert((U == get_number_inner_type(T) && hermitian) || (U == T && ! hermitian));
        }
        assert(upper || char_eq(uplo, #char "l"));
        #if hermitian {
            assert(char_eq(trans, #char "n", #char "c"));
        } else {
            assert(char_eq(trans, #char "n", #char "t"));
        }
        assert(n >= 0);
        assert(k >= 0);
        assert(lda >= max(1,nrowa));
        assert(ldc >= max(1,n));
    }

    // Quick return if possible.
    if n == 0 || ((alpha == zeroU || k == 0) && beta == oneU) then return;

    // And when  alpha == zeroU
    if alpha == zeroU {
        if upper {
            if beta == zeroU {
                for j:1..n {
                    for i:1..j {
                        fset(c, i, j, zero_t);
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
                        fset(c, i, j, zero_t);
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

    temp : T = ---;
    rtemp : U = ---;

    // Start the operations.
    if char_eq(trans, #char "n") {
        // Form  C := alpha*A*A**H + beta*C.
        if upper {
            for j:1..n {
                if beta == zeroU {
                    for i:1..j {
                        fset(c, i, j, zero_t);
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
                    if fget(a, j, l) != zero_t {
                        temp = alpha*conj(fget(a, j, l));
                        for i:1..j-1 {
                            fset(c, i, j, fget(c.*, i, j) + temp*fget(a, i, l));
                        }
                        fset(c, j, j, re(fget(c.*, j, j)) + re(temp*fget(a, j, l)));
                    }
                }
            }
        } else {
            for j:1..n {
                if beta == zeroU {
                    for i:j..n {
                        fset(c, i, j, zero_t);
                    }
                } else if beta != oneU {
                    fset(c, j, j, beta*re(fget(c.*, j, j)));
                    for i:j+1..n {
                        fset(c, i, j, beta*fget(c.*, i, j));
                    }
                } else {
                    fset(c, j, j, re(fget(c.*, j, j)));
                }
                for l:1..k {
                    if fget(a, j, l) != zero_t {
                        temp = alpha*conj(fget(a, j, l));
                        fset(c, j, j, re(fget(c.*, j, j)) + re(temp*fget(a, j, l)));
                        for i:j+1..n {
                            fset(c, i, j, fget(c.*, i, j) + temp*fget(a, i, l));
                        }
                    }
                }
            }
        }
    } else {
        // Form  C := alpha*A**H*A + beta*C.
        if upper {
            for j:1..n {
                for i:1..j-1 {
                    temp = zero_t;
                    for l:1..k {
                        temp += conj(fget(a, l, i))*fget(a, l, j);
                    }
                    if beta == zeroU {
                        fset(c, i, j, alpha*temp);
                    } else {
                        fset(c, i, j, alpha*temp + beta*fget(c.*, i, j));
                    }
                }
                #if hermitian {
                    rtemp = zeroU;
                    for l:1..k {
                        rtemp += real(conj(fget(a, l, j))*fget(a, l, j));
                    }
                    if beta == zeroU {
                        fset(c, j, j, cpx(alpha*rtemp, zeroU));
                    } else {
                        fset(c, j, j, cpx(alpha*rtemp + beta*real(fget(c.*, j, j)), zeroU));
                    }
                } else {
                    temp = zero_t;
                    for l:1..k {
                        temp += conj(fget(a, l, j))*fget(a, l, j);
                    }
                    if beta == zeroU {
                        fset(c, j, j, alpha*temp);
                    } else {
                        fset(c, j, j, alpha*temp + beta*re(fget(c.*, j, j)));
                    }
                }
            }
        } else {
            for j:1..n {
                #if hermitian {
                    rtemp = zeroU;
                    for l:1..k {
                        rtemp += real(conj(fget(a, l, j))*fget(a, l, j));
                    }
                    if beta == zeroU {
                        fset(c, j, j, cpx(alpha*rtemp, zeroU));
                    } else {
                        fset(c, j, j, cpx(alpha*rtemp + beta*real(fget(c.*, j, j)), zeroU));
                    }
                } else {
                    temp = zero_t;
                    for l:1..k {
                        temp += conj(fget(a, l, j))*fget(a, l, j);
                    }
                    if beta == zeroU {
                        fset(c, j, j, alpha*temp);
                    } else {
                        fset(c, j, j, alpha*temp + beta*re(fget(c.*, j, j)));
                    }
                }
                for i:j+1..n {
                    temp = zero_t;
                    for l:1..k {
                        temp += conj(fget(a, l, i))*fget(a, l, j);
                    }
                    if beta == zeroU {
                        fset(c, i, j, alpha*temp);
                    } else {
                        fset(c, i, j, alpha*temp + beta*fget(c.*, i, j));
                    }
                }
            }
        }
    }
}
