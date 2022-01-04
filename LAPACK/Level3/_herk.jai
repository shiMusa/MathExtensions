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

ssyrk :: _herk;
dsyrk :: _herk;
cherk :: _herk;
zherk :: _herk;

_herk :: (
    $$uplo: u8,
    $$trans: u8,
    $$n: int,
    $$k: int,
    $$alpha: $U,
    $$a: $A/ArrayTypedD,
    $$lda: int,
    $$beta: U,
    $$c: *$C/ArrayType2D,
    $$ldc: int
) #modify {
    if is_non_real_number(U) { U = null; return; }    
} {
    T :: A.data_type;
    zero :: #run zero(T); // complex
    one  :: #run one(T);  // complex
    rzero :: #run zero(U);  // real
    rone  :: #run one(U);  // real

    // Set nrowa as the number of rows of a.
    nrowa := ifx char_eq(trans, #char "n") then n else k;
    upper := char_eq(uplo, #char "u");
    
    // Test the input parameters.
    #if CHECKS {
        #run {
            assert(T == C.data_type);
            assert(U == get_number_inner_type(T));
        }
        assert(upper || char_eq(uplo, #char "l"));
        assert(char_eq(trans, #char "n", #char "c"));
        assert(n >= 0);
        assert(k >= 0);
        assert(lda >= max(1,nrowa));
        assert(ldc >= max(1,n));
    }

    // Quick return if possible.
    if n == 0 || ((alpha == rzero || k == 0) && beta == rone) then return;

    // And when  alpha == rzero
    if alpha == rzero {
        if upper {
            if beta == rzero {
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
                    set(c,j,j, beta*real_nc(get(c,j,j)));
                }
            }
        } else {
            if beta == rzero {
                for j:0..n-1 {
                    for i:j..n-1 {
                        set(c,i,j, zero);
                    }
                }
            } else {
                for j:0..n-1 {
                    set(c,j,j, beta*real_nc(get(c,j,j)));
                    for i:j+1..n-1 {
                        set(c,i,j, beta*get(c,i,j));
                    }
                }
            }
        }
    }

    temp : T = ---;
    rtemp : U = ---;

    // Start the operations.
    if char_eq(trans, #char "n") {
        // Form  C := alpha*A*A**H + beta*C.
        if upper {
            for j:0..n-1 {
                if beta == rzero {
                    for i:0..j {
                        set(c,i,j, zero);
                    }
                } else if beta != rone {
                    for i:0..j-1 {
                        set(c,i,j, beta*get(c,i,j));
                    }
                    set(c,j,j, beta*real_nc(get(c,j,j)));
                } else {
                    set(c,j,j, real_nc(get(c,j,j)));
                }
                for l:0..k {
                    if get(a,j,l) != zero {
                        temp = alpha*conjugate(get(a,j,l));
                        for i:0..j-1 {
                            set(c,i,j, get(c,i,j) + temp*get(a,i,l));
                        }
                        set(c,j,j, real_nc(get(c,j,j)) + real_nc(temp*get(a,i,l)));
                    }
                }
            }
        } else {
            for j:0..n-1 {
                if beta == rzero {
                    for i:j..n-1 {
                        set(c,i,j, zero);
                    }
                } else if beta != rone {
                    set(c,j,j, beta*real_nc(get(c,j,j)));
                    for i:j+1..n-1 {
                        set(c,i,j, beta*get(c,i,j));
                    }
                } else {
                    set(c,j,j) = real_nc(get(c,j,j));
                }
                for l:0..k {
                    if get(a,j,l) != zero {
                        temp = alpha*conjugate(get(a,j,l));
                        set(c,j,j, real_nc(get(c,j,j)) + real_nc(temp*get(a,j,l)));
                        for i:j+1..n-1 {
                            set(c,i,j, get(c,i,j) + temp*get(a,i,l));
                        }
                    }
                }
            }
        }
    } else {
        // Form  C := alpha*A**H*A + beta*C.
        if upper {
            for j:0..n-1 {
                for i:0..j-1 {
                    temp = zero;
                    for l:0..k {
                        temp += conjugate(get(a,l,i))*get(a,l,j);
                    }
                    if beta == rzero {
                        set(c,i,j, alpha*temp);
                    } else {
                        set(c,i,j, alpha*temp + beta*get(c,i,j));
                    }
                }
                rtemp = rzero;
                for l:0..k {
                    rtemp += real(conjugate(get(a,l,j))*get(a,l,j));
                }
                if beta == rzero {
                    set(c,j,j, alpha*rtemp);
                } else {
                    set(c,j,j, alpha*rtemp + beta*real_nc(get(c,j,j)));
                }
            }
        } else {
            for j:0..n-1 {
                rtemp = rzero;
                for l:0..k {
                    rtemp += real(conjugate(get(a,l,j))*get(a,l,j));
                }
                if beta == rzero {
                    set(c,j,j, alpha*rtemp);
                } else {
                    set(c,j,j, alpha*rtemp + beta*real_nc(get(c,j,j)));
                }
                for i:j+1..n-1 {
                    temp = zero;
                    for l:0..k {
                        temp += conjugate(get(a,l,i))*get(a,l,j);
                    }
                    if beta == rzero {
                        set(c,i,j, alpha*temp);
                    } else {
                        set(c,i,j, alpha*temp + beta*get(c,i,j));
                    }
                }
            }
        }
    }
}