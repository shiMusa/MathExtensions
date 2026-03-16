// ? LAPACK *hemm() and *symm()
// ? http://www.netlib.org/lapack/explore-html/db/def/group__complex__blas__level3_gad2d1853a142397404eae974b6574ece3.html#gad2d1853a142397404eae974b6574ece3

// ? hemm = hermitian matrix * matrix
// ? symm = symmetric matrix * matrix

// ? original comment
// CHEMM  performs one of the matrix-matrix operations
//
//    C := alpha*A*B + beta*C,
//
// or
//
//    C := alpha*B*A + beta*C,
//
// where alpha and beta are scalars, A is an hermitian matrix and  B and
// C are m by n matrices.

// Level 3 Blas routine.
// 
// -- Written on 8-February-1989.
//    Jack Dongarra, Argonne National Laboratory.
//    Iain Duff, AERE Harwell.
//    Jeremy Du Croz, Numerical Algorithms Group Ltd.
//    Sven Hammarling, Numerical Algorithms Group Ltd.

ssymm :: #bake_arguments _hemm(hermitian=false);
dsymm :: #bake_arguments _hemm(hermitian=false);
csymm :: #bake_arguments _hemm(hermitian=false);
zsymm :: #bake_arguments _hemm(hermitian=false);
chemm :: #bake_arguments _hemm(hermitian=true);
zhemm :: #bake_arguments _hemm(hermitian=true);

_hemm :: (
    $$side: u8,
    $$uplo: u8,
    $$m: int,
    $$n: int,
    alpha: $T,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    beta: T,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $hermitian: bool = true
) {
    zero_t :: #run zero(T);
    one_t  :: #run one(T);

    // Set nrowa as the number of rows of a.
    nrowa := ifx char_eq(side, #char "l") then m else n;
    upper := char_eq(uplo, #char "u");
    
    // Test the input parameters.
    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == B.data_type);
            assert(T == C.data_type);
        }
        assert(char_eq(side, #char "l", #char "r"));
        assert(upper || char_eq(uplo, #char "l"));
        assert(m >= 0);
        assert(n >= 0);
        assert(lda >= max(1,nrowa));
        assert(ldb >= max(1,m));
        assert(ldc >= max(1,m));
    }

    // Quick return if possible.
    if m == 0 || n == 0 || (alpha == zero_t && beta == one_t) then return;

    // And when  alpha == zero
    if alpha == zero_t {
        if beta == zero_t {
            for j:1..n {
                for i:1..m {
                    fset(c, i, j, zero_t);
                }
            }
        } else {
            for j:1..n {
                for i:1..m {
                    fset(c, i, j, beta*fget(c.*, i, j));
                }
            }
        }
        return;
    }

    temp1, temp2 : T = ---;

    // Start the operations
    if char_eq(side, #char "l") {
        // Form  C := alpha*A*B + beta*C
        if upper {
            for j:1..n {
                for i:1..m {
                    temp1 = alpha*fget(b.*, i, j);
                    temp2 = zero_t;
                    for k:1..i-1 {
                        fset(c, k, j, fget(c.*, k, j) + temp1*fget(a, k, i));
                        #if hermitian {
                            temp2 += fget(b.*, k, j)*conjugate(fget(a, k, i));
                        } else {
                            temp2 += fget(b.*, k, j)*fget(a, k, i);
                        }
                    }
                    if beta == zero_t {
                        #if hermitian {
                            fset(c, i, j, temp1*real_nc(fget(a, i, i)) + alpha*temp2);
                        } else {
                            fset(c, i, j, temp1*fget(a, i, i) + alpha*temp2);
                        }
                    } else {
                        #if hermitian {
                            fset(c, i, j, beta*fget(c.*, i, j) + temp1*real_nc(fget(a, i, i)) + alpha*temp2);
                        } else {
                            fset(c, i, j, beta*fget(c.*, i, j) + temp1*fget(a, i, i) + alpha*temp2);
                        }
                    }
                }
            }
        } else {
            for j:1..n {
                for <i:1..m {
                    temp1 = alpha*fget(b.*, i, j);
                    temp2 = zero_t;
                    for k:i+1..m {
                        fset(c, k, j, fget(c.*, k, j) + temp1*fget(a, k, i));
                        #if hermitian {
                            temp2 += fget(b.*, k, j)*conjugate(fget(a, k, i));
                        } else {
                            temp2 += fget(b.*, k, j)*fget(a, k, i);
                        }
                    }
                    if beta == zero_t {
                        #if hermitian {
                            fset(c, i, j, temp1*real_nc(fget(a, i, i)) + alpha*temp2);
                        } else {
                            fset(c, i, j, temp1*fget(a, i, i) + alpha*temp2);
                        }
                    } else {
                        #if hermitian {
                            fset(c, i, j, beta*fget(c.*, i, j) + temp1*real_nc(fget(a, i, i)) + alpha*temp2);
                        } else {
                            fset(c, i, j, beta*fget(c.*, i, j) + temp1*fget(a, i, i) + alpha*temp2);
                        }
                    }
                }
            }
        }
    } else {
        // Form  C := alpha*B*A + beta*C
        for j:1..n {
            #if hermitian {
                temp1 = alpha*real_nc(fget(a, j, j));
            } else {
                temp1 = alpha*fget(a, j, j);
            }
            if beta == zero_t {
                for i:1..m {
                    fset(c, i, j, temp1*fget(b.*, i, j));
                }
            } else {
                for i:1..m {
                    fset(c, i, j, beta*fget(c.*, i, j) + temp1*fget(b.*, i, j));
                }
            }
            for k:1..j-1 {
                if upper {
                    temp1 = alpha*fget(a, k, j);
                } else {
                    #if hermitian {
                        temp1 = alpha*conjugate(fget(a, j, k));
                    } else {
                        temp1 = alpha*fget(a, j, k);
                    }
                }
                for i:1..m {
                    fset(c, i, j, fget(c.*, i, j) + temp1*fget(b.*, i, k));
                }
            }
            for k:j+1..n {
                if upper {
                    #if hermitian {
                        temp1 = alpha*conjugate(fget(a, j, k));
                    } else {
                        temp1 = alpha*fget(a, j, k);
                    }
                } else {
                    temp1 = alpha*fget(a, k, j);
                }
                for i:1..m {
                    fset(c, i, j, fget(c.*, i, j) + temp1*fget(b.*, i, k));
                }
            }
        }
    }
}
