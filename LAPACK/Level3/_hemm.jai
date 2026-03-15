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

ssymm :: #bake_arguments _hemm(hermitian=true);
dsymm :: #bake_arguments _hemm(hermitian=true);
chemm :: #bake_arguments _hemm(hermitian=true);
zhemm :: #bake_arguments _hemm(hermitian=true);
csymm :: #bake_arguments _hemm(hermitian=false);
zsymm :: #bake_arguments _hemm(hermitian=false);

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
    zero :: #run zero(T);
    one  :: #run one(T);

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
    if m == 0 || n == 0 || (alpha == zero && beta == one) then return;

    // And when  alpha == zero
    if alpha == zero {
        if beta == zero {
            for j:1..n {
                for i:1..m {
                    fset(c, i, j, zero);
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
                    temp2 = zero;
                    for k:1..i-1 {
                        fset(c, k, j, fget(c.*, k, j) + temp1*fget(a, k, i));
                        temp2 += fget(b.*, k, j)*conj(fget(a, k, i));
                    }
                    if beta == zero {
                        fset(c, i, j, temp1*re(fget(a, i, i)) + alpha*temp2);
                    } else {
                        fset(c, i, j, beta*fget(c.*, i, j) + temp1*re(fget(a, i, i)) + alpha*temp2);
                    }
                }
            }
        } else {
            for j:1..n {
                for <i:1..m {
                    temp1 = alpha*fget(b.*, i, j);
                    temp2 = zero;
                    for k:i+1..m {
                        fset(c, k, j, fget(c.*, k, j) + temp1*fget(a, k, i));
                        temp2 += fget(b.*, k, j)*conj(fget(a, k, i));
                    }
                    if beta == zero {
                        fset(c, i, j, temp1*re(fget(a, i, i)) + alpha*temp2);
                    } else {
                        fset(c, i, j, beta*fget(c.*, i, j) + temp1*re(fget(a, i, i)) + alpha*temp2);
                    }
                }
            }
        }
    } else {
        // Form  C := alpha*B*A + beta*C
        for j:1..n {
            temp1 = alpha*re(fget(a, j, j));
            if beta == zero {
                for i:1..m {
                    fset(c, i, j, temp1*fget(b.*, i, j));
                }
            } else {
                for i:1..m {
                    fset(c, i, j, beta*fget(c.*, i, j) + temp1*fget(b.*, i, j));
                }
            }
            for k:1..j-1 {
                temp1 = ifx upper then alpha*fget(a, k, j) else alpha*conj(fget(a, j, k));
                for i:1..m {
                    fset(c, i, j, fget(c.*, i, j) + temp1*fget(b.*, i, k));
                }
            }
            for k:j+1..n {
                temp1 = ifx upper then alpha*conj(fget(a, j, k)) else alpha*fget(a, k, j);
                for i:1..m {
                    fset(c, i, j, fget(c.*, i, j) + temp1*fget(b.*, i, k));
                }
            }
        }
    }
}
