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

ssymm :: _hemm;
dsymm :: _hemm;
chemm :: _hemm;
zhemm :: _hemm;

_hemm :: (
    $$side: u8,
    $$uplo: u8,
    $$m: int,
    $$n: int,
    $$alpha: $T,
    $$a: $A/ArrayType2D,
    $$lda: int,
    $$b: *$B/ArrayType2D,
    $$ldb: int,
    $$beta: T,
    $$c: *$C/ArrayType2D,
    $$ldc: int
) {
    zero :: #run zero(T);
    one  :: #run one(T);

    // Set nrowa as the number of rows of a.
    nrowa := ifx char_eq(side, #char "l") then m else k;
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
            for j:0..n-1 {
                for i:0..m-1 {
                    set(c,i,j, zero);
                }
            }
        } else {
            for j:0..n-1 {
                for i:0..m-1 {
                    set(c,i,j, beta*get(c,i,j));
                }
            }
        }
    }

    temp1, temp2 : T = ---;

    // Start the operations
    if char_eq(side, #char "l") {
        // Form  C := alpha*A*B + beta*C
        if upper {
            for j:0..n-1 {
                for i:0..m-1 {
                    temp1 = alpha*get(b,i,j);
                    temp2 = zero;
                    for k:0..i-1 {
                        set(c,k,j, get(c,k,j) + temp1*get(a,k,i));
                        temp2 += get(b,k,j)*conjugate(get(a,k,i));
                    }
                    if beta == zero {
                        set(c,i,j, temp1*real_nc(get(a,i,i)) + alpha*temp2);
                    } else {
                        set(c,i,j, beta*get(c,i,j) + temp1*real_nc(get(a,i,j)) + alpha*temp2);
                    }
                }
            }
        } else {
            for j:0..n-1 {
                for <i:m-1..0{
                    temp1 = alpha*get(b,i,j);
                    temp2 = zero;
                    for k:i+1..m-1 {
                        set(c,k,j, get(c,k,j) + temp1*get(a,k,i));
                        temp2 += get(b,k,j)*conjugate(get(a,k,i));
                    }
                    if beta == zero {
                        set(c,i,j, temp1*real_nc(get(a,i,i)) + alpha*temp2);
                    } else {
                        set(c,i,j, beta*get(c,i,j) + temp1*real_nc(get(a,i,i)) + alpha*temp2);
                    }
                }
            }
        }
    } else {
        // Form  C := alpha*B*A + beta*C
        for j:0..n-1 {
            temp1 = alpha*real_nc(get(a,j,j));
            if beta == zero {
                for i:0..m-1 {
                    set(c,i,j, temp1*get(b,i,j));
                }
            } else {
                for i:0..m-1 {
                    set(c,i,j, beta*get(c,i,j) + temp1*get(b,i,j));
                }
            }
            for k:0..j-1 {
                temp1 = ifx upper then alpha*get(a,k,j) else alpha*conjugate(get(a,j,k));
                for i:0..m-1 {
                    set(c,i,j, get(c,i,j) + temp1*get(b,i,k));
                }
            }
            for k:j+1..n-1 {
                temp1 = ifx upper then alpha*conjugate(get(a,j,k)) else alpha*get(a,k,j);
                for i:0..m-1 {
                    set(c,i,j, get(c,i,j) + temp1*get(b,i,k));
                }
            }
        }
    }
}