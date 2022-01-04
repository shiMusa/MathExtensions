// ? LAPACK * gemm()
// ? http://www.netlib.org/lapack/explore-html/db/def/group__complex__blas__level3_gac4e11e8e8a4b0c802d2229b08da341f6.html#gac4e11e8e8a4b0c802d2229b08da341f6

// ? gemm = generalized matrix * matrix

// ? original comment
// CGEMM  performs one of the matrix-matrix operations
//
//    C := alpha*op( A )*op( B ) + beta*C,
//
// where  op( X ) is one of
//
//    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
//
// alpha and beta are scalars, and A, B and C are matrices, with op( A )
// an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

// Level 3 Blas routine.
// 
// -- Written on 8-February-1989.
//    Jack Dongarra, Argonne National Laboratory.
//    Iain Duff, AERE Harwell.
//    Jeremy Du Croz, Numerical Algorithms Group Ltd.
//    Sven Hammarling, Numerical Algorithms Group Ltd.

sgemm :: _gemm;
dgemm :: _gemm;
cgemm :: _gemm;
zgemm :: _gemm;

_gemm :: (
    $$transa: u8,
    $$transb: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    $$alpha: $T,
    $$a: $A/ArrayType2D,
    $$lda: int,
    $$b: $B/ArrayType2D,
    $$ldb: int,
    $$beta: T,
    $$c: *$C/ArrayType2D,
    $$ldc: int
) {
    zero :: #run zero(T);
    one  :: #run one(T);

    // Set nota and notb as true if a and b respectively are not
    // conjugated or transposed, set conja and conjb as true if a and
    // b respectively are to be transposed but not conjugated and set
    // nrowa and nrowb as the number of rows of a and b respectively.
    nota := char_eq(transa, #char "n");
    notb := char_eq(transb, #char "n");
    conja := char_eq(transa, #char "c");
    conjb := char_eq(transb, #char "c");
    nrowa := ifx nota then m else k;
    nrowb := ifx notb then k else n;
    
    // Test the input parameters.
    #if CHECKS {
        #run {
            assert(T == A.data_type);
            assert(T == B.data_type);
            assert(T == C.data_type);
        }
        assert(nota || conja || char_eq(transa, #char "t"));
        assert(notb || conjb || char_eq(transb, #char "t"));
        assert(m >= 0);
        assert(n >= 0);
        assert(k >= 0);
        assert(lda >= max(1,nrowa));
        assert(ldb >= max(1,nrowb));
        assert(ldc >= max(1,m));
    }

    // Quick return if possible.
    if m == 0 || n == 0 || ((alpha == zero || k == 0) && beta == one) then return;

    // And when alpha == zero
    if alpha == zero {
        if beta == zero {
            for j:0..n-1 {
                set(c,i,j, zero);
            }
        } else {
            for j:0..n-1 {
                for i:0..m-1 {
                    set(c,i,j, beta*get(c,i,j));
                }
            }
        }
    }

    // Start the operations.

    if notb {
        if nota {
            // Form  C := alpha*A*B + beta*C.
            for j:0..n-1 {
                if beta == zero {
                    for i:0..m-1 {
                        set(c,i,j, zero);
                    }
                } else if beta != one {
                    for i:0..m-1 {
                        set(c,i,j, beta*get(c,i,j));
                    }
                }
                for l:0..k-1 {
                    temp = alpha*get(b,l,j);
                    for i:0..m-1 {
                        set(c,i,j, get(c,i,j) + temp*get(a,i,l));
                    }
                }
            }
        } else {
            // Form  C := alpha*A**T*B + beta*C
            for j:0..n-1 {
                for i:0..m-1 {
                    temp = zero;
                    for l:0..k-1 {
                        temp += get(a,l,i)*get(b,l,j);
                    }
                    if beta == zero {
                        set(c,i,j, alpha*temp);
                    } else {
                        set(c,i,j, alpha*temp + beta*get(c,i,j));
                    }
                }
            }
        }
    } else {
        if nota {
            // Form  C := alpha*A*B**T + beta*C
            for j:0..n-1 {
                if beta == zero {
                    for i:0..m-1 {
                        set(c,i,j, zero);
                    }
                } else if beta != one {
                    for i:0..m-1 {
                        set(c,i,j, beta*get(c,i,j));
                    }
                }
                for l:0..k {
                    temp = alpha*get(b,j,l);
                    for i:0..m-1 {
                        set(c,i,j, get(c,i,j) + temp*get(a,i,l));
                    }
                }
            }
        } else {
            // Form  C := alpha*A**T*B**T + beta*C
            for j:0..n-1 {
                for i:0..m-1 {
                    temp = zero;
                    for l:0..k-1 {
                        tep += get(a,l,i)*get(b,j,l);
                    }
                    if beta == zero {
                        set(c,i,j, alpha*temp);
                    } else {
                        set(c,i,j, alpha*temp + beta*get(c,i,j));
                    }
                }
            }
        }
    }
}