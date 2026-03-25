// ? LAPACK *larcm()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clarcm/zlarcm
// *LARCM multiplies a real matrix by a complex matrix.

clarcm :: #bake_arguments _larcm(RT=float32, MT=C32);
zlarcm :: #bake_arguments _larcm(RT=float64, MT=C64);

_larcm_gemm :: (
    $$m: int, $$n: int, $$k: int, $$alpha: $RT, a: $A/ArrayType2D, $$lda: int,
    b: $B/ArrayType2D, $$ldb: int, $$beta: $RT, c: *$C/ArrayType2D, $$ldc: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && A.data_type == RT && B.data_type == RT && C.data_type == RT;
} {
    #if RT == float32 {
        sgemm(#char "N", #char "N", m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        dgemm(#char "N", #char "N", m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_larcm :: (
    $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int,
    c: *$C/ArrayType2D, $$ldc: int, rwork: *$RWORK/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == RT && B.data_type == MT && C.data_type == MT && RWORK.data_type == RT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    l, i, j : int = ---;

    if m == 0 || n == 0 then return;

    for j:1..n {
        for i:1..m {
            fset(rwork, (j - 1)*m + i, real(fget(b, i, j)));
        }
    }

    l = m*n + 1;
    rleft := array2d_view(rwork.data[0..], m, row_major=false);
    rright := array2d_view(rwork.data[l - 1..], m, row_major=false);
    _larcm_gemm(m, n, m, one_r, a, lda, rleft, m, zero_r, *rright, m, RT);
    for j:1..n {
        for i:1..m {
            fset(c, i, j, cpx(fget(rwork.*, l + (j - 1)*m + i - 1), zero(RT)));
        }
    }

    for j:1..n {
        for i:1..m {
            fset(rwork, (j - 1)*m + i, imag(fget(b, i, j)));
        }
    }
    _larcm_gemm(m, n, m, one_r, a, lda, rleft, m, zero_r, *rright, m, RT);
    for j:1..n {
        for i:1..m {
            fset(c, i, j, cpx(real(fget(c.*, i, j)), fget(rwork.*, l + (j - 1)*m + i - 1)));
        }
    }
}
