// ? LAPACK *lacrm()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clacrm/zlacrm
// *LACRM multiplies a complex matrix by a square real matrix.

clacrm :: #bake_arguments _lacrm(RT=float32, MT=C32);
zlacrm :: #bake_arguments _lacrm(RT=float64, MT=C64);

#scope_file
_lacrm_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lacrm_gemm :: inline (
    $$m: int,
    $$n: int,
    $$k: int,
    alpha: $RT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: RT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           B.data_type == RT &&
           C.data_type == RT;
} {
    #if RT == float32 {
        sgemm(#char "N", #char "N", m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        dgemm(#char "N", #char "N", m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_lacrm :: (
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    rwork: *$RWORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == RT &&
           C.data_type == MT &&
           RWORK.data_type == RT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    i, j, l : int = ---;

    if m == 0 || n == 0 then return;

    work_a := array2d_view(_lacrm_raw_slice(rwork.data, m*n), m, row_major=false);
    l = m*n + 1;
    work_c := array2d_view(_lacrm_raw_slice(rwork.data + (l - 1), m*n), m, row_major=false);

    for j:1..n {
        for i:1..m {
            fset(*work_a, i, j, real(fget(a, i, j)));
        }
    }
    _lacrm_gemm(m, n, n, one_r, work_a, m, b, ldb, zero_r, *work_c, m, RT);
    for j:1..n {
        for i:1..m {
            fset(c, i, j, ncast(MT, fget(*work_c, i, j)));
        }
    }

    for j:1..n {
        for i:1..m {
            fset(*work_a, i, j, imag(fget(a, i, j)));
        }
    }
    _lacrm_gemm(m, n, n, one_r, work_a, m, b, ldb, zero_r, *work_c, m, RT);
    for j:1..n {
        for i:1..m {
            fset(c, i, j, cpx(real(fget(c.*, i, j)), fget(*work_c, i, j)));
        }
    }
}
