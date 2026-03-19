// ? LAPACK *getrf()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgetrf/cgetrf/dgetrf/zgetrf
// *GETRF computes the LU factorization of a general M-by-N matrix
// using partial pivoting with row interchanges.

sgetrf :: #bake_arguments _getrf(MT=float32);
dgetrf :: #bake_arguments _getrf(MT=float64);
cgetrf :: #bake_arguments _getrf(MT=C32);
zgetrf :: #bake_arguments _getrf(MT=C64);

#scope_file
_getrf_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    view: [] A.data_type;
    view.data = a.data + (col - 1)*lda + (row - 1);
    view.count = count;
    return view;
}

_getrf_laswp :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$k1: int,
    $$k2: int,
    ipiv: $IPIV/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int;
} {
    #if MT == float32 {
        slaswp(n, a, lda, k1, k2, ipiv, incx);
    } else if MT == float64 {
        dlaswp(n, a, lda, k1, k2, ipiv, incx);
    } else if MT == C32 {
        claswp(n, a, lda, k1, k2, ipiv, incx);
    } else {
        zlaswp(n, a, lda, k1, k2, ipiv, incx);
    }
}

_getrf_trsm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == float32 {
        strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == float64 {
        dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_getrf_gemm :: (
    $$transa: u8,
    $$transb: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: MT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
} {
    #if MT == float32 {
        sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == float64 {
        dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_getrf_getrf2 :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int;
} {
    #if MT == float32 {
        sgetrf2(m, n, a, lda, ipiv, info);
    } else if MT == float64 {
        dgetrf2(m, n, a, lda, ipiv, info);
    } else if MT == C32 {
        cgetrf2(m, n, a, lda, ipiv, info);
    } else {
        zgetrf2(m, n, a, lda, ipiv, info);
    }
}

_getrf :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int;
} {
    one_t :: #run one(MT);
    routine_name :: #run if MT == float32 {
        "SGETRF"
    } else if MT == float64 {
        "DGETRF"
    } else if MT == C32 {
        "CGETRF"
    } else {
        "ZGETRF"
    };

    i, iinfo, j, jb, nb : int = ---;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if m == 0 || n == 0 then return;

    nb = ilaenv(1, routine_name, #char " ", m, n, -1, -1);
    if nb <= 1 || nb >= min(m, n) {
        _getrf_getrf2(m, n, a, lda, ipiv, info, MT);
    } else {
        j = 1;
        while j <= min(m, n) {
            jb = min(min(m, n) - j + 1, nb);

            panel_raw := _getrf_matrix_storage(a.*, lda, j, j, m - j + 1, jb);
            panel_view := array2d_view(panel_raw, lda, row_major=false);
            ipiv_panel: [] int = ipiv.*;
            ipiv_panel.data += j - 1;
            ipiv_panel.count = min(m, n) - j + 1;
            iinfo = 0;
            _getrf_getrf2(m - j + 1, jb, *panel_view, lda, *ipiv_panel, *iinfo, MT);

            if info.* == 0 && iinfo > 0 {
                info.* = iinfo + j - 1;
            }
            for i:j..min(m, j + jb - 1) {
                fset(ipiv, i, j - 1 + fget(ipiv.*, i));
            }

            if j > 1 {
                left_raw := _getrf_matrix_storage(a.*, lda, 1, 1, m, j - 1);
                left_view := array2d_view(left_raw, lda, row_major=false);
                _getrf_laswp(j - 1, *left_view, lda, j, j + jb - 1, ipiv.*, 1, MT);
            }

            if j + jb <= n {
                right_raw := _getrf_matrix_storage(a.*, lda, 1, j + jb, m, n - j - jb + 1);
                right_view := array2d_view(right_raw, lda, row_major=false);
                _getrf_laswp(n - j - jb + 1, *right_view, lda, j, j + jb - 1, ipiv.*, 1, MT);

                a11_raw := _getrf_matrix_storage(a.*, lda, j, j, jb, jb);
                a11_view := array2d_view(a11_raw, lda, row_major=false);
                a12_raw := _getrf_matrix_storage(a.*, lda, j, j + jb, jb, n - j - jb + 1);
                a12_view := array2d_view(a12_raw, lda, row_major=false);
                _getrf_trsm(#char "L", #char "L", #char "N", #char "U", jb, n - j - jb + 1, one_t, *a11_view, lda, *a12_view, lda, MT);

                if j + jb <= m {
                    a21_raw := _getrf_matrix_storage(a.*, lda, j + jb, j, m - j - jb + 1, jb);
                    a21_view := array2d_view(a21_raw, lda, row_major=false);
                    a22_raw := _getrf_matrix_storage(a.*, lda, j + jb, j + jb, m - j - jb + 1, n - j - jb + 1);
                    a22_view := array2d_view(a22_raw, lda, row_major=false);
                    _getrf_gemm(#char "N", #char "N", m - j - jb + 1, n - j - jb + 1, jb, -one_t, *a21_view, lda, *a12_view, lda, one_t, *a22_view, lda, MT);
                }
            }
            j += nb;
        }
    }
}
