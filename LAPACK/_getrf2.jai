// ? LAPACK *getrf2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgetrf2/cgetrf2/dgetrf2/zgetrf2
// *GETRF2 computes the recursive LU factorization of a general
// M-by-N matrix using partial pivoting with row interchanges.

sgetrf2 :: #bake_arguments _getrf2(RT=float32, MT=float32);
dgetrf2 :: #bake_arguments _getrf2(RT=float64, MT=float64);
cgetrf2 :: #bake_arguments _getrf2(RT=float32, MT=C32);
zgetrf2 :: #bake_arguments _getrf2(RT=float64, MT=C64);

#scope_file
_getrf2_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    view: [] A.data_type;
    view.data = a.data + (col - 1)*lda + (row - 1);
    view.count = count;
    return view;
}

_getrf2_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_getrf2_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $MT: Type
) -> int #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return isamax(n, x, 1) + 1;
    } else if MT == float64 {
        return idamax(n, x, 1) + 1;
    } else if MT == C32 {
        return icamax(n, x, 1) + 1;
    } else {
        return izamax(n, x, 1) + 1;
    }
}

_getrf2_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_getrf2_laswp :: (
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

_getrf2_trsm :: (
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

_getrf2_gemm :: (
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

_getrf2_scal_real :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        sscal(n, alpha, x, 1);
    } else {
        dscal(n, alpha, x, 1);
    }
}

_getrf2_scal_complex :: (
    $$n: int,
    alpha: $MT,
    x: *$X/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        cscal(n, alpha, x, 1);
    } else {
        zscal(n, alpha, x, 1);
    }
}

_getrf2_swap_values :: inline (x: *$T, y: *T) {
    temp := x.*;
    x.* = y.*;
    y.* = temp;
}

_getrf2 :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           IPIV.data_type == int;
} {
    one_t  :: #run ifx MT == RT then ncast(MT, one(RT)) else cpx(one(RT), zero(RT));
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == float32 {
        "SGETRF2"
    } else if MT == float64 {
        "DGETRF2"
    } else if MT == C32 {
        "CGETRF2"
    } else {
        "ZGETRF2"
    };

    i, iinfo, n1, n2 : int = ---;
    sfmin : RT = ---;

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

    if m == 1 {
        fset(ipiv, 1, 1);
        if fget(a.*, 1, 1) == zero_t {
            info.* = 1;
        }
    } else if n == 1 {
        sfmin = _getrf2_lamch_safe_minimum(RT);
        col_view := _getrf2_matrix_storage(a.*, lda, 1, 1, m, 1);
        i = _getrf2_iamax_1based(m, col_view, MT);
        fset(ipiv, 1, i);
        if fget(a.*, i, 1) != zero_t {
            if i != 1 {
                _getrf2_swap_values(a.data, a.data + (i - 1));
            }
            if m > 1 {
                tail_view: [] MT;
                tail_view.data = a.data + 1;
                tail_view.count = m - 1;
                #if MT == RT {
                    if abs(fget(a.*, 1, 1)) >= sfmin {
                        _getrf2_scal_real(m - 1, one(RT) / fget(a.*, 1, 1), *tail_view, RT);
                    } else {
                        for i:1..m - 1 {
                            fset(a, i + 1, 1, fget(a.*, i + 1, 1) / fget(a.*, 1, 1));
                        }
                    }
                } else {
                    _getrf2_scal_complex(m - 1, one_t / fget(a.*, 1, 1), *tail_view, MT);
                }
            }
        } else {
            info.* = 1;
        }
    } else {
        n1 = min(m, n) / 2;
        n2 = n - n1;

        left_info : int = 0;
        right_info : int = 0;

        left_raw := _getrf2_matrix_storage(a.*, lda, 1, 1, m, n1);
        left_view := array2d_view(left_raw, lda, row_major=false);
        _getrf2(m, n1, *left_view, lda, ipiv, *left_info, RT, MT);
        iinfo = left_info;
        if info.* == 0 && iinfo > 0 {
            info.* = iinfo;
        }

        right_raw := _getrf2_matrix_storage(a.*, lda, 1, n1 + 1, m, n2);
        right_view := array2d_view(right_raw, lda, row_major=false);
        _getrf2_laswp(n2, *right_view, lda, 1, n1, ipiv.*, 1, MT);

        a11_raw := _getrf2_matrix_storage(a.*, lda, 1, 1, n1, n1);
        a11_view := array2d_view(a11_raw, lda, row_major=false);
        a12_raw := _getrf2_matrix_storage(a.*, lda, 1, n1 + 1, n1, n2);
        a12_view := array2d_view(a12_raw, lda, row_major=false);
        _getrf2_trsm(#char "L", #char "L", #char "N", #char "U", n1, n2, one_t, *a11_view, lda, *a12_view, lda, MT);

        if m > n1 {
            a21_raw := _getrf2_matrix_storage(a.*, lda, n1 + 1, 1, m - n1, n1);
            a21_view := array2d_view(a21_raw, lda, row_major=false);
            a22_raw := _getrf2_matrix_storage(a.*, lda, n1 + 1, n1 + 1, m - n1, n2);
            a22_view := array2d_view(a22_raw, lda, row_major=false);
            _getrf2_gemm(#char "N", #char "N", m - n1, n2, n1, -one_t, *a21_view, lda, *a12_view, lda, one_t, *a22_view, lda, MT);

            ipiv_tail := _getrf2_subvector(ipiv.*, n1 + 1, min(m, n) - n1);
            _getrf2(m - n1, n2, *a22_view, lda, *ipiv_tail, *right_info, RT, MT);
            iinfo = right_info;
            if info.* == 0 && iinfo > 0 {
                info.* = iinfo + n1;
            }
            for i:n1 + 1..min(m, n) {
                fset(ipiv, i, fget(ipiv.*, i) + n1);
            }

            a21_apply_raw := _getrf2_matrix_storage(a.*, lda, 1, 1, m, n1);
            a21_apply_view := array2d_view(a21_apply_raw, lda, row_major=false);
            _getrf2_laswp(n1, *a21_apply_view, lda, n1 + 1, min(m, n), ipiv.*, 1, MT);
        }
    }
}
