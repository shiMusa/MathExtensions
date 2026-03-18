// ? LAPACK *geqrt3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeqrt3/cgeqrt3/dgeqrt3/zgeqrt3
// *GEQRT3 recursively computes a compact-WY QR factorization.

sgeqrt3 :: #bake_arguments _geqrt3(RT=float32, MT=float32);
dgeqrt3 :: #bake_arguments _geqrt3(RT=float64, MT=float64);
cgeqrt3 :: #bake_arguments _geqrt3(RT=float32, MT=C32);
zgeqrt3 :: #bake_arguments _geqrt3(RT=float64, MT=C64);

#scope_file
_geqrt3_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geqrt3_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _geqrt3_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geqrt3_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_geqrt3_larfg :: (
    $$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *$MT, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == float32 {
        slarfg(n, alpha, x, incx, tau);
    } else if MT == float64 {
        dlarfg(n, alpha, x, incx, tau);
    } else if MT == C32 {
        clarfg(n, alpha, x, incx, tau);
    } else {
        zlarfg(n, alpha, x, incx, tau);
    }
}

_geqrt3_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int,
    beta: $MT, c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT && C.data_type == MT;
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

_geqrt3_trmm :: (
    $$side: u8, $$uplo: u8, $$transa: u8, $$diag: u8, $$m: int, $$n: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT;
} {
    #if MT == float32 {
        strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == float64 {
        dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == C32 {
        ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_geqrt3 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, t: *$T/ArrayType2D, $$ldt: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && T.data_type == MT;
} {
    one_t :: #run if MT == float32 || MT == float64 { one(MT) } else { cpx(one(RT), zero(RT)) };
    trans_t :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };

    routine_name :: #run if MT == float32 {
        "SGEQRT3"
    } else if MT == float64 {
        "DGEQRT3"
    } else if MT == C32 {
        "CGEQRT3"
    } else {
        "ZGEQRT3"
    };

    i, i1, j, j1, n1, n2, iinfo : int = ---;

    info.* = 0;
    if n < 0 {
        info.* = -2;
    } else if m < n {
        info.* = -1;
    } else if lda < max(1, m) {
        info.* = -4;
    } else if ldt < max(1, n) {
        info.* = -6;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 1 {
        tail := _geqrt3_raw_slice(_geqrt3_matrix_entry_ptr(a, lda, min(2, m), 1), max(0, m - 1));
        _geqrt3_larfg(m, _geqrt3_matrix_entry_ptr(a, lda, 1, 1), *tail, 1, _geqrt3_matrix_entry_ptr(t, ldt, 1, 1), MT);
        return;
    }

    n1 = n / 2;
    n2 = n - n1;
    j1 = min(n1 + 1, n);
    i1 = min(n + 1, m);

    _geqrt3(m, n1, a, lda, t, ldt, *iinfo, RT, MT);

    for j:1..n2 {
        for i:1..n1 {
            fset(t, i, j + n1, fget(a.*, i, j + n1));
        }
    }

    a11_raw := _geqrt3_matrix_storage(a.*, lda, 1, 1, n1, n1);
    a11 := array2d_view(a11_raw, lda, row_major=false);
    t12_raw := _geqrt3_matrix_storage(t.*, ldt, 1, j1, n1, n2);
    t12 := array2d_view(t12_raw, ldt, row_major=false);
    _geqrt3_trmm(#char "L", #char "L", trans_t, #char "U", n1, n2, one_t, a11.*, lda, *t12, ldt, MT);

    a21_raw := _geqrt3_matrix_storage(a.*, lda, j1, 1, m - n1, n1);
    a21 := array2d_view(a21_raw, lda, row_major=false);
    a22_raw := _geqrt3_matrix_storage(a.*, lda, j1, j1, m - n1, n2);
    a22 := array2d_view(a22_raw, lda, row_major=false);
    _geqrt3_gemm(trans_t, #char "N", n1, n2, m - n1, one_t, a21.*, lda, a22.*, lda, one_t, *t12, ldt, MT);

    t11_raw := _geqrt3_matrix_storage(t.*, ldt, 1, 1, n1, n1);
    t11 := array2d_view(t11_raw, ldt, row_major=false);
    _geqrt3_trmm(#char "L", #char "U", trans_t, #char "N", n1, n2, one_t, t11.*, ldt, *t12, ldt, MT);

    minus_one_t := -one_t;
    _geqrt3_gemm(#char "N", #char "N", m - n1, n2, n1, minus_one_t, a21.*, lda, t12.*, ldt, one_t, *a22, lda, MT);

    _geqrt3_trmm(#char "L", #char "L", #char "N", #char "U", n1, n2, one_t, a11.*, lda, *t12, ldt, MT);

    for j:1..n2 {
        for i:1..n1 {
            fset(a, i, j + n1, fget(a.*, i, j + n1) - fget(t.*, i, j + n1));
        }
    }

    a22_fact_raw := _geqrt3_matrix_storage(a.*, lda, j1, j1, m - n1, n2);
    a22_fact := array2d_view(a22_fact_raw, lda, row_major=false);
    t22_raw := _geqrt3_matrix_storage(t.*, ldt, j1, j1, n2, n2);
    t22 := array2d_view(t22_raw, ldt, row_major=false);
    _geqrt3(m - n1, n2, *a22_fact, lda, *t22, ldt, *iinfo, RT, MT);

    for i:1..n1 {
        for j:1..n2 {
            #if MT == C32 || MT == C64 {
                fset(t, i, j + n1, conjugate(fget(a.*, j + n1, i)));
            } else {
                fset(t, i, j + n1, fget(a.*, j + n1, i));
            }
        }
    }

    a22_now_raw := _geqrt3_matrix_storage(a.*, lda, j1, j1, n2, n2);
    a22_now := array2d_view(a22_now_raw, lda, row_major=false);
    _geqrt3_trmm(#char "R", #char "L", #char "N", #char "U", n1, n2, one_t, a22_now.*, lda, *t12, ldt, MT);

    if m - n > 0 {
        a31_raw := _geqrt3_matrix_storage(a.*, lda, i1, 1, m - n, n1);
        a31 := array2d_view(a31_raw, lda, row_major=false);
        a32_raw := _geqrt3_matrix_storage(a.*, lda, i1, j1, m - n, n2);
        a32 := array2d_view(a32_raw, lda, row_major=false);
        _geqrt3_gemm(trans_t, #char "N", n1, n2, m - n, one_t, a31.*, lda, a32.*, lda, one_t, *t12, ldt, MT);
    }

    _geqrt3_trmm(#char "L", #char "U", #char "N", #char "N", n1, n2, -one_t, t11.*, ldt, *t12, ldt, MT);
    _geqrt3_trmm(#char "R", #char "U", #char "N", #char "N", n1, n2, one_t, t22.*, ldt, *t12, ldt, MT);
}
