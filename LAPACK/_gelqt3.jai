// ? LAPACK *gelqt3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgelqt3/cgelqt3/dgelqt3/zgelqt3
// *GELQT3 recursively computes a compact-WY LQ factorization.

sgelqt3 :: #bake_arguments _gelqt3(RT=float32, MT=float32);
dgelqt3 :: #bake_arguments _gelqt3(RT=float64, MT=float64);
cgelqt3 :: #bake_arguments _gelqt3(RT=float32, MT=C32);
zgelqt3 :: #bake_arguments _gelqt3(RT=float64, MT=C64);

#scope_file
_gelqt3_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gelqt3_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gelqt3_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gelqt3_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_gelqt3_larfg :: (
    $$n: int,
    alpha: *$MT,
    x: *$X/ArrayType1D,
    $$incx: int,
    tau: *$MT,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
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

_gelqt3_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int,
    beta: MT, c: *$C/ArrayType2D, $$ldc: int, $MT: Type
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

_gelqt3_trmm :: (
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

_gelqt3_conj_tau :: inline (x: $MT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == C32 || MT == C64 {
        return conjugate(x);
    } else {
        return x;
    }
}

_gelqt3 :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           T.data_type == MT;
} {
    one_t :: #run if MT == float32 || MT == float64 { one(MT) } else { cpx(one(RT), zero(RT)) };
    zero_t :: #run zero(MT);
    trans_t :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
    routine_name :: #run if MT == float32 {
        "SGELQT3"
    } else if MT == float64 {
        "DGELQT3"
    } else if MT == C32 {
        "CGELQT3"
    } else {
        "ZGELQT3"
    };

    i, i1, j, j1, m1, m2, iinfo : int = ---;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < m {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    } else if ldt < max(1, m) {
        info.* = -6;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if m == 1 {
        tail := _gelqt3_raw_slice(_gelqt3_matrix_entry_ptr(a, lda, 1, min(2, n)), max(0, n - 1));
        _gelqt3_larfg(n, _gelqt3_matrix_entry_ptr(a, lda, 1, 1), *tail, lda, _gelqt3_matrix_entry_ptr(t, ldt, 1, 1), MT);
        #if MT == C32 || MT == C64 {
            fset(t, 1, 1, conjugate(fget(t.*, 1, 1)));
        }
        return;
    }

    m1 = m / 2;
    m2 = m - m1;
    i1 = min(m1 + 1, m);
    j1 = min(m + 1, n);

    _gelqt3(m1, n, a, lda, t, ldt, *iinfo, RT, MT);

    for i:1..m2 {
        for j:1..m1 {
            fset(t, i + m1, j, fget(a.*, i + m1, j));
        }
    }

    a11_raw := _gelqt3_matrix_storage(a.*, lda, 1, 1, m1, m1);
    a11 := array2d_view(a11_raw, lda, row_major=false);
    t21_raw := _gelqt3_matrix_storage(t.*, ldt, i1, 1, m2, m1);
    t21 := array2d_view(t21_raw, ldt, row_major=false);
    _gelqt3_trmm(#char "R", #char "U", trans_t, #char "U", m2, m1, one_t, a11.*, lda, *t21, ldt, MT);

    a22_raw := _gelqt3_matrix_storage(a.*, lda, i1, i1, m2, n - m1);
    a22 := array2d_view(a22_raw, lda, row_major=false);
    a12_raw := _gelqt3_matrix_storage(a.*, lda, 1, i1, m1, n - m1);
    a12 := array2d_view(a12_raw, lda, row_major=false);
    _gelqt3_gemm(#char "N", trans_t, m2, m1, n - m1, one_t, a22.*, lda, a12.*, lda, one_t, *t21, ldt, MT);

    t11_raw := _gelqt3_matrix_storage(t.*, ldt, 1, 1, m1, m1);
    t11 := array2d_view(t11_raw, ldt, row_major=false);
    _gelqt3_trmm(#char "R", #char "U", #char "N", #char "N", m2, m1, one_t, t11.*, ldt, *t21, ldt, MT);

    minus_one_t := -one_t;
    _gelqt3_gemm(#char "N", #char "N", m2, n - m1, m1, minus_one_t, t21.*, ldt, a12.*, lda, one_t, *a22, lda, MT);

    _gelqt3_trmm(#char "R", #char "U", #char "N", #char "U", m2, m1, one_t, a11.*, lda, *t21, ldt, MT);

    for i:1..m2 {
        for j:1..m1 {
            fset(a, i + m1, j, fget(a.*, i + m1, j) - fget(t.*, i + m1, j));
            fset(t, i + m1, j, zero_t);
        }
    }

    a22_fact_raw := _gelqt3_matrix_storage(a.*, lda, i1, i1, m2, n - m1);
    a22_fact := array2d_view(a22_fact_raw, lda, row_major=false);
    t22_raw := _gelqt3_matrix_storage(t.*, ldt, i1, i1, m2, m2);
    t22 := array2d_view(t22_raw, ldt, row_major=false);
    _gelqt3(m2, n - m1, *a22_fact, lda, *t22, ldt, *iinfo, RT, MT);

    for i:1..m2 {
        for j:1..m1 {
            fset(t, j, i + m1, fget(a.*, j, i + m1));
        }
    }

    t12_raw := _gelqt3_matrix_storage(t.*, ldt, 1, i1, m1, m2);
    t12 := array2d_view(t12_raw, ldt, row_major=false);
    a22_now_raw := _gelqt3_matrix_storage(a.*, lda, i1, i1, m2, m2);
    a22_now := array2d_view(a22_now_raw, lda, row_major=false);
    _gelqt3_trmm(#char "R", #char "U", trans_t, #char "U", m1, m2, one_t, a22_now.*, lda, *t12, ldt, MT);

    if n - m > 0 {
        a13_raw := _gelqt3_matrix_storage(a.*, lda, 1, j1, m1, n - m);
        a13 := array2d_view(a13_raw, lda, row_major=false);
        a23_raw := _gelqt3_matrix_storage(a.*, lda, i1, j1, m2, n - m);
        a23 := array2d_view(a23_raw, lda, row_major=false);
        _gelqt3_gemm(#char "N", trans_t, m1, m2, n - m, one_t, a13.*, lda, a23.*, lda, one_t, *t12, ldt, MT);
    }

    _gelqt3_trmm(#char "L", #char "U", #char "N", #char "N", m1, m2, -one_t, t11.*, ldt, *t12, ldt, MT);
    _gelqt3_trmm(#char "R", #char "U", #char "N", #char "N", m1, m2, one_t, t22.*, ldt, *t12, ldt, MT);
}
