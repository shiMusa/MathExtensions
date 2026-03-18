// ? LAPACK *gelq2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgelq2/cgelq2/dgelq2/zgelq2
// *GELQ2 computes an LQ factorization using an unblocked algorithm.

sgelq2 :: #bake_arguments _gelq2(RT=float32, MT=float32);
dgelq2 :: #bake_arguments _gelq2(RT=float64, MT=float64);
cgelq2 :: #bake_arguments _gelq2(RT=float32, MT=C32);
zgelq2 :: #bake_arguments _gelq2(RT=float64, MT=C64);

#scope_file
_gelq2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gelq2_matrix_row :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _gelq2_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_gelq2_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gelq2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gelq2_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_gelq2_vector_entry_ptr :: inline (x: *$X/ArrayType1D, $$index: int) -> *X.data_type {
    return x.data + (index - 1);
}

_gelq2_larfg :: (
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

_gelq2_larf1f :: (
    $$side: u8,
    $$m: int,
    $$n: int,
    v: $V/ArrayType1D,
    $$incv: int,
    tau: $MT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        slarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == float64 {
        dlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == C32 {
        clarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else {
        zlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    }
}

_gelq2_lacgv :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        clacgv(n, x, incx);
    } else {
        zlacgv(n, x, incx);
    }
}

_gelq2 :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGELQ2"
    } else if MT == float64 {
        "DGELQ2"
    } else if MT == C32 {
        "CGELQ2"
    } else {
        "ZGELQ2"
    };

    i, k : int = ---;

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

    k = min(m, n);
    for i:1..k {
        row := _gelq2_matrix_row(a.*, lda, i, i, n - i + 1);
        #if MT == C32 || MT == C64 {
            _gelq2_lacgv(n - i + 1, *row, lda, MT);
        }

        tail := _gelq2_matrix_row(a.*, lda, i, min(i + 1, n), max(0, n - i));
        _gelq2_larfg(n - i + 1, _gelq2_matrix_entry_ptr(a, lda, i, i), *tail, lda, _gelq2_vector_entry_ptr(tau, i), MT);

        if i < m {
            v := _gelq2_matrix_row(a.*, lda, i, i, n - i + 1);
            c_raw := _gelq2_matrix_storage(a.*, lda, i + 1, i, m - i, n - i + 1);
            c_view := array2d_view(c_raw, lda, row_major=false);
            _gelq2_larf1f(#char "R", m - i, n - i + 1, v, lda, fget(tau.*, i), *c_view, lda, work, MT);
        }

        #if MT == C32 || MT == C64 {
            row_restore := _gelq2_matrix_row(a.*, lda, i, i, n - i + 1);
            _gelq2_lacgv(n - i + 1, *row_restore, lda, MT);
        }
    }
}
