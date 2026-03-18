// ? LAPACK *gerq2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgerq2/cgerq2/dgerq2/zgerq2
// *GERQ2 computes an RQ factorization using an unblocked algorithm.

sgerq2 :: #bake_arguments _gerq2(RT=float32, MT=float32);
dgerq2 :: #bake_arguments _gerq2(RT=float64, MT=float64);
cgerq2 :: #bake_arguments _gerq2(RT=float32, MT=C32);
zgerq2 :: #bake_arguments _gerq2(RT=float64, MT=C64);

#scope_file
_gerq2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gerq2_matrix_row :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _gerq2_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_gerq2_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gerq2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gerq2_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_gerq2_vector_entry_ptr :: inline (x: *$X/ArrayType1D, $$index: int) -> *X.data_type {
    return x.data + (index - 1);
}

_gerq2_larfg :: (
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

_gerq2_larf1l :: (
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
        slarf1l(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == float64 {
        dlarf1l(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == C32 {
        clarf1l(side, m, n, v, incv, tau, c, ldc, work);
    } else {
        zlarf1l(side, m, n, v, incv, tau, c, ldc, work);
    }
}

_gerq2_lacgv :: (
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

_gerq2 :: (
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
        "SGERQ2"
    } else if MT == float64 {
        "DGERQ2"
    } else if MT == C32 {
        "CGERQ2"
    } else {
        "ZGERQ2"
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
    for <i:1..k {
        mi := m - k + i;
        ni := n - k + i;

        row := _gerq2_matrix_row(a.*, lda, mi, 1, ni);
        #if MT == C32 || MT == C64 {
            _gerq2_lacgv(ni, *row, lda, MT);
        }

        tail := _gerq2_matrix_row(a.*, lda, mi, 1, max(0, ni - 1));
        _gerq2_larfg(ni, _gerq2_matrix_entry_ptr(a, lda, mi, ni), *tail, lda, _gerq2_vector_entry_ptr(tau, i), MT);

        if mi > 1 {
            v := _gerq2_matrix_row(a.*, lda, mi, 1, ni);
            c_raw := _gerq2_matrix_storage(a.*, lda, 1, 1, mi - 1, ni);
            c_view := array2d_view(c_raw, lda, row_major=false);
            _gerq2_larf1l(#char "R", mi - 1, ni, v, lda, fget(tau.*, i), *c_view, lda, work, MT);
        }

        #if MT == C32 || MT == C64 {
            row_restore := _gerq2_matrix_row(a.*, lda, mi, 1, max(0, ni - 1));
            _gerq2_lacgv(ni - 1, *row_restore, lda, MT);
        }
    }
}
