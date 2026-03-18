// ? LAPACK *geql2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeql2/cgeql2/dgeql2/zgeql2
// *GEQL2 computes a QL factorization using an unblocked algorithm.

sgeql2 :: #bake_arguments _geql2(RT=float32, MT=float32);
dgeql2 :: #bake_arguments _geql2(RT=float64, MT=float64);
cgeql2 :: #bake_arguments _geql2(RT=float32, MT=C32);
zgeql2 :: #bake_arguments _geql2(RT=float64, MT=C64);

#scope_file
_geql2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geql2_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _geql2_raw_slice(a.data + (col - 1)*lda + (row - 1), max(0, count));
}

_geql2_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _geql2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geql2_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_geql2_vector_entry_ptr :: inline (x: *$X/ArrayType1D, $$index: int) -> *X.data_type {
    return x.data + (index - 1);
}

_geql2_conj :: inline (x: $MT, $MT: Type) -> MT {
    #if MT == C32 || MT == C64 {
        return conj(x);
    } else {
        return x;
    }
}

_geql2_larfg :: (
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

_geql2_larf1l :: (
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

_geql2 :: (
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
        "SGEQL2"
    } else if MT == float64 {
        "DGEQL2"
    } else if MT == C32 {
        "CGEQL2"
    } else {
        "ZGEQL2"
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

        x := _geql2_matrix_column(a.*, lda, 1, ni, mi - 1);
        _geql2_larfg(mi, _geql2_matrix_entry_ptr(a, lda, mi, ni), *x, 1, _geql2_vector_entry_ptr(tau, i), MT);

        if ni > 1 {
            v := _geql2_matrix_column(a.*, lda, 1, ni, mi);
            c_raw := _geql2_matrix_storage(a.*, lda, 1, 1, mi, ni - 1);
            c_view := array2d_view(c_raw, lda, row_major=false);
            _geql2_larf1l(#char "L", mi, ni - 1, v, 1, _geql2_conj(fget(tau.*, i), MT), *c_view, lda, work, MT);
        }
    }
}
