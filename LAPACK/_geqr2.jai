// ? LAPACK *geqr2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeqr2/cgeqr2/dgeqr2/zgeqr2
// *GEQR2 computes a QR factorization using an unblocked algorithm.

sgeqr2 :: #bake_arguments _geqr2(RT=float32, MT=float32);
dgeqr2 :: #bake_arguments _geqr2(RT=float64, MT=float64);
cgeqr2 :: #bake_arguments _geqr2(RT=float32, MT=C32);
zgeqr2 :: #bake_arguments _geqr2(RT=float64, MT=C64);

#scope_file
_geqr2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geqr2_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _geqr2_raw_slice(a.data + (col - 1)*lda + (row - 1), max(0, count));
}

_geqr2_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _geqr2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geqr2_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_geqr2_vector_entry_ptr :: inline (x: *$X/ArrayType1D, $$index: int) -> *X.data_type {
    return x.data + (index - 1);
}

_geqr2_conj :: inline (x: $MT, $MT: Type) -> MT {
    #if MT == C32 || MT == C64 {
        return conjugate(x);
    } else {
        return x;
    }
}

_geqr2_larfg :: (
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

_geqr2_larf1f :: (
    $$side: u8, $$m: int, $$n: int, v: $V/ArrayType1D, $$incv: int, tau: $MT,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT && C.data_type == MT && WORK.data_type == MT;
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

_geqr2 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGEQR2"
    } else if MT == float64 {
        "DGEQR2"
    } else if MT == C32 {
        "CGEQR2"
    } else {
        "ZGEQR2"
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
        tail := _geqr2_matrix_column(a.*, lda, min(i + 1, m), i, max(0, m - i));
        _geqr2_larfg(m - i + 1, _geqr2_matrix_entry_ptr(a, lda, i, i), *tail, 1, _geqr2_vector_entry_ptr(tau, i), MT);
        if i < n {
            v := _geqr2_matrix_column(a.*, lda, i, i, m - i + 1);
            c_raw := _geqr2_matrix_storage(a.*, lda, i, i + 1, m - i + 1, n - i);
            c_view := array2d_view(c_raw, lda, row_major=false);
            _geqr2_larf1f(#char "L", m - i + 1, n - i, v, 1, _geqr2_conj(fget(tau.*, i), MT), *c_view, lda, work, MT);
        }
    }
}
