// ? LAPACK *gehd2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgehd2/cgehd2/dgehd2/zgehd2
// *GEHD2 reduces a general matrix to upper Hessenberg form
// using an unblocked algorithm.

sgehd2 :: #bake_arguments _gehd2(RT=float32, MT=float32);
dgehd2 :: #bake_arguments _gehd2(RT=float64, MT=float64);
cgehd2 :: #bake_arguments _gehd2(RT=float32, MT=C32);
zgehd2 :: #bake_arguments _gehd2(RT=float64, MT=C64);

#scope_file
_gehd2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gehd2_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _gehd2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gehd2_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gehd2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gehd2_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_gehd2_vector_entry_ptr :: inline (x: *$X/ArrayType1D, $$index: int) -> *X.data_type {
    return x.data + (index - 1);
}

_gehd2_larfg :: (
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

_gehd2_larf1f :: (
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

_gehd2 :: (
    $$n: int,
    $$ilo: int,
    $$ihi: int,
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
        "SGEHD2"
    } else if MT == float64 {
        "DGEHD2"
    } else if MT == C32 {
        "CGEHD2"
    } else {
        "ZGEHD2"
    };

    i : int = ---;

    info.* = 0;
    if n < 0 {
        info.* = -1;
    } else if ilo < 1 || ilo > max(1, n) {
        info.* = -2;
    } else if ihi < min(ilo, n) || ihi > n {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    for i:ilo..ihi-1 {
        col_tail := _gehd2_matrix_column(a.*, lda, min(i + 2, n), i, max(0, ihi - i - 1));
        _gehd2_larfg(ihi - i, _gehd2_matrix_entry_ptr(a, lda, i + 1, i), *col_tail, 1, _gehd2_vector_entry_ptr(tau, i), MT);

        v := _gehd2_matrix_column(a.*, lda, i + 1, i, ihi - i);

        right_raw := _gehd2_matrix_storage(a.*, lda, 1, i + 1, ihi, ihi - i);
        right_view := array2d_view(right_raw, lda, row_major=false);
        _gehd2_larf1f(#char "R", ihi, ihi - i, v, 1, fget(tau.*, i), *right_view, lda, work, MT);

        left_raw := _gehd2_matrix_storage(a.*, lda, i + 1, i + 1, ihi - i, n - i);
        left_view := array2d_view(left_raw, lda, row_major=false);
        #if MT == C32 || MT == C64 {
            _gehd2_larf1f(#char "L", ihi - i, n - i, v, 1, conjugate(fget(tau.*, i)), *left_view, lda, work, MT);
        } else {
            _gehd2_larf1f(#char "L", ihi - i, n - i, v, 1, fget(tau.*, i), *left_view, lda, work, MT);
        }
    }
}
