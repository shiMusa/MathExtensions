// ? LAPACK *geqrt2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeqrt2/cgeqrt2/dgeqrt2/zgeqrt2
// *GEQRT2 computes a compact-WY QR factorization.

sgeqrt2 :: #bake_arguments _geqrt2(RT=float32, MT=float32);
dgeqrt2 :: #bake_arguments _geqrt2(RT=float64, MT=float64);
cgeqrt2 :: #bake_arguments _geqrt2(RT=float32, MT=C32);
zgeqrt2 :: #bake_arguments _geqrt2(RT=float64, MT=C64);

#scope_file
_geqrt2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geqrt2_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _geqrt2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geqrt2_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _geqrt2_raw_slice(a.data + (col - 1)*lda + (row - 1), max(0, count));
}

_geqrt2_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_geqrt2_larfg :: (
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

_geqrt2_gemv :: (
    $$trans: u8, $$m: int, $$n: int, alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    x: $X/ArrayType1D, $$incx: int, beta: $MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == float64 {
        dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == C32 {
        cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_geqrt2_ger_update :: (
    $$m: int, $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int,
    a: *$A/ArrayType2D, $$lda: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT && Y.data_type == MT && A.data_type == MT;
} {
    #if MT == float32 {
        sger(m, n, alpha, x, incx, y, incy, a, lda);
    } else if MT == float64 {
        dger(m, n, alpha, x, incx, y, incy, a, lda);
    } else if MT == C32 {
        cgerc(m, n, alpha, x, incx, y, incy, a, lda);
    } else {
        zgerc(m, n, alpha, x, incx, y, incy, a, lda);
    }
}

_geqrt2_trmv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int,
    a: $A/ArrayType2D, $$lda: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && X.data_type == MT;
} {
    #if MT == float32 {
        strmv(uplo, trans, diag, n, a, lda, x, incx);
    } else if MT == float64 {
        dtrmv(uplo, trans, diag, n, a, lda, x, incx);
    } else if MT == C32 {
        ctrmv(uplo, trans, diag, n, a, lda, x, incx);
    } else {
        ztrmv(uplo, trans, diag, n, a, lda, x, incx);
    }
}

_geqrt2 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, t: *$T/ArrayType2D, $$ldt: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && T.data_type == MT;
} {
    one_t :: #run if MT == float32 || MT == float64 { one(MT) } else { cpx(one(RT), zero(RT)) };
    zero_t :: #run zero(MT);
    gemv_trans :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };

    routine_name :: #run if MT == float32 {
        "SGEQRT2"
    } else if MT == float64 {
        "DGEQRT2"
    } else if MT == C32 {
        "CGEQRT2"
    } else {
        "ZGEQRT2"
    };

    aii, alpha : MT = ---;
    i, k : int = ---;

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

    k = min(m, n);
    for i:1..k {
        tail := _geqrt2_matrix_column(a.*, lda, min(i + 1, m), i, max(0, m - i));
        _geqrt2_larfg(m - i + 1, _geqrt2_matrix_entry_ptr(a, lda, i, i), *tail, 1, _geqrt2_matrix_entry_ptr(t, ldt, i, 1), MT);
        if i < n {
            aii = fget(a.*, i, i);
            fset(a, i, i, one_t);

            a_right_raw := _geqrt2_matrix_storage(a.*, lda, i, i + 1, m - i + 1, n - i);
            a_right := array2d_view(a_right_raw, lda, row_major=false);
            v := _geqrt2_matrix_column(a.*, lda, i, i, m - i + 1);
            w := _geqrt2_raw_slice(_geqrt2_matrix_entry_ptr(t, ldt, 1, n), max(0, n - i));
            _geqrt2_gemv(gemv_trans, m - i + 1, n - i, one_t, a_right.*, lda, v, 1, zero_t, *w, 1, MT);

            #if MT == C32 || MT == C64 {
                alpha = -conjugate(fget(t.*, i, 1));
            } else {
                alpha = -fget(t.*, i, 1);
            }
            _geqrt2_ger_update(m - i + 1, n - i, alpha, v, 1, w, 1, *a_right, lda, RT, MT);
            fset(a, i, i, aii);
        }
    }

    for i:2..n {
        aii = fget(a.*, i, i);
        fset(a, i, i, one_t);

        a_left_raw := _geqrt2_matrix_storage(a.*, lda, i, 1, m - i + 1, i - 1);
        a_left := array2d_view(a_left_raw, lda, row_major=false);
        v := _geqrt2_matrix_column(a.*, lda, i, i, m - i + 1);
        tcol := _geqrt2_raw_slice(_geqrt2_matrix_entry_ptr(t, ldt, 1, i), i - 1);
        alpha = -fget(t.*, i, 1);
        _geqrt2_gemv(gemv_trans, m - i + 1, i - 1, alpha, a_left.*, lda, v, 1, zero_t, *tcol, 1, MT);
        fset(a, i, i, aii);

        t_tri_raw := _geqrt2_matrix_storage(t.*, ldt, 1, 1, i - 1, i - 1);
        t_tri := array2d_view(t_tri_raw, ldt, row_major=false);
        _geqrt2_trmv(#char "U", #char "N", #char "N", i - 1, t_tri.*, ldt, *tcol, 1, MT);

        fset(t, i, i, fget(t.*, i, 1));
        fset(t, i, 1, zero_t);
    }
}
