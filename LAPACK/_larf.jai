// ? LAPACK *larf()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarf/dlarf/clarf/zlarf
// *LARF applies an elementary reflector to a matrix.

slarf :: #bake_arguments _larf(RT=float32, MT=float32);
dlarf :: #bake_arguments _larf(RT=float64, MT=float64);
clarf :: #bake_arguments _larf(RT=float32, MT=C32);
zlarf :: #bake_arguments _larf(RT=float64, MT=C64);

#scope_file
_larf_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_larf_matrix_storage :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _larf_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_larf_nonzero :: inline (x: $MT, $MT: Type) -> bool #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    return x != zero(MT);
}

_larf_last_nonzero_col :: (
    $$m: int, $$n: int, c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           C.data_type == MT;
} -> int {
    if m <= 0 || n <= 0 then return 0;

    for <j:1..n {
        for i:1..m {
            if _larf_nonzero(fget(c.*, i, j), MT) then return j;
        }
    }
    return 0;
}

_larf_last_nonzero_row :: (
    $$m: int, $$n: int, c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           C.data_type == MT;
} -> int {
    if m <= 0 || n <= 0 then return 0;

    for <i:1..m {
        for j:1..n {
            if _larf_nonzero(fget(c.*, i, j), MT) then return i;
        }
    }
    return 0;
}

_larf_gemv :: (
    $$trans: u8, $$m: int, $$n: int, alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    x: $X/ArrayType1D, $$incx: int, beta: MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
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

_larf_ger :: (
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

_larf :: (
    $$side: u8, $$m: int, $$n: int, v: $V/ArrayType1D, $$incv: int, tau: $MT,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    zero_t   :: #run zero(MT);
    one_t    :: #run one(MT);
    minus_tau :: MT = -tau;

    applleft := lsame(side, #char "L");
    trans_op :: u8 = if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
    i, lastv, lastc : int = ---;

    lastv = 0;
    lastc = 0;

    if tau != zero_t {
        if applleft {
            lastv = m;
        } else {
            lastv = n;
        }

        if incv > 0 {
            i = 1 + (lastv - 1)*incv;
        } else {
            i = 1;
        }

        while lastv > 0 && fget(v, i) == zero_t {
            lastv -= 1;
            i -= incv;
        }

        if applleft {
            lastc = _larf_last_nonzero_col(lastv, n, c, ldc, MT);
        } else {
            lastc = _larf_last_nonzero_row(m, lastv, c, ldc, MT);
        }
    }

    if applleft {
        if lastv > 0 {
            c_block_raw := _larf_matrix_storage(c, ldc, 1, 1, lastv, lastc);
            c_block := array2d_view(c_block_raw, ldc, row_major=false);
            _larf_gemv(trans_op, lastv, lastc, one_t, c_block.*, ldc, v, incv, zero_t, work, 1, MT);
            _larf_ger(lastv, lastc, minus_tau, v, incv, work.*, 1, *c_block, ldc, RT, MT);
        }
    } else {
        if lastv > 0 {
            c_block_raw := _larf_matrix_storage(c, ldc, 1, 1, lastc, lastv);
            c_block := array2d_view(c_block_raw, ldc, row_major=false);
            _larf_gemv(#char "N", lastc, lastv, one_t, c_block.*, ldc, v, incv, zero_t, work, 1, MT);
            _larf_ger(lastc, lastv, minus_tau, work.*, 1, v, incv, *c_block, ldc, RT, MT);
        }
    }
}
