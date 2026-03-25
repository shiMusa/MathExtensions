// ? LAPACK *larf1f()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarf1f/dlarf1f/clarf1f/zlarf1f
// *LARF1F applies an elementary reflector assuming v(1)=1.

slarf1f :: #bake_arguments _larf1f(RT=float32, MT=float32);
dlarf1f :: #bake_arguments _larf1f(RT=float64, MT=float64);
clarf1f :: #bake_arguments _larf1f(RT=float32, MT=C32);
zlarf1f :: #bake_arguments _larf1f(RT=float64, MT=C64);

#scope_file
_larf1f_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_larf1f_matrix_storage :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _larf1f_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_larf1f_matrix_row :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _larf1f_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_larf1f_nonzero :: inline (x: $MT, $MT: Type) -> bool #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    return x != zero(MT);
}

_larf1f_last_nonzero_col :: (
    $$m: int, $$n: int, c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           C.data_type == MT;
} -> int {
    if m <= 0 || n <= 0 then return 0;
    for <j:1..n {
        for i:1..m {
            if _larf1f_nonzero(fget(c.*, i, j), MT) then return j;
        }
    }
    return 0;
}

_larf1f_last_nonzero_row :: (
    $$m: int, $$n: int, c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           C.data_type == MT;
} -> int {
    if m <= 0 || n <= 0 then return 0;
    for <i:1..m {
        for j:1..n {
            if _larf1f_nonzero(fget(c.*, i, j), MT) then return i;
        }
    }
    return 0;
}

_larf1f_gemv :: (
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

_larf1f_ger :: (
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

_larf1f_scal :: (
    $$n: int, alpha: $MT, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sscal(n, alpha, x, incx);
    } else if MT == float64 {
        dscal(n, alpha, x, incx);
    } else if MT == C32 {
        cscal(n, alpha, x, incx);
    } else {
        zscal(n, alpha, x, incx);
    }
}

_larf1f_axpy :: (
    $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        saxpy(n, alpha, x, incx, y, incy);
    } else if MT == float64 {
        daxpy(n, alpha, x, incx, y, incy);
    } else if MT == C32 {
        caxpy(n, alpha, x, incx, y, incy);
    } else {
        zaxpy(n, alpha, x, incx, y, incy);
    }
}

_larf1f :: (
    $$side: u8, $$m: int, $$n: int, v: $V/ArrayType1D, $$incv: int, tau: $MT,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    zero_t :: #run zero(MT);
    one_t  :: #run one(MT);

    applleft := lsame(side, #char "L");
    trans_op :: u8 = if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
    lastv, lastc, i : int = ---;

    lastv = 1;
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

        while lastv > 1 && fget(v, i) == zero_t {
            lastv -= 1;
            i -= incv;
        }

        if applleft {
            lastc = _larf1f_last_nonzero_col(lastv, n, c, ldc, MT);
        } else {
            lastc = _larf1f_last_nonzero_row(m, lastv, c, ldc, MT);
        }
    }

    if lastc == 0 then return;

    if applleft {
        if lastv == 1 {
            row1 := _larf1f_matrix_row(c, ldc, 1, 1, lastc);
            _larf1f_scal(lastc, one_t - tau, *row1, ldc, MT);
        } else {
            c21_raw := _larf1f_matrix_storage(c, ldc, 2, 1, lastv - 1, lastc);
            c21 := array2d_view(c21_raw, ldc, row_major=false);
            v2 := _larf1f_raw_slice(v.data + incv, 1 + (lastv - 2)*abs(incv));
            _larf1f_gemv(trans_op, lastv - 1, lastc, one_t, c21.*, ldc, v2, incv, zero_t, work, 1, MT);

            c1 := _larf1f_matrix_row(c, ldc, 1, 1, lastc);
            if MT == float32 || MT == float64 {
                _larf1f_axpy(lastc, one_t, c1, ldc, work, 1, MT);
                _larf1f_axpy(lastc, -tau, work.*, 1, *c1, ldc, MT);
            } else {
                for j:1..lastc {
                    fset(work, j, fget(work.*, j) + conjugate(fget(c.*, 1, j)));
                    fset(c, 1, j, fget(c.*, 1, j) - tau*conjugate(fget(work.*, j)));
                }
            }

            _larf1f_ger(lastv - 1, lastc, -tau, v2, incv, work.*, 1, *c21, ldc, RT, MT);
        }
    } else {
        if lastv == 1 {
            col1 := _larf1f_raw_slice(c.data, lastc);
            _larf1f_scal(lastc, one_t - tau, *col1, 1, MT);
        } else {
            c12_raw := _larf1f_matrix_storage(c, ldc, 1, 2, lastc, lastv - 1);
            c12 := array2d_view(c12_raw, ldc, row_major=false);
            v2 := _larf1f_raw_slice(v.data + incv, 1 + (lastv - 2)*abs(incv));
            _larf1f_gemv(#char "N", lastc, lastv - 1, one_t, c12.*, ldc, v2, incv, zero_t, work, 1, MT);

            col1 := _larf1f_raw_slice(c.data, lastc);
            _larf1f_axpy(lastc, one_t, col1, 1, work, 1, MT);
            _larf1f_axpy(lastc, -tau, work.*, 1, *col1, 1, MT);

            _larf1f_ger(lastc, lastv - 1, -tau, work.*, 1, v2, incv, *c12, ldc, RT, MT);
        }
    }
}
