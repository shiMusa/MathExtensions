// ? LAPACK *gehrd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgehrd/cgehrd/dgehrd/zgehrd
// *GEHRD reduces a general matrix to upper Hessenberg form.

sgehrd :: #bake_arguments _gehrd(RT=float32, MT=float32);
dgehrd :: #bake_arguments _gehrd(RT=float64, MT=float64);
cgehrd :: #bake_arguments _gehrd(RT=float32, MT=C32);
zgehrd :: #bake_arguments _gehrd(RT=float64, MT=C64);

#scope_file
NBMAX :: 64;
LDT   :: NBMAX + 1;
TSIZE :: LDT * NBMAX;

_gehrd_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gehrd_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gehrd_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gehrd_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gehrd_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _gehrd_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gehrd_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_gehrd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gehrd_set_work1 :: inline (work: *$WORK/ArrayType1D, $$value: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           WORK.data_type == MT;
} {
    fset(work, 1, _gehrd_complex_from_real(ncast(RT, value), RT, MT));
}

_gehrd_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
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

_gehrd_gemm :: (
    $$transa: u8,
    $$transb: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: MT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
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

_gehrd_trmm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
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

_gehrd_lahr2 :: (
    $$n: int,
    $$k: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
    t: *$T/ArrayType2D,
    $$ldt: int,
    y: *$Y/ArrayType2D,
    $$ldy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           T.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        slahr2(n, k, nb, a, lda, tau, t, ldt, y, ldy);
    } else if MT == float64 {
        dlahr2(n, k, nb, a, lda, tau, t, ldt, y, ldy);
    } else if MT == C32 {
        clahr2(n, k, nb, a, lda, tau, t, ldt, y, ldy);
    } else {
        zlahr2(n, k, nb, a, lda, tau, t, ldt, y, ldy);
    }
}

_gehrd_larfb :: (
    $$side: u8,
    $$trans: u8,
    $$direct: u8,
    $$storev: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    v: $V/ArrayType2D,
    $$ldv: int,
    t: $T/ArrayType2D,
    $$ldt: int,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType2D,
    $$ldwork: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT &&
           T.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        slarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    } else if MT == float64 {
        dlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    } else if MT == C32 {
        clarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    } else {
        zlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
    }
}

_gehrd_gehd2 :: (
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
    #if MT == float32 {
        sgehd2(n, ilo, ihi, a, lda, tau, work, info);
    } else if MT == float64 {
        dgehd2(n, ilo, ihi, a, lda, tau, work, info);
    } else if MT == C32 {
        cgehd2(n, ilo, ihi, a, lda, tau, work, info);
    } else {
        zgehd2(n, ilo, ihi, a, lda, tau, work, info);
    }
}

_gehrd :: (
    $$n: int,
    $$ilo: int,
    $$ihi: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
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
    one_t :: #run one(MT);
    zero_t :: #run zero(MT);
    minus_one_t :: #run -one(MT);
    routine_name :: #run if MT == float32 {
        "SGEHRD"
    } else if MT == float64 {
        "DGEHRD"
    } else if MT == C32 {
        "CGEHRD"
    } else {
        "ZGEHRD"
    };
    trans_t :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };

    lquery : bool = ---;
    i, ib, iinfo, iwt, j, ldwork, lwkopt, nb, nbmin, nh, nx : int = ---;
    ei : MT = ---;

    info.* = 0;
    lquery = lwork == -1;
    if n < 0 {
        info.* = -1;
    } else if ilo < 1 || ilo > max(1, n) {
        info.* = -2;
    } else if ihi < min(ilo, n) || ihi > n {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if lwork < max(1, n) && !lquery {
        info.* = -8;
    }

    nh = ihi - ilo + 1;
    if info.* == 0 {
        if nh <= 1 {
            lwkopt = 1;
        } else {
            nb = min(NBMAX, ilaenv(1, routine_name, " "s, n, ilo, ihi, -1));
            lwkopt = n*nb + TSIZE;
        }
        _gehrd_set_work1(work, lwkopt, RT, MT);
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    for i:1..ilo-1 {
        fset(tau, i, zero_t);
    }
    for i:max(1, ihi)..n-1 {
        fset(tau, i, zero_t);
    }

    if nh <= 1 {
        _gehrd_set_work1(work, 1, RT, MT);
        return;
    }

    nb = min(NBMAX, ilaenv(1, routine_name, " "s, n, ilo, ihi, -1));
    nbmin = 2;
    if nb > 1 && nb < nh {
        nx = max(nb, ilaenv(3, routine_name, " "s, n, ilo, ihi, -1));
        if nx < nh {
            if lwork < lwkopt {
                nbmin = max(2, ilaenv(2, routine_name, " "s, n, ilo, ihi, -1));
                if lwork >= (n*nbmin + TSIZE) {
                    nb = (lwork - TSIZE) / n;
                } else {
                    nb = 1;
                }
            }
        }
    } else {
        nx = nh;
    }
    ldwork = n;

    i = ilo;
    if !(nb < nbmin || nb >= nh) {
        iwt = 1 + n*nb;
        while i <= ihi - 1 - nx {
            ib = min(nb, ihi - i);

            a_panel_raw := _gehrd_matrix_storage(a.*, lda, 1, i, ihi, ib);
            a_panel_view := array2d_view(a_panel_raw, lda, row_major=false);
            tau_view := _gehrd_subvector(tau.*, i, ib);
            t_raw := _gehrd_subvector(work.*, iwt, LDT*ib);
            t_view := array2d_view(t_raw, LDT, row_major=false);
            y_raw := _gehrd_subvector(work.*, 1, ldwork*ib);
            y_view := array2d_view(y_raw, ldwork, row_major=false);
            _gehrd_lahr2(ihi, i, ib, *a_panel_view, lda, tau_view, *t_view, LDT, *y_view, ldwork, MT);

            ei = fget(a.*, i + ib, i + ib - 1);
            fset(a, i + ib, i + ib - 1, one_t);

            b_raw := _gehrd_matrix_storage(a.*, lda, i + ib, i, ihi - i - ib + 1, ib);
            b_view := array2d_view(b_raw, lda, row_major=false);
            c_raw := _gehrd_matrix_storage(a.*, lda, 1, i + ib, ihi, ihi - i - ib + 1);
            c_view := array2d_view(c_raw, lda, row_major=false);
            _gehrd_gemm(#char "N", trans_t, ihi, ihi - i - ib + 1, ib, minus_one_t, y_view, ldwork, b_view, lda, one_t, c_view, lda, MT);

            if ib > 1 && i > 0 {
                v_raw := _gehrd_matrix_storage(a.*, lda, i + 1, i, ib - 1, ib - 1);
                v_view := array2d_view(v_raw, lda, row_major=false);
                y_top_raw := _gehrd_subvector(work.*, 1, ldwork*(ib - 1));
                y_top_view := array2d_view(y_top_raw, ldwork, row_major=false);
                _gehrd_trmm(#char "R", #char "L", trans_t, #char "U", i, ib - 1, one_t, v_view, lda, y_top_view, ldwork, MT);

                for j:0..ib-2 {
                    y_col := _gehrd_matrix_column(y_view.*, ldwork, 1, j + 1, i);
                    a_col := _gehrd_matrix_column(a.*, lda, 1, i + j + 1, i);
                    _gehrd_axpy(i, minus_one_t, y_col, 1, *a_col, 1, MT);
                }
            }

            v_apply_raw := _gehrd_matrix_storage(a.*, lda, i + 1, i, ihi - i, ib);
            v_apply_view := array2d_view(v_apply_raw, lda, row_major=false);
            c_apply_raw := _gehrd_matrix_storage(a.*, lda, i + 1, i + ib, ihi - i, n - i - ib + 1);
            c_apply_view := array2d_view(c_apply_raw, lda, row_major=false);
            work_apply_raw := _gehrd_subvector(work.*, 1, ldwork*ib);
            work_apply_view := array2d_view(work_apply_raw, ldwork, row_major=false);
            _gehrd_larfb(#char "L", trans_t, #char "F", #char "C", ihi - i, n - i - ib + 1, ib, v_apply_view, lda, t_view, LDT, c_apply_view, lda, work_apply_view, ldwork, MT);

            fset(a, i + ib, i + ib - 1, ei);
            i += nb;
        }
    }

    _gehrd_gehd2(n, i, ihi, a, lda, tau, work, *iinfo, RT, MT);
    _gehrd_set_work1(work, lwkopt, RT, MT);
}
