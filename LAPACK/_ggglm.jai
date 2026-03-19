// ? LAPACK *ggglm()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggglm/cggglm/dggglm/zggglm
// *GGGLM solves the Gauss-Markov linear model using a generalized QR
// factorization of (A, B).

sggglm :: #bake_arguments _ggglm(RT=float32, MT=float32);
dggglm :: #bake_arguments _ggglm(RT=float64, MT=float64);
cggglm :: #bake_arguments _ggglm(RT=float32, MT=C32);
zggglm :: #bake_arguments _ggglm(RT=float64, MT=C64);

#scope_file
_ggglm_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_ggglm_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_ggglm_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _ggglm_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_ggglm_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_ggglm_query_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_ggglm_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_ggglm_ggqrf :: (
    $$n: int, $$m: int, $$p: int,
    a: *$A/ArrayType2D, $$lda: int, taua: *$TAUA/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, taub: *$TAUB/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && TAUA.data_type == MT &&
           B.data_type == MT && TAUB.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sggqrf(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
    } else if MT == float64 {
        dggqrf(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
    } else if MT == C32 {
        cggqrf(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
    } else {
        zggqrf(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
    }
}

_ggglm_unmqr :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$k: int,
    a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT &&
           C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == C32 {
        cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_ggglm_unmrq :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$k: int,
    a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT &&
           C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sormrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dormrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == C32 {
        cunmrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_ggglm_trtrs :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int, $$nrhs: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT;
} {
    #if MT == float32 {
        strtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
    } else if MT == float64 {
        dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
    } else if MT == C32 {
        ctrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
    } else {
        ztrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
    }
}

_ggglm_copy :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        scopy(n, x, incx, y, incy);
    } else if MT == float64 {
        dcopy(n, x, incx, y, incy);
    } else if MT == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_ggglm_gemv :: (
    $$trans: u8, $$m: int, $$n: int, alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    x: $X/ArrayType1D, $$incx: int, beta: MT, y: *$Y/ArrayType1D, $$incy: int,
    $MT: Type
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

_ggglm :: (
    $$n: int, $$m: int, $$p: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    d: *$D/ArrayType1D, x: *$X/ArrayType1D, y: *$Y/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && B.data_type == MT &&
           D.data_type == MT && X.data_type == MT && Y.data_type == MT &&
           WORK.data_type == MT;
} {
    zero_t :: #run _ggglm_complex_from_real(zero(RT), RT, MT);
    one_t  :: #run _ggglm_complex_from_real(one(RT), RT, MT);
    trans_q :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
    routine_name :: #run if MT == float32 { "SGGGLM" } else if MT == float64 { "DGGGLM" } else if MT == C32 { "CGGGLM" } else { "ZGGGLM" };
    geqrf_name :: #run if MT == float32 { "SGEQRF" } else if MT == float64 { "DGEQRF" } else if MT == C32 { "CGEQRF" } else { "ZGEQRF" };
    gerqf_name :: #run if MT == float32 { "SGERQF" } else if MT == float64 { "DGERQF" } else if MT == C32 { "CGERQF" } else { "ZGERQF" };
    ormqr_name :: #run if MT == float32 { "SORMQR" } else if MT == float64 { "DORMQR" } else if MT == C32 { "CUNMQR" } else { "ZUNMQR" };
    ormrq_name :: #run if MT == float32 { "SORMRQ" } else if MT == float64 { "DORMRQ" } else if MT == C32 { "CUNMRQ" } else { "ZUNMRQ" };

    lquery : bool = ---;
    i, lopt, lwkmin, lwkopt, nb, nb1, nb2, nb3, nb4, np : int = ---;
    ierr : int = ---;

    info.* = 0;
    np = min(n, p);
    lquery = lwork == -1;
    if n < 0 { info.* = -1; }
    else if m < 0 || m > n { info.* = -2; }
    else if p < 0 || p < n - m { info.* = -3; }
    else if lda < max(1, n) { info.* = -5; }
    else if ldb < max(1, n) { info.* = -7; }

    if info.* == 0 {
        if n == 0 {
            lwkmin = 1;
            lwkopt = 1;
        } else {
            nb1 = ilaenv(1, geqrf_name, " "s, n, m, -1, -1);
            nb2 = ilaenv(1, gerqf_name, " "s, n, m, -1, -1);
            nb3 = ilaenv(1, ormqr_name, " "s, n, m, p, -1);
            nb4 = ilaenv(1, ormrq_name, " "s, n, m, p, -1);
            nb = max(nb1, nb2, nb3, nb4);
            lwkmin = m + n + p;
            lwkopt = m + np + max(n, p)*nb;
        }
        fset(work, 1, _ggglm_complex_from_real(_ggglm_roundup_lwork(lwkopt, RT), RT, MT));
        if lwork < lwkmin && !lquery {
            info.* = -12;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        for i:1..m fset(x, i, zero_t);
        for i:1..p fset(y, i, zero_t);
        return;
    }

    taua := _ggglm_subvector(work.*, 1, max(0, m));
    taub := _ggglm_subvector(work.*, m + 1, max(0, np));
    work_tail := _ggglm_subvector(work.*, m + np + 1, max(0, lwork - m - np));
    _ggglm_ggqrf(n, m, p, a, lda, taua, b, ldb, taub, work_tail, lwork - m - np, info, MT);
    lopt = _ggglm_query_to_int(fget(work.*, m + np + 1), RT, MT);

    d_raw := _ggglm_raw_slice(d.data, max(1, n));
    d_mat := array2d_view(d_raw, max(1, n), row_major=false);
    work_tail = _ggglm_subvector(work.*, m + np + 1, max(0, lwork - m - np));
    _ggglm_unmqr(#char "L", trans_q, n, 1, m, a, lda, taua, d_mat, max(1, n), work_tail, lwork - m - np, info, RT, MT);
    lopt = max(lopt, _ggglm_query_to_int(fget(work.*, m + np + 1), RT, MT));

    if n > m {
        b22_raw := _ggglm_matrix_storage(b.*, ldb, m + 1, m + p - n + 1, n - m, n - m);
        b22 := array2d_view(b22_raw, ldb, row_major=false);
        d2_raw := _ggglm_raw_slice(d.data + m, max(1, n - m));
        d2_mat := array2d_view(d2_raw, max(1, n - m), row_major=false);
        _ggglm_trtrs(#char "U", #char "N", #char "N", n - m, 1, b22, ldb, d2_mat, n - m, info, MT);
        if info.* > 0 {
            info.* = 1;
            return;
        }
        d2 := _ggglm_subvector(d.*, m + 1, n - m);
        y2 := _ggglm_subvector(y.*, m + p - n + 1, n - m);
        _ggglm_copy(n - m, d2, 1, y2, 1, MT);
    }

    for i:1..m + p - n fset(y, i, zero_t);

    if n > m {
        t12_raw := _ggglm_matrix_storage(b.*, ldb, 1, m + p - n + 1, m, n - m);
        t12 := array2d_view(t12_raw, ldb, row_major=false);
        y2 := _ggglm_subvector(y.*, m + p - n + 1, n - m);
        d1 := _ggglm_subvector(d.*, 1, m);
        _ggglm_gemv(#char "N", m, n - m, -one_t, t12, ldb, y2, 1, one_t, d1, 1, MT);
    }

    if m > 0 {
        r11_raw := _ggglm_matrix_storage(a.*, lda, 1, 1, m, m);
        r11 := array2d_view(r11_raw, lda, row_major=false);
        d1_raw := _ggglm_raw_slice(d.data, max(1, m));
        d1_mat := array2d_view(d1_raw, max(1, m), row_major=false);
        _ggglm_trtrs(#char "U", #char "N", #char "N", m, 1, r11, lda, d1_mat, m, info, MT);
        if info.* > 0 {
            info.* = 2;
            return;
        }
        d1 := _ggglm_subvector(d.*, 1, m);
        _ggglm_copy(m, d1, 1, x, 1, MT);
    }

    b_rq_row := max(1, n - p + 1);
    b_rq_raw := _ggglm_matrix_storage(b.*, ldb, b_rq_row, 1, np, p);
    b_rq := array2d_view(b_rq_raw, ldb, row_major=false);
    y_raw := _ggglm_raw_slice(y.data, max(1, p));
    y_mat := array2d_view(y_raw, max(1, p), row_major=false);
    work_tail = _ggglm_subvector(work.*, m + np + 1, max(0, lwork - m - np));
    _ggglm_unmrq(#char "L", trans_q, p, 1, np, b_rq, ldb, taub, y_mat, max(1, p), work_tail, lwork - m - np, info, RT, MT);
    fset(work, 1, _ggglm_complex_from_real(cast(RT) (m + np + max(lopt, _ggglm_query_to_int(fget(work.*, m + np + 1), RT, MT))), RT, MT));
}
