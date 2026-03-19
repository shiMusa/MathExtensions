// ? LAPACK *gglse()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgglse/cgglse/dgglse/zgglse
// *GGLSE solves the linear equality-constrained least squares problem.

sgglse :: #bake_arguments _gglse(RT=float32, MT=float32);
dgglse :: #bake_arguments _gglse(RT=float64, MT=float64);
cgglse :: #bake_arguments _gglse(RT=float32, MT=C32);
zgglse :: #bake_arguments _gglse(RT=float64, MT=C64);

#scope_file
_gglse_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gglse_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gglse_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gglse_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gglse_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gglse_query_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_gglse_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_gglse_ggrqf :: (
    $$p: int, $$m: int, $$n: int,
    b: *$B/ArrayType2D, $$ldb: int, taub: *$TAUB/ArrayType1D,
    a: *$A/ArrayType2D, $$lda: int, taua: *$TAUA/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           B.data_type == MT && TAUB.data_type == MT &&
           A.data_type == MT && TAUA.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sggrqf(p, m, n, b, ldb, taub, a, lda, taua, work, lwork, info);
    } else if MT == float64 {
        dggrqf(p, m, n, b, ldb, taub, a, lda, taua, work, lwork, info);
    } else if MT == C32 {
        cggrqf(p, m, n, b, ldb, taub, a, lda, taua, work, lwork, info);
    } else {
        zggrqf(p, m, n, b, ldb, taub, a, lda, taua, work, lwork, info);
    }
}

_gglse_unmqr :: (
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

_gglse_unmrq :: (
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

_gglse_trtrs :: (
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

_gglse_copy :: (
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

_gglse_axpy :: (
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

_gglse_gemv :: (
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

_gglse_trmv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int,
    a: $A/ArrayType2D, $$lda: int, x: *$X/ArrayType1D, $$incx: int,
    $MT: Type
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

_gglse :: (
    $$m: int, $$n: int, $$p: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    c: *$C/ArrayType1D, d: *$D/ArrayType1D, x: *$X/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && B.data_type == MT &&
           C.data_type == MT && D.data_type == MT &&
           X.data_type == MT && WORK.data_type == MT;
} {
    one_t  :: #run _gglse_complex_from_real(one(RT), RT, MT);
    trans_q :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
    routine_name :: #run if MT == float32 { "SGGLSE" } else if MT == float64 { "DGGLSE" } else if MT == C32 { "CGGLSE" } else { "ZGGLSE" };
    geqrf_name :: #run if MT == float32 { "SGEQRF" } else if MT == float64 { "DGEQRF" } else if MT == C32 { "CGEQRF" } else { "ZGEQRF" };
    gerqf_name :: #run if MT == float32 { "SGERQF" } else if MT == float64 { "DGERQF" } else if MT == C32 { "CGERQF" } else { "ZGERQF" };
    ormqr_name :: #run if MT == float32 { "SORMQR" } else if MT == float64 { "DORMQR" } else if MT == C32 { "CUNMQR" } else { "ZUNMQR" };
    ormrq_name :: #run if MT == float32 { "SORMRQ" } else if MT == float64 { "DORMRQ" } else if MT == C32 { "CUNMRQ" } else { "ZUNMRQ" };

    lquery : bool = ---;
    lopt, lwkmin, lwkopt, mn, nb, nb1, nb2, nb3, nb4, nr : int = ---;

    info.* = 0;
    mn = min(m, n);
    lquery = lwork == -1;
    if m < 0 { info.* = -1; }
    else if n < 0 { info.* = -2; }
    else if p < 0 || p > n || p < n - m { info.* = -3; }
    else if lda < max(1, m) { info.* = -5; }
    else if ldb < max(1, p) { info.* = -7; }

    if info.* == 0 {
        if n == 0 {
            lwkmin = 1;
            lwkopt = 1;
        } else {
            nb1 = ilaenv(1, geqrf_name, " "s, m, n, -1, -1);
            nb2 = ilaenv(1, gerqf_name, " "s, m, n, -1, -1);
            nb3 = ilaenv(1, ormqr_name, " "s, m, n, p, -1);
            nb4 = ilaenv(1, ormrq_name, " "s, m, n, p, -1);
            nb = max(nb1, nb2, nb3, nb4);
            lwkmin = m + n + p;
            lwkopt = p + mn + max(m, n)*nb;
        }
        fset(work, 1, _gglse_complex_from_real(_gglse_roundup_lwork(lwkopt, RT), RT, MT));
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

    if n == 0 then return;

    taub := _gglse_subvector(work.*, 1, max(0, p));
    taua := _gglse_subvector(work.*, p + 1, max(0, mn));
    work_tail := _gglse_subvector(work.*, p + mn + 1, max(0, lwork - p - mn));
    _gglse_ggrqf(p, m, n, b, ldb, taub, a, lda, taua, work_tail, lwork - p - mn, info, MT);
    lopt = _gglse_query_to_int(fget(work.*, p + mn + 1), RT, MT);

    c_raw := _gglse_raw_slice(c.data, max(1, m));
    c_mat := array2d_view(c_raw, max(1, m), row_major=false);
    work_tail = _gglse_subvector(work.*, p + mn + 1, max(0, lwork - p - mn));
    _gglse_unmqr(#char "L", trans_q, m, 1, mn, a, lda, taua, c_mat, max(1, m), work_tail, lwork - p - mn, info, RT, MT);
    lopt = max(lopt, _gglse_query_to_int(fget(work.*, p + mn + 1), RT, MT));

    if p > 0 {
        t12_raw := _gglse_matrix_storage(b.*, ldb, 1, n - p + 1, p, p);
        t12 := array2d_view(t12_raw, ldb, row_major=false);
        d_raw := _gglse_raw_slice(d.data, max(1, p));
        d_mat := array2d_view(d_raw, max(1, p), row_major=false);
        _gglse_trtrs(#char "U", #char "N", #char "N", p, 1, t12, ldb, d_mat, p, info, MT);
        if info.* > 0 {
            info.* = 1;
            return;
        }

        x2 := _gglse_subvector(x.*, n - p + 1, p);
        d_vec := _gglse_subvector(d.*, 1, p);
        _gglse_copy(p, d_vec, 1, x2, 1, MT);

        if n - p > 0 {
            a12_raw := _gglse_matrix_storage(a.*, lda, 1, n - p + 1, n - p, p);
            a12 := array2d_view(a12_raw, lda, row_major=false);
            c1 := _gglse_subvector(c.*, 1, n - p);
            _gglse_gemv(#char "N", n - p, p, -one_t, a12, lda, d_vec, 1, one_t, c1, 1, MT);
        }
    }

    if n > p {
        r11_raw := _gglse_matrix_storage(a.*, lda, 1, 1, n - p, n - p);
        r11 := array2d_view(r11_raw, lda, row_major=false);
        c1_raw := _gglse_raw_slice(c.data, max(1, n - p));
        c1_mat := array2d_view(c1_raw, max(1, n - p), row_major=false);
        _gglse_trtrs(#char "U", #char "N", #char "N", n - p, 1, r11, lda, c1_mat, n - p, info, MT);
        if info.* > 0 {
            info.* = 2;
            return;
        }

        c1 := _gglse_subvector(c.*, 1, n - p);
        x1 := _gglse_subvector(x.*, 1, n - p);
        _gglse_copy(n - p, c1, 1, x1, 1, MT);
    }

    if m < n {
        nr = m + p - n;
        if nr > 0 {
            a22_raw := _gglse_matrix_storage(a.*, lda, n - p + 1, m + 1, nr, n - m);
            a22 := array2d_view(a22_raw, lda, row_major=false);
            d_tail := _gglse_subvector(d.*, nr + 1, n - m);
            c_tail := _gglse_subvector(c.*, n - p + 1, nr);
            _gglse_gemv(#char "N", nr, n - m, -one_t, a22, lda, d_tail, 1, one_t, c_tail, 1, MT);
        }
    } else {
        nr = p;
    }

    if nr > 0 {
        a_res_raw := _gglse_matrix_storage(a.*, lda, n - p + 1, n - p + 1, nr, nr);
        a_res := array2d_view(a_res_raw, lda, row_major=false);
        d_head := _gglse_subvector(d.*, 1, nr);
        _gglse_trmv(#char "U", #char "N", #char "N", nr, a_res, lda, d_head, 1, MT);
        c_tail := _gglse_subvector(c.*, n - p + 1, nr);
        _gglse_axpy(nr, -one_t, d_head, 1, c_tail, 1, MT);
    }

    x_raw := _gglse_raw_slice(x.data, max(1, n));
    x_mat := array2d_view(x_raw, max(1, n), row_major=false);
    work_tail = _gglse_subvector(work.*, p + mn + 1, max(0, lwork - p - mn));
    _gglse_unmrq(#char "L", trans_q, n, 1, p, b, ldb, taub, x_mat, n, work_tail, lwork - p - mn, info, RT, MT);
    fset(work, 1, _gglse_complex_from_real(cast(RT) (p + mn + max(lopt, _gglse_query_to_int(fget(work.*, p + mn + 1), RT, MT))), RT, MT));
}
