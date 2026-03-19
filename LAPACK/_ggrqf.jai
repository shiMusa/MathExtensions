// ? LAPACK *ggrqf()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggrqf/cggrqf/dggrqf/zggrqf
// *GGRQF computes a generalized RQ factorization of A and B.

sggrqf :: #bake_arguments _ggrqf(RT=float32, MT=float32);
dggrqf :: #bake_arguments _ggrqf(RT=float64, MT=float64);
cggrqf :: #bake_arguments _ggrqf(RT=float32, MT=C32);
zggrqf :: #bake_arguments _ggrqf(RT=float64, MT=C64);

#scope_file
_ggrqf_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_ggrqf_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _ggrqf_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_ggrqf_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_ggrqf_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_ggrqf_query_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_ggrqf_gerqf :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sgerqf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == float64 {
        dgerqf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == C32 {
        cgerqf(m, n, a, lda, tau, work, lwork, info);
    } else {
        zgerqf(m, n, a, lda, tau, work, lwork, info);
    }
}

_ggrqf_unmrq :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$k: int,
    a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT && C.data_type == MT && WORK.data_type == MT;
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

_ggrqf_geqrf :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == float64 {
        dgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == C32 {
        cgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else {
        zgeqrf(m, n, a, lda, tau, work, lwork, info);
    }
}

_ggrqf :: (
    $$m: int, $$p: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, taua: *$TAUA/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, taub: *$TAUB/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAUA.data_type == MT &&
           B.data_type == MT && TAUB.data_type == MT &&
           WORK.data_type == MT;
} {
    trans_q :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
    routine_name :: #run if MT == float32 { "SGGRQF" } else if MT == float64 { "DGGRQF" } else if MT == C32 { "CGGRQF" } else { "ZGGRQF" };
    gerqf_name :: #run if MT == float32 { "SGERQF" } else if MT == float64 { "DGERQF" } else if MT == C32 { "CGERQF" } else { "ZGERQF" };
    geqrf_name :: #run if MT == float32 { "SGEQRF" } else if MT == float64 { "DGEQRF" } else if MT == C32 { "CGEQRF" } else { "ZGEQRF" };
    ormrq_name :: #run if MT == float32 { "SORMRQ" } else if MT == float64 { "DORMRQ" } else if MT == C32 { "CUNMRQ" } else { "ZUNMRQ" };

    lquery : bool = ---;
    k, lopt, lwkopt, nb, nb1, nb2, nb3 : int = ---;

    info.* = 0;
    nb1 = ilaenv(1, gerqf_name, " "s, m, n, -1, -1);
    nb2 = ilaenv(1, geqrf_name, " "s, p, n, -1, -1);
    nb3 = ilaenv(1, ormrq_name, " "s, m, n, p, -1);
    nb = max(nb1, nb2, nb3);
    lwkopt = max(1, max(n, m, p)*nb);
    fset(work, 1, _ggrqf_complex_from_real(_ggrqf_roundup_lwork(lwkopt, RT), RT, MT));

    lquery = lwork == -1;
    if m < 0 { info.* = -1; }
    else if p < 0 { info.* = -2; }
    else if n < 0 { info.* = -3; }
    else if lda < max(1, m) { info.* = -5; }
    else if ldb < max(1, p) { info.* = -8; }
    else if lwork < max(1, m, p, n) && !lquery { info.* = -11; }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    _ggrqf_gerqf(m, n, a, lda, taua, work, lwork, info, MT);
    lopt = _ggrqf_query_to_int(fget(work.*, 1), RT, MT);

    k = min(m, n);
    a_rq_row := max(1, m - n + 1);
    a_rq_raw := _ggrqf_matrix_storage(a.*, lda, a_rq_row, 1, k, n);
    a_rq := array2d_view(a_rq_raw, lda, row_major=false);
    _ggrqf_unmrq(#char "R", trans_q, p, n, k, *a_rq, lda, taua, b, ldb, work, lwork, info, RT, MT);
    lopt = max(lopt, _ggrqf_query_to_int(fget(work.*, 1), RT, MT));

    _ggrqf_geqrf(p, n, b, ldb, taub, work, lwork, info, MT);
    fset(work, 1, _ggrqf_complex_from_real(_ggrqf_roundup_lwork(max(lopt, _ggrqf_query_to_int(fget(work.*, 1), RT, MT)), RT), RT, MT));
}
