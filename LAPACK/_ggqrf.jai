// ? LAPACK *ggqrf()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggqrf/cggqrf/dggqrf/zggqrf
// *GGQRF computes a generalized QR factorization of A and B.

sggqrf :: #bake_arguments _ggqrf(RT=float32, MT=float32);
dggqrf :: #bake_arguments _ggqrf(RT=float64, MT=float64);
cggqrf :: #bake_arguments _ggqrf(RT=float32, MT=C32);
zggqrf :: #bake_arguments _ggqrf(RT=float64, MT=C64);

#scope_file
_ggqrf_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_ggqrf_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_ggqrf_query_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_ggqrf_geqrf :: (
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

_ggqrf_unmqr :: (
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
        sormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == C32 {
        cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_ggqrf_gerqf :: (
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

_ggqrf :: (
    $$n: int, $$m: int, $$p: int,
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
    routine_name :: #run if MT == float32 { "SGGQRF" } else if MT == float64 { "DGGQRF" } else if MT == C32 { "CGGQRF" } else { "ZGGQRF" };
    geqrf_name :: #run if MT == float32 { "SGEQRF" } else if MT == float64 { "DGEQRF" } else if MT == C32 { "CGEQRF" } else { "ZGEQRF" };
    gerqf_name :: #run if MT == float32 { "SGERQF" } else if MT == float64 { "DGERQF" } else if MT == C32 { "CGERQF" } else { "ZGERQF" };
    ormqr_name :: #run if MT == float32 { "SORMQR" } else if MT == float64 { "DORMQR" } else if MT == C32 { "CUNMQR" } else { "ZUNMQR" };

    lquery : bool = ---;
    lopt, lwkopt, nb, nb1, nb2, nb3 : int = ---;

    info.* = 0;
    nb1 = ilaenv(1, geqrf_name, " "s, n, m, -1, -1);
    nb2 = ilaenv(1, gerqf_name, " "s, n, p, -1, -1);
    nb3 = ilaenv(1, ormqr_name, " "s, n, m, p, -1);
    nb = max(nb1, nb2, nb3);
    lwkopt = max(1, max(n, m, p)*nb);
    fset(work, 1, _ggqrf_complex_from_real(_ggqrf_roundup_lwork(lwkopt, RT), RT, MT));

    lquery = lwork == -1;
    if n < 0 { info.* = -1; }
    else if m < 0 { info.* = -2; }
    else if p < 0 { info.* = -3; }
    else if lda < max(1, n) { info.* = -5; }
    else if ldb < max(1, n) { info.* = -8; }
    else if lwork < max(1, n, m, p) && !lquery { info.* = -11; }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    _ggqrf_geqrf(n, m, a, lda, taua, work, lwork, info, MT);
    lopt = _ggqrf_query_to_int(fget(work.*, 1), RT, MT);

    _ggqrf_unmqr(#char "L", trans_q, n, p, min(n, m), a, lda, taua, b, ldb, work, lwork, info, RT, MT);
    lopt = max(lopt, _ggqrf_query_to_int(fget(work.*, 1), RT, MT));

    _ggqrf_gerqf(n, p, b, ldb, taub, work, lwork, info, MT);
    fset(work, 1, _ggqrf_complex_from_real(_ggqrf_roundup_lwork(max(lopt, _ggqrf_query_to_int(fget(work.*, 1), RT, MT)), RT), RT, MT));
}
