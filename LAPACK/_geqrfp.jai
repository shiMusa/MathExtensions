// ? LAPACK *geqrfp()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeqrfp/cgeqrfp/dgeqrfp/zgeqrfp
// *GEQRFP computes a QR factorization with nonnegative diagonal entries.

sgeqrfp :: #bake_arguments _geqrfp(RT=float32, MT=float32);
dgeqrfp :: #bake_arguments _geqrfp(RT=float64, MT=float64);
cgeqrfp :: #bake_arguments _geqrfp(RT=float32, MT=C32);
zgeqrfp :: #bake_arguments _geqrfp(RT=float64, MT=C64);

#scope_file
_geqrfp_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geqrfp_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_geqrfp_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _geqrfp_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geqrfp_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_geqrfp_roundup_lwork :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_geqrfp_set_work1 :: inline (work: *$WORK/ArrayType1D, $$value: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           WORK.data_type == MT;
} {
    #if MT == float32 || MT == float64 {
        fset(work, 1, ncast(MT, value));
    } else {
        fset(work, 1, _geqrfp_complex_from_real(ncast(RT, value), RT, MT));
    }
}

_geqrfp_larft :: (
    $$direct: u8, $$storev: u8, $$n: int, $$k: int,
    v: $V/ArrayType2D, $$ldv: int, tau: $TAU/ArrayType1D, t: *$T/ArrayType2D, $$ldt: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT && TAU.data_type == MT && T.data_type == MT;
} {
    #if MT == float32 {
        slarft(direct, storev, n, k, v, ldv, tau, t, ldt);
    } else if MT == float64 {
        dlarft(direct, storev, n, k, v, ldv, tau, t, ldt);
    } else if MT == C32 {
        clarft(direct, storev, n, k, v, ldv, tau, t, ldt);
    } else {
        zlarft(direct, storev, n, k, v, ldv, tau, t, ldt);
    }
}

_geqrfp_larfb :: (
    $$side: u8, $$trans: u8, $$direct: u8, $$storev: u8,
    $$m: int, $$n: int, $$k: int,
    v: $V/ArrayType2D, $$ldv: int, t: $T/ArrayType2D, $$ldt: int,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType2D, $$ldwork: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT && T.data_type == MT && C.data_type == MT && WORK.data_type == MT;
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

_geqrfp_geqr2p :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == float32 {
        sgeqr2p(m, n, a, lda, tau, work, info);
    } else if MT == float64 {
        dgeqr2p(m, n, a, lda, tau, work, info);
    } else if MT == C32 {
        cgeqr2p(m, n, a, lda, tau, work, info);
    } else {
        zgeqr2p(m, n, a, lda, tau, work, info);
    }
}

_geqrfp :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGEQRFP"
    } else if MT == float64 {
        "DGEQRFP"
    } else if MT == C32 {
        "CGEQRFP"
    } else {
        "ZGEQRFP"
    };
    env_name :: #run if MT == float32 {
        "SGEQRF"
    } else if MT == float64 {
        "DGEQRF"
    } else if MT == C32 {
        "CGEQRF"
    } else {
        "ZGEQRF"
    };

    lquery : bool = ---;
    i, ib, iinfo, iws, k, ldwork, lwkmin, lwkopt, nb, nbmin, nx : int = ---;

    info.* = 0;
    nb = ilaenv(1, env_name, " "s, m, n, -1, -1);
    k = min(m, n);
    if k == 0 {
        lwkmin = 1;
        lwkopt = 1;
    } else {
        lwkmin = n;
        lwkopt = n*nb;
    }

    #if MT == float32 || MT == float64 {
        fset(work, 1, _geqrfp_roundup_lwork(ncast(RT, lwkopt), RT));
    } else {
        fset(work, 1, _geqrfp_complex_from_real(_geqrfp_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
    }

    lquery = lwork == -1;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    } else if lwork < lwkmin && !lquery {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if k == 0 {
        _geqrfp_set_work1(work, 1, RT, MT);
        return;
    }

    nbmin = 2;
    nx = 0;
    iws = lwkmin;
    if nb > 1 && nb < k {
        nx = max(0, ilaenv(3, env_name, " "s, m, n, -1, -1));
        if nx < k {
            ldwork = n;
            iws = ldwork*nb;
            if lwork < iws {
                nb = lwork / ldwork;
                nbmin = max(2, ilaenv(2, env_name, " "s, m, n, -1, -1));
            }
        }
    }

    if nb >= nbmin && nb < k && nx < k {
        i = 1;
        while i <= k - nx {
            ib = min(k - i + 1, nb);

            a_block_raw := _geqrfp_matrix_storage(a.*, lda, i, i, m - i + 1, ib);
            a_block := array2d_view(a_block_raw, lda, row_major=false);
            tau_block := _geqrfp_subvector(tau.*, i, ib);
            _geqrfp_geqr2p(m - i + 1, ib, *a_block, lda, *tau_block, work, *iinfo, RT, MT);

            if i + ib <= n {
                t_block_raw := _geqrfp_raw_slice(work.data, ldwork*ib);
                t_block := array2d_view(t_block_raw, ldwork, row_major=false);
                _geqrfp_larft(#char "F", #char "C", m - i + 1, ib, a_block.*, lda, tau_block, *t_block, ldwork, MT);

                c_raw := _geqrfp_matrix_storage(a.*, lda, i, i + ib, m - i + 1, n - i - ib + 1);
                c_view := array2d_view(c_raw, lda, row_major=false);
                work2_raw := _geqrfp_raw_slice(work.data + ib, ldwork*(n - i - ib + 1));
                work2_view := array2d_view(work2_raw, ldwork, row_major=false);
                #if MT == float32 || MT == float64 {
                    _geqrfp_larfb(#char "L", #char "T", #char "F", #char "C", m - i + 1, n - i - ib + 1, ib,
                                  a_block.*, lda, t_block.*, ldwork, *c_view, lda, *work2_view, ldwork, MT);
                } else {
                    _geqrfp_larfb(#char "L", #char "C", #char "F", #char "C", m - i + 1, n - i - ib + 1, ib,
                                  a_block.*, lda, t_block.*, ldwork, *c_view, lda, *work2_view, ldwork, MT);
                }
            }

            i += nb;
        }
    } else {
        i = 1;
    }

    if i <= k {
        a_rem_raw := _geqrfp_matrix_storage(a.*, lda, i, i, m - i + 1, n - i + 1);
        a_rem := array2d_view(a_rem_raw, lda, row_major=false);
        tau_rem := _geqrfp_subvector(tau.*, i, k - i + 1);
        _geqrfp_geqr2p(m - i + 1, n - i + 1, *a_rem, lda, *tau_rem, work, *iinfo, RT, MT);
    }

    #if MT == float32 || MT == float64 {
        fset(work, 1, _geqrfp_roundup_lwork(ncast(RT, iws), RT));
    } else {
        fset(work, 1, _geqrfp_complex_from_real(_geqrfp_roundup_lwork(ncast(RT, iws), RT), RT, MT));
    }
}
