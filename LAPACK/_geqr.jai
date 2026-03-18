// ? LAPACK *geqr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeqr/cgeqr/dgeqr/zgeqr
// *GEQR computes a QR factorization with the modern T-array interface.

sgeqr :: #bake_arguments _geqr(RT=float32, MT=float32);
dgeqr :: #bake_arguments _geqr(RT=float64, MT=float64);
cgeqr :: #bake_arguments _geqr(RT=float32, MT=C32);
zgeqr :: #bake_arguments _geqr(RT=float64, MT=C64);

#scope_file
_geqr_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geqr_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_geqr_roundup_lwork :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_geqr_set_query_values :: inline (
    t: *$T/ArrayType1D, $$tsize_value: int, $$mb: int, $$nb: int,
    work: *$WORK/ArrayType1D, $$lwork_value: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           T.data_type == MT && WORK.data_type == MT;
} {
    fset(t, 1, _geqr_complex_from_real(ncast(RT, tsize_value), RT, MT));
    fset(t, 2, _geqr_complex_from_real(ncast(RT, mb), RT, MT));
    fset(t, 3, _geqr_complex_from_real(ncast(RT, nb), RT, MT));
    fset(work, 1, _geqr_complex_from_real(_geqr_roundup_lwork(ncast(RT, lwork_value), RT), RT, MT));
}

_geqr_geqrt :: (
    $$m: int, $$n: int, $$nb: int, a: *$A/ArrayType2D, $$lda: int,
    t: *$T/ArrayType1D, work: *$WORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && T.data_type == MT && WORK.data_type == MT;
} {
    t_view := array2d_view(t.*, nb, row_major=false);
    #if MT == float32 {
        sgeqrt(m, n, nb, a, lda, *t_view, nb, work, info);
    } else if MT == float64 {
        dgeqrt(m, n, nb, a, lda, *t_view, nb, work, info);
    } else if MT == C32 {
        cgeqrt(m, n, nb, a, lda, *t_view, nb, work, info);
    } else {
        zgeqrt(m, n, nb, a, lda, *t_view, nb, work, info);
    }
}

_geqr_latsqr :: (
    $$m: int, $$n: int, $$mb: int, $$nb: int, a: *$A/ArrayType2D, $$lda: int,
    t: *$T/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && T.data_type == MT && WORK.data_type == MT;
} {
    t_view := array2d_view(t.*, nb, row_major=false);
    #if MT == float32 {
        slatsqr(m, n, mb, nb, a, lda, *t_view, nb, work, lwork, info);
    } else if MT == float64 {
        dlatsqr(m, n, mb, nb, a, lda, *t_view, nb, work, lwork, info);
    } else if MT == C32 {
        clatsqr(m, n, mb, nb, a, lda, *t_view, nb, work, lwork, info);
    } else {
        zlatsqr(m, n, mb, nb, a, lda, *t_view, nb, work, lwork, info);
    }
}

_geqr :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    t: *$T/ArrayType1D, $$tsize: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT && T.data_type == MT && WORK.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGEQR"
    } else if MT == float64 {
        "DGEQR"
    } else if MT == C32 {
        "CGEQR"
    } else {
        "ZGEQR"
    };
    env_name :: #run if MT == float32 {
        "SGEQR "
    } else if MT == float64 {
        "DGEQR "
    } else if MT == C32 {
        "CGEQR "
    } else {
        "ZGEQR "
    };

    lquery, lminws, mint, minw : bool = ---;
    mb, nb, mintsz, nblcks, lwmin, lwreq : int = ---;

    info.* = 0;
    lquery = (tsize == -1 || tsize == -2 || lwork == -1 || lwork == -2);

    mint = false;
    minw = false;
    if tsize == -2 || lwork == -2 {
        if tsize != -1 { mint = true; }
        if lwork != -1 { minw = true; }
    }

    if min(m, n) > 0 {
        mb = ilaenv(1, env_name, " "s, m, n, 1, -1);
        nb = ilaenv(1, env_name, " "s, m, n, 2, -1);
    } else {
        mb = m;
        nb = 1;
    }
    if mb > m || mb <= n { mb = m; }
    if nb > min(m, n) || nb < 1 { nb = 1; }
    mintsz = n + 5;

    if mb > n && m > n {
        if mod(m - n, mb - n) == 0 {
            nblcks = (m - n) / (mb - n);
        } else {
            nblcks = (m - n) / (mb - n) + 1;
        }
    } else {
        nblcks = 1;
    }

    lwmin = max(1, n);
    lwreq = max(1, n*nb);
    lminws = false;
    if (tsize < max(1, nb*n*nblcks + 5) || lwork < lwreq) &&
       (lwork >= n) && (tsize >= mintsz) && !lquery {
        if tsize < max(1, nb*n*nblcks + 5) {
            lminws = true;
            nb = 1;
            mb = m;
        }
        if lwork < lwreq {
            lminws = true;
            nb = 1;
        }
    }

    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    } else if tsize < max(1, nb*n*nblcks + 5) && !lquery && !lminws {
        info.* = -6;
    } else if lwork < lwreq && !lquery && !lminws {
        info.* = -8;
    }

    if info.* == 0 {
        _geqr_set_query_values(t, ifx mint then mintsz else nb*n*nblcks + 5, mb, nb,
                               work, ifx minw then lwmin else lwreq, RT, MT);
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if min(m, n) == 0 {
        return;
    }

    t_payload := _geqr_raw_slice(t.data + 5, max(0, t.count - 5));
    t_payload_view : [] MT = t_payload;
    if (m <= n) || (mb <= n) || (mb >= m) {
        _geqr_geqrt(m, n, nb, a, lda, *t_payload_view, work, info, RT, MT);
    } else {
        _geqr_latsqr(m, n, mb, nb, a, lda, *t_payload_view, work, lwork, info, RT, MT);
    }

    fset(work, 1, _geqr_complex_from_real(_geqr_roundup_lwork(ncast(RT, lwreq), RT), RT, MT));
}
