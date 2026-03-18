// ? LAPACK *gelq()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgelq/cgelq/dgelq/zgelq
// *GELQ computes an LQ factorization with the modern T-array interface.

sgelq :: #bake_arguments _gelq(RT=float32, MT=float32);
dgelq :: #bake_arguments _gelq(RT=float64, MT=float64);
cgelq :: #bake_arguments _gelq(RT=float32, MT=C32);
zgelq :: #bake_arguments _gelq(RT=float64, MT=C64);

#scope_file
_gelq_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gelq_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_gelq_roundup_lwork :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_gelq_set_query_values :: inline (
    t: *$T/ArrayType1D,
    $$tsize_value: int,
    $$mb: int,
    $$nb: int,
    work: *$WORK/ArrayType1D,
    $$lwork_value: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    fset(t, 1, _gelq_complex_from_real(ncast(RT, tsize_value), RT, MT));
    fset(t, 2, _gelq_complex_from_real(ncast(RT, mb), RT, MT));
    fset(t, 3, _gelq_complex_from_real(ncast(RT, nb), RT, MT));
    fset(work, 1, _gelq_complex_from_real(_gelq_roundup_lwork(ncast(RT, lwork_value), RT), RT, MT));
}

_gelq_gelqt :: (
    $$m: int,
    $$n: int,
    $$mb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType1D,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    t_view := array2d_view(t.*, mb, row_major=false);
    #if MT == float32 {
        sgelqt(m, n, mb, a, lda, *t_view, mb, work, info);
    } else if MT == float64 {
        dgelqt(m, n, mb, a, lda, *t_view, mb, work, info);
    } else if MT == C32 {
        cgelqt(m, n, mb, a, lda, *t_view, mb, work, info);
    } else {
        zgelqt(m, n, mb, a, lda, *t_view, mb, work, info);
    }
}

_gelq_laswlq :: (
    $$m: int,
    $$n: int,
    $$mb: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    t_view := array2d_view(t.*, mb, row_major=false);
    #if MT == float32 {
        slaswlq(m, n, mb, nb, a, lda, *t_view, mb, work, lwork, info);
    } else if MT == float64 {
        dlaswlq(m, n, mb, nb, a, lda, *t_view, mb, work, lwork, info);
    } else if MT == C32 {
        claswlq(m, n, mb, nb, a, lda, *t_view, mb, work, lwork, info);
    } else {
        zlaswlq(m, n, mb, nb, a, lda, *t_view, mb, work, lwork, info);
    }
}

_gelq :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType1D,
    $$tsize: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    routine_name :: #run if MT == float32 {
        "SGELQ"
    } else if MT == float64 {
        "DGELQ"
    } else if MT == C32 {
        "CGELQ"
    } else {
        "ZGELQ"
    };
    env_name :: #run if MT == float32 {
        "SGELQ "
    } else if MT == float64 {
        "DGELQ "
    } else if MT == C32 {
        "CGELQ "
    } else {
        "ZGELQ "
    };

    lquery, lminws, mint, minw : bool = ---;
    mb, nb, mintsz, nblcks, lwmin, lwopt, lwreq : int = ---;

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
        mb = 1;
        nb = n;
    }
    if mb > min(m, n) || mb < 1 { mb = 1; }
    if nb > n || nb <= m { nb = n; }
    mintsz = m + 5;

    if nb > m && n > m {
        if mod(n - m, nb - m) == 0 {
            nblcks = (n - m) / (nb - m);
        } else {
            nblcks = (n - m) / (nb - m) + 1;
        }
    } else {
        nblcks = 1;
    }

    if (n <= m) || (nb <= m) || (nb >= n) {
        lwmin = max(1, n);
        lwopt = max(1, mb*n);
    } else {
        lwmin = max(1, m);
        lwopt = max(1, mb*m);
    }

    lminws = false;
    if (tsize < max(1, mb*m*nblcks + 5) || lwork < lwopt) &&
       (lwork >= lwmin) && (tsize >= mintsz) && !lquery {
        if tsize < max(1, mb*m*nblcks + 5) {
            lminws = true;
            nb = n;
            mb = 1;
        }
        if lwork < lwopt {
            lminws = true;
            mb = 1;
        }
    }

    if (n <= m) || (nb <= m) || (nb >= n) {
        lwreq = max(1, mb*n);
    } else {
        lwreq = max(1, mb*m);
    }

    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    } else if tsize < max(1, mb*m*nblcks + 5) && !lquery && !lminws {
        info.* = -6;
    } else if lwork < lwreq && !lquery && !lminws {
        info.* = -8;
    }

    if info.* == 0 {
        _gelq_set_query_values(t, ifx mint then mintsz else mb*m*nblcks + 5, mb, nb,
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

    t_payload := _gelq_raw_slice(t.data + 5, max(0, t.count - 5));
    t_payload_view : [] MT = t_payload;
    if (n <= m) || (nb <= m) || (nb >= n) {
        _gelq_gelqt(m, n, mb, a, lda, *t_payload_view, work, info, RT, MT);
    } else {
        _gelq_laswlq(m, n, mb, nb, a, lda, *t_payload_view, work, lwork, info, RT, MT);
    }

    fset(work, 1, _gelq_complex_from_real(_gelq_roundup_lwork(ncast(RT, lwreq), RT), RT, MT));
}
