// ? LAPACK *getsqrhrt()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgetsqrhrt/cgetsqrhrt/dgetsqrhrt/zgetsqrhrt
// *GETSQRHRT computes a blocked QR factorization by TSQR followed by
// Householder reconstruction.

sgetsqrhrt :: #bake_arguments _getsqrhrt(RT=float32, MT=float32);
dgetsqrhrt :: #bake_arguments _getsqrhrt(RT=float64, MT=float64);
cgetsqrhrt :: #bake_arguments _getsqrhrt(RT=float32, MT=C32);
zgetsqrhrt :: #bake_arguments _getsqrhrt(RT=float64, MT=C64);

#scope_file
_getsqrhrt_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_getsqrhrt_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_getsqrhrt_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_getsqrhrt_roundup_lwork :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_getsqrhrt_ceil_div_nonneg :: inline ($$num: int, $$den: int) -> int {
    if num <= 0 {
        return 0;
    }
    return (num + den - 1) / den;
}

_getsqrhrt_latsqr :: (
    $$m: int,
    $$n: int,
    $$mb: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        slatsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    } else if MT == float64 {
        dlatsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    } else if MT == C32 {
        clatsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    } else {
        zlatsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    }
}

_getsqrhrt_orgtsqr_row :: (
    $$m: int,
    $$n: int,
    $$mb: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sorgtsqr_row(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    } else if MT == float64 {
        dorgtsqr_row(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    } else if MT == C32 {
        cungtsqr_row(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    } else {
        zungtsqr_row(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
    }
}

_getsqrhrt_orhr_col :: (
    $$m: int,
    $$n: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
    d: *$D/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           T.data_type == MT &&
           D.data_type == MT;
} {
    #if MT == float32 {
        sorhr_col(m, n, nb, a, lda, t, ldt, d, info);
    } else if MT == float64 {
        dorhr_col(m, n, nb, a, lda, t, ldt, d, info);
    } else if MT == C32 {
        cunhr_col(m, n, nb, a, lda, t, ldt, d, info);
    } else {
        zunhr_col(m, n, nb, a, lda, t, ldt, d, info);
    }
}

_getsqrhrt :: (
    $$m: int,
    $$n: int,
    $$mb1: int,
    $$nb1: int,
    $$nb2: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    t: *$T/ArrayType2D,
    $$ldt: int,
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
        "SGETSQRHRT"
    } else if MT == float64 {
        "DGETSQRHRT"
    } else if MT == C32 {
        "CGETSQRHRT"
    } else {
        "ZGETSQRHRT"
    };

    minus_one_t :: #run -_getsqrhrt_complex_from_real(one(RT), RT, MT);

    iinfo, ldwt, lw1, lw2, lworkopt, lwt, nb1local, nb2local, num_all_row_blocks : int = ---;
    lquery : bool = ---;

    info.* = 0;
    lquery = lwork == -1;

    if m < 0 {
        info.* = -1;
    } else if n < 0 || m < n {
        info.* = -2;
    } else if mb1 <= n {
        info.* = -3;
    } else if nb1 < 1 {
        info.* = -4;
    } else if nb2 < 1 {
        info.* = -5;
    } else if lda < max(1, m) {
        info.* = -7;
    } else if ldt < max(1, min(nb2, n)) {
        info.* = -9;
    } else {
        if lwork < n*n + 1 && !lquery {
            info.* = -11;
        } else {
            nb1local = min(nb1, n);
            num_all_row_blocks = max(1, _getsqrhrt_ceil_div_nonneg(m - n, mb1 - n));
            lwt = num_all_row_blocks*n*nb1local;
            ldwt = nb1local;
            lw1 = nb1local*n;
            lw2 = nb1local*max(nb1local, n - nb1local);
            lworkopt = max(lwt + lw1, max(lwt + n*n + lw2, lwt + n*n + n));
            lworkopt = max(1, lworkopt);
            if lwork < lworkopt && !lquery {
                info.* = -11;
            }
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        fset(work, 1, _getsqrhrt_complex_from_real(_getsqrhrt_roundup_lwork(ncast(RT, lworkopt), RT), RT, MT));
        return;
    }

    if min(m, n) == 0 {
        fset(work, 1, _getsqrhrt_complex_from_real(_getsqrhrt_roundup_lwork(ncast(RT, lworkopt), RT), RT, MT));
        return;
    }

    nb2local = min(nb2, n);

    tsqr_t_data := _getsqrhrt_subvector(work.*, 1, lwt);
    tsqr_t_view := array2d_view(tsqr_t_data, ldwt, row_major=false);
    tsqr_work := _getsqrhrt_subvector(work.*, lwt + 1, lw1);

    iinfo = 0;
    _getsqrhrt_latsqr(m, n, mb1, nb1local, a, lda, *tsqr_t_view, ldwt, *tsqr_work, lw1, *iinfo, MT);

    r_tsqr := _getsqrhrt_subvector(work.*, lwt + 1, n*n);
    for j:1..n {
        for i:1..j {
            fset(r_tsqr, n*(j - 1) + i, fget(a, i, j));
        }
    }

    org_work := _getsqrhrt_subvector(work.*, lwt + n*n + 1, lw2);
    _getsqrhrt_orgtsqr_row(m, n, mb1, nb1local, a, lda, *tsqr_t_view, ldwt, *org_work, lw2, *iinfo, MT);

    diag_work := _getsqrhrt_subvector(work.*, lwt + n*n + 1, n);
    _getsqrhrt_orhr_col(m, n, nb2local, a, lda, t, ldt, *diag_work, *iinfo, MT);

    for i:1..n {
        if fget(diag_work, i) == minus_one_t {
            for j:i..n {
                fset(a, i, j, -fget(r_tsqr, n*(j - 1) + i));
            }
        } else {
            for j:i..n {
                fset(a, i, j, fget(r_tsqr, n*(j - 1) + i));
            }
        }
    }

    fset(work, 1, _getsqrhrt_complex_from_real(_getsqrhrt_roundup_lwork(ncast(RT, lworkopt), RT), RT, MT));
}
