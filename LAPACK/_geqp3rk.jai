// ? LAPACK *geqp3rk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeqp3rk/cgeqp3rk/dgeqp3rk/zgeqp3rk
// *GEQP3RK computes a truncated QR factorization with column pivoting.

sgeqp3rk :: (
    $$m: int, $$n: int, $$nrhs: int, $$kmax: int, abstol: float32, reltol: float32,
    a: *$A/ArrayType2D, $$lda: int, k: *int, maxc2nrmk: *float32, relmaxc2nrmk: *float32,
    jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, info: *int
) {
    _geqp3rk_real(m, n, nrhs, kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk, jpiv, tau, work, lwork, iwork, info, RT=float32);
}

dgeqp3rk :: (
    $$m: int, $$n: int, $$nrhs: int, $$kmax: int, abstol: float64, reltol: float64,
    a: *$A/ArrayType2D, $$lda: int, k: *int, maxc2nrmk: *float64, relmaxc2nrmk: *float64,
    jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, info: *int
) {
    _geqp3rk_real(m, n, nrhs, kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk, jpiv, tau, work, lwork, iwork, info, RT=float64);
}

cgeqp3rk :: (
    $$m: int, $$n: int, $$nrhs: int, $$kmax: int, abstol: float32, reltol: float32,
    a: *$A/ArrayType2D, $$lda: int, k: *int, maxc2nrmk: *float32, relmaxc2nrmk: *float32,
    jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, info: *int
) {
    _geqp3rk_complex(m, n, nrhs, kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk, jpiv, tau, work, lwork, rwork, iwork, info, RT=float32, MT=C32);
}

zgeqp3rk :: (
    $$m: int, $$n: int, $$nrhs: int, $$kmax: int, abstol: float64, reltol: float64,
    a: *$A/ArrayType2D, $$lda: int, k: *int, maxc2nrmk: *float64, relmaxc2nrmk: *float64,
    jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, info: *int
) {
    _geqp3rk_complex(m, n, nrhs, kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk, jpiv, tau, work, lwork, rwork, iwork, info, RT=float64, MT=C64);
}

#scope_file
_geqp3rk_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geqp3rk_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_geqp3rk_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _geqp3rk_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geqp3rk_roundup_lwork :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_geqp3rk_zero_mt :: inline ($RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return zero(MT);
    } else {
        return cpx(zero(RT), zero(RT));
    }
}

_geqp3rk_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && (MT == C32 || MT == C64);
} {
    return cpx(x, zero(RT));
}

_geqp3rk_set_work1 :: inline (work: *$WORK/ArrayType1D, value: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           WORK.data_type == MT;
} {
    #if MT == float32 || MT == float64 {
        fset(work, 1, ncast(MT, value));
    } else {
        fset(work, 1, _geqp3rk_complex_from_real(value, RT, MT));
    }
}

_geqp3rk_lamch :: ($$cmach: u8, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(cmach);
    } else {
        return dlamch(cmach);
    }
}

_geqp3rk_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_geqp3rk_amax_index :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return isamax(n, x, incx);
    } else {
        return idamax(n, x, incx);
    }
}

_geqp3rk_nrm2_real :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return snrm2(n, x, incx);
    } else {
        return dnrm2(n, x, incx);
    }
}

_geqp3rk_nrm2_complex :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        return scnrm2(n, x, incx);
    } else {
        return dznrm2(n, x, incx);
    }
}

_geqp3rk_laqp3rk_real :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$nb: int, $$abstol: $RT, $$reltol: $RT,
    $$kp1: int, $$maxc2nrm: $RT, a: *$A/ArrayType2D, $$lda: int, done: *bool, kb: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, auxv: *$AUXV/ArrayType1D,
    f: *$F/ArrayType2D, $$ldf: int, iwork: *$IWORK/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPIV.data_type == int &&
           TAU.data_type == RT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           AUXV.data_type == RT &&
           F.data_type == RT &&
           IWORK.data_type == int;
} {
    #if RT == float32 {
        slaqp3rk(m, n, nrhs, ioffset, nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb, maxc2nrmk, relmaxc2nrmk,
                 jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info);
    } else {
        dlaqp3rk(m, n, nrhs, ioffset, nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb, maxc2nrmk, relmaxc2nrmk,
                 jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info);
    }
}

_geqp3rk_laqp2rk_real :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$kmax_: int, $$abstol: $RT, $$reltol: $RT,
    $$kp1: int, $$maxc2nrm: $RT, a: *$A/ArrayType2D, $$lda: int, kf: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, work: *$WORK/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPIV.data_type == int &&
           TAU.data_type == RT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           WORK.data_type == RT;
} {
    #if RT == float32 {
        slaqp2rk(m, n, nrhs, ioffset, kmax_, abstol, reltol, kp1, maxc2nrm, a, lda, kf, maxc2nrmk, relmaxc2nrmk,
                 jpiv, tau, vn1, vn2, work, info);
    } else {
        dlaqp2rk(m, n, nrhs, ioffset, kmax_, abstol, reltol, kp1, maxc2nrm, a, lda, kf, maxc2nrmk, relmaxc2nrmk,
                 jpiv, tau, vn1, vn2, work, info);
    }
}

_geqp3rk_laqp3rk_complex :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$nb: int, $$abstol: $RT, $$reltol: $RT,
    $$kp1: int, $$maxc2nrm: $RT, a: *$A/ArrayType2D, $$lda: int, done: *bool, kb: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, auxv: *$AUXV/ArrayType1D,
    f: *$F/ArrayType2D, $$ldf: int, iwork: *$IWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           JPIV.data_type == int &&
           TAU.data_type == MT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           AUXV.data_type == MT &&
           F.data_type == MT &&
           IWORK.data_type == int;
} {
    #if MT == C32 {
        claqp3rk(m, n, nrhs, ioffset, nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb, maxc2nrmk, relmaxc2nrmk,
                 jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info);
    } else {
        zlaqp3rk(m, n, nrhs, ioffset, nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb, maxc2nrmk, relmaxc2nrmk,
                 jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info);
    }
}

_geqp3rk_laqp2rk_complex :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$kmax_: int, $$abstol: $RT, $$reltol: $RT,
    $$kp1: int, $$maxc2nrm: $RT, a: *$A/ArrayType2D, $$lda: int, kf: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, work: *$WORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           JPIV.data_type == int &&
           TAU.data_type == MT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           WORK.data_type == MT;
} {
    #if MT == C32 {
        claqp2rk(m, n, nrhs, ioffset, kmax_, abstol, reltol, kp1, maxc2nrm, a, lda, kf, maxc2nrmk, relmaxc2nrmk,
                 jpiv, tau, vn1, vn2, work, info);
    } else {
        zlaqp2rk(m, n, nrhs, ioffset, kmax_, abstol, reltol, kp1, maxc2nrm, a, lda, kf, maxc2nrmk, relmaxc2nrmk,
                 jpiv, tau, vn1, vn2, work, info);
    }
}

_geqp3rk_real :: (
    $$m: int, $$n: int, $$nrhs: int, $$kmax: int, abstol: $RT, reltol: $RT,
    a: *$A/ArrayType2D, $$lda: int, k: *int, maxc2nrmk: *RT, relmaxc2nrmk: *RT,
    jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D,
    $$lwork: int, iwork: *$IWORK/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPIV.data_type == int &&
           TAU.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 then "SGEQP3RK" else "DGEQP3RK";
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run one_r + one_r;

    lquery, done : bool = ---;
    iinfo, ioffset, iws, j, jb, jbf, jmaxb, jmax, jmaxc2nrm, kp1, kf, lwkopt, minmn, n_sub, nb, nbmin, nx : int = ---;
    eps, hugeval, maxc2nrm, safmin : RT = ---;

    info.* = 0;
    lquery = lwork == -1;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if kmax < 0 {
        info.* = -4;
    } else if _geqp3rk_isnan(abstol, RT) {
        info.* = -5;
    } else if _geqp3rk_isnan(reltol, RT) {
        info.* = -6;
    } else if lda < max(1, m) {
        info.* = -8;
    }

    if info.* == 0 {
        minmn = min(m, n);
        if minmn == 0 {
            iws = 1;
            lwkopt = 1;
        } else {
            iws = 3*n + nrhs - 1;
            nb = ilaenv(1, routine_name, " "s, m, n, -1, -1);
            lwkopt = 2*n + nb*(n + nrhs + 1);
        }
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, RT);
        if lwork < iws && !lquery {
            info.* = -15;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if minmn == 0 {
        k.* = 0;
        maxc2nrmk.* = zero_r;
        relmaxc2nrmk.* = zero_r;
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, RT);
        return;
    }

    for j:1..n {
        fset(jpiv, j, j);
    }

    for j:1..n {
        col := _geqp3rk_raw_slice(a.data + (j - 1)*lda, m);
        norm := _geqp3rk_nrm2_real(m, col, 1, RT);
        fset(work, j, norm);
        fset(work, n + j, norm);
    }

    norms := _geqp3rk_subvector(work.*, 1, n);
    kp1 = _geqp3rk_amax_index(n, norms, 1, RT);
    maxc2nrm = fget(work.*, kp1);

    if _geqp3rk_isnan(maxc2nrm, RT) {
        k.* = 0;
        info.* = kp1;
        maxc2nrmk.* = maxc2nrm;
        relmaxc2nrmk.* = maxc2nrm;
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, RT);
        return;
    }

    if maxc2nrm == zero_r {
        k.* = 0;
        maxc2nrmk.* = zero_r;
        relmaxc2nrmk.* = zero_r;
        for j:1..minmn {
            fset(tau, j, zero_r);
        }
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, RT);
        return;
    }

    hugeval = _geqp3rk_lamch(#char "o", RT);
    if maxc2nrm > hugeval {
        info.* = n + kp1;
    }

    if kmax == 0 {
        k.* = 0;
        maxc2nrmk.* = maxc2nrm;
        relmaxc2nrmk.* = one_r;
        for j:1..minmn {
            fset(tau, j, zero_r);
        }
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, RT);
        return;
    }

    eps = _geqp3rk_lamch(#char "e", RT);
    if abstol >= zero_r {
        safmin = _geqp3rk_lamch(#char "s", RT);
        abstol = max(abstol, two_r*safmin);
    }
    if reltol >= zero_r {
        reltol = max(reltol, eps);
    }

    jmax = min(kmax, minmn);
    if maxc2nrm <= abstol || one_r <= reltol {
        k.* = 0;
        maxc2nrmk.* = maxc2nrm;
        relmaxc2nrmk.* = one_r;
        for j:1..minmn {
            fset(tau, j, zero_r);
        }
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, RT);
        return;
    }

    nbmin = 2;
    nx = 0;
    if nb > 1 && nb < minmn {
        nx = max(0, ilaenv(3, routine_name, " "s, m, n, -1, -1));
        if nx < minmn {
            if lwork < lwkopt {
                nb = (lwork - 2*n) / (n + 1);
                nbmin = max(2, ilaenv(2, routine_name, " "s, m, n, -1, -1));
            }
        }
    }

    done = false;
    j = 1;
    jmaxb = min(kmax, minmn - nx);
    if nb >= nbmin && nb < jmax && jmaxb > 0 {
        while j <= jmaxb {
            jb = min(nb, jmaxb - j + 1);
            n_sub = n - j + 1;
            ioffset = j - 1;

            a_sub_raw := _geqp3rk_matrix_storage(a.*, lda, 1, j, m, n_sub + nrhs);
            a_sub := array2d_view(a_sub_raw, lda, row_major=false);
            jpiv_sub := _geqp3rk_subvector(jpiv.*, j, n_sub);
            tau_sub := _geqp3rk_subvector(tau.*, j, n_sub);
            vn1 := _geqp3rk_subvector(work.*, j, n_sub);
            vn2 := _geqp3rk_subvector(work.*, n + j, n_sub);
            auxv := _geqp3rk_subvector(work.*, 2*n + 1, jb);
            f_raw := _geqp3rk_raw_slice(work.data + (2*n + jb), max(0, (n_sub + nrhs)*jb));
            f_view := array2d_view(f_raw, n_sub + nrhs, row_major=false);

            _geqp3rk_laqp3rk_real(m, n_sub, nrhs, ioffset, jb, abstol, reltol, kp1, maxc2nrm, *a_sub, lda, *done, *jbf,
                                  maxc2nrmk, relmaxc2nrmk, *jpiv_sub, *tau_sub, *vn1, *vn2, *auxv, *f_view, n_sub + nrhs,
                                  iwork, *iinfo, RT);

            if iinfo > n_sub && info.* == 0 {
                info.* = 2*ioffset + iinfo;
            }

            if done {
                k.* = ioffset + jbf;
                if iinfo <= n_sub && iinfo > 0 {
                    info.* = ioffset + iinfo;
                }
                _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, RT);
                return;
            }

            j += jbf;
        }
    }

    if j <= jmax {
        n_sub = n - j + 1;
        ioffset = j - 1;

        a_sub_raw := _geqp3rk_matrix_storage(a.*, lda, 1, j, m, n_sub + nrhs);
        a_sub := array2d_view(a_sub_raw, lda, row_major=false);
        jpiv_sub := _geqp3rk_subvector(jpiv.*, j, n_sub);
        tau_sub := _geqp3rk_subvector(tau.*, j, n_sub);
        vn1 := _geqp3rk_subvector(work.*, j, n_sub);
        vn2 := _geqp3rk_subvector(work.*, n + j, n_sub);
        work_apply := _geqp3rk_subvector(work.*, 2*n + 1, max(0, n_sub + nrhs - 1));

        _geqp3rk_laqp2rk_real(m, n_sub, nrhs, ioffset, jmax - j + 1, abstol, reltol, kp1, maxc2nrm, *a_sub, lda, *kf,
                              maxc2nrmk, relmaxc2nrmk, *jpiv_sub, *tau_sub, *vn1, *vn2, *work_apply, *iinfo, RT);

        k.* = j - 1 + kf;
        if iinfo > n_sub && info.* == 0 {
            info.* = 2*ioffset + iinfo;
        } else if iinfo <= n_sub && iinfo > 0 {
            info.* = ioffset + iinfo;
        }
    } else {
        k.* = jmax;
        if k.* < minmn {
            rem_norms := _geqp3rk_subvector(work.*, k.* + 1, n - k.*);
            jmaxc2nrm = k.* + _geqp3rk_amax_index(n - k.*, rem_norms, 1, RT);
            maxc2nrmk.* = fget(work.*, jmaxc2nrm);
            if k.* == 0 {
                relmaxc2nrmk.* = one_r;
            } else {
                relmaxc2nrmk.* = maxc2nrmk.* / maxc2nrm;
            }

            for j:k.* + 1..minmn {
                fset(tau, j, zero_r);
            }
        }
    }

    _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, RT);
}

_geqp3rk_complex :: (
    $$m: int, $$n: int, $$nrhs: int, $$kmax: int, abstol: $RT, reltol: $RT,
    a: *$A/ArrayType2D, $$lda: int, k: *int, maxc2nrmk: *RT, relmaxc2nrmk: *RT,
    jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D,
    $$lwork: int, rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           JPIV.data_type == int &&
           TAU.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 then "CGEQP3RK" else "ZGEQP3RK";
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run one_r + one_r;
    zero_mt :: #run _geqp3rk_zero_mt(RT, MT);

    lquery, done : bool = ---;
    iinfo, ioffset, j, jb, jbf, jmaxb, jmax, jmaxc2nrm, kp1, kf, lwkopt, minmn, n_sub, nb, nbmin, nx : int = ---;
    iws : int = ---;
    eps, hugeval, maxc2nrm, safmin : RT = ---;

    info.* = 0;
    lquery = lwork == -1;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if kmax < 0 {
        info.* = -4;
    } else if _geqp3rk_isnan(abstol, RT) {
        info.* = -5;
    } else if _geqp3rk_isnan(reltol, RT) {
        info.* = -6;
    } else if lda < max(1, m) {
        info.* = -8;
    }

    if info.* == 0 {
        minmn = min(m, n);
        if minmn == 0 {
            iws = 1;
            lwkopt = 1;
        } else {
            iws = n + nrhs - 1;
            nb = ilaenv(1, routine_name, " "s, m, n, -1, -1);
            lwkopt = 2*n + nb*(n + nrhs + 1);
        }
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT);
        if lwork < iws && !lquery {
            info.* = -15;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if minmn == 0 {
        k.* = 0;
        maxc2nrmk.* = zero_r;
        relmaxc2nrmk.* = zero_r;
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT);
        return;
    }

    for j:1..n {
        fset(jpiv, j, j);
    }

    for j:1..n {
        col := _geqp3rk_raw_slice(a.data + (j - 1)*lda, m);
        norm := _geqp3rk_nrm2_complex(m, col, 1, RT, MT);
        fset(rwork, j, norm);
        fset(rwork, n + j, norm);
    }

    norms := _geqp3rk_subvector(rwork.*, 1, n);
    kp1 = _geqp3rk_amax_index(n, norms, 1, RT);
    maxc2nrm = fget(rwork.*, kp1);

    if _geqp3rk_isnan(maxc2nrm, RT) {
        k.* = 0;
        info.* = kp1;
        maxc2nrmk.* = maxc2nrm;
        relmaxc2nrmk.* = maxc2nrm;
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT);
        return;
    }

    if maxc2nrm == zero_r {
        k.* = 0;
        maxc2nrmk.* = zero_r;
        relmaxc2nrmk.* = zero_r;
        for j:1..minmn {
            fset(tau, j, zero_mt);
        }
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT);
        return;
    }

    hugeval = _geqp3rk_lamch(#char "o", RT);
    if maxc2nrm > hugeval {
        info.* = n + kp1;
    }

    if kmax == 0 {
        k.* = 0;
        maxc2nrmk.* = maxc2nrm;
        relmaxc2nrmk.* = one_r;
        for j:1..minmn {
            fset(tau, j, zero_mt);
        }
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT);
        return;
    }

    eps = _geqp3rk_lamch(#char "e", RT);
    if abstol >= zero_r {
        safmin = _geqp3rk_lamch(#char "s", RT);
        abstol = max(abstol, two_r*safmin);
    }
    if reltol >= zero_r {
        reltol = max(reltol, eps);
    }

    jmax = min(kmax, minmn);
    if maxc2nrm <= abstol || one_r <= reltol {
        k.* = 0;
        maxc2nrmk.* = maxc2nrm;
        relmaxc2nrmk.* = one_r;
        for j:1..minmn {
            fset(tau, j, zero_mt);
        }
        _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT);
        return;
    }

    nbmin = 2;
    nx = 0;
    if nb > 1 && nb < minmn {
        nx = max(0, ilaenv(3, routine_name, " "s, m, n, -1, -1));
        if nx < minmn {
            if lwork < lwkopt {
                nb = (lwork - 2*n) / (n + 1);
                nbmin = max(2, ilaenv(2, routine_name, " "s, m, n, -1, -1));
            }
        }
    }

    done = false;
    j = 1;
    jmaxb = min(kmax, minmn - nx);
    if nb >= nbmin && nb < jmax && jmaxb > 0 {
        while j <= jmaxb {
            jb = min(nb, jmaxb - j + 1);
            n_sub = n - j + 1;
            ioffset = j - 1;

            a_sub_raw := _geqp3rk_matrix_storage(a.*, lda, 1, j, m, n_sub + nrhs);
            a_sub := array2d_view(a_sub_raw, lda, row_major=false);
            jpiv_sub := _geqp3rk_subvector(jpiv.*, j, n_sub);
            tau_sub := _geqp3rk_subvector(tau.*, j, n_sub);
            vn1 := _geqp3rk_subvector(rwork.*, j, n_sub);
            vn2 := _geqp3rk_subvector(rwork.*, n + j, n_sub);
            auxv := _geqp3rk_subvector(work.*, 1, jb);
            f_raw := _geqp3rk_raw_slice(work.data + jb, max(0, (n_sub + nrhs)*jb));
            f_view := array2d_view(f_raw, n_sub + nrhs, row_major=false);

            _geqp3rk_laqp3rk_complex(m, n_sub, nrhs, ioffset, jb, abstol, reltol, kp1, maxc2nrm, *a_sub, lda, *done, *jbf,
                                     maxc2nrmk, relmaxc2nrmk, *jpiv_sub, *tau_sub, *vn1, *vn2, *auxv, *f_view, n_sub + nrhs,
                                     iwork, *iinfo, RT, MT);

            if iinfo > n_sub && info.* == 0 {
                info.* = 2*ioffset + iinfo;
            }

            if done {
                k.* = ioffset + jbf;
                if iinfo <= n_sub && iinfo > 0 {
                    info.* = ioffset + iinfo;
                }
                _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT);
                return;
            }

            j += jbf;
        }
    }

    if j <= jmax {
        n_sub = n - j + 1;
        ioffset = j - 1;

        a_sub_raw := _geqp3rk_matrix_storage(a.*, lda, 1, j, m, n_sub + nrhs);
        a_sub := array2d_view(a_sub_raw, lda, row_major=false);
        jpiv_sub := _geqp3rk_subvector(jpiv.*, j, n_sub);
        tau_sub := _geqp3rk_subvector(tau.*, j, n_sub);
        vn1 := _geqp3rk_subvector(rwork.*, j, n_sub);
        vn2 := _geqp3rk_subvector(rwork.*, n + j, n_sub);
        work_apply := _geqp3rk_subvector(work.*, 1, max(0, n_sub + nrhs - 1));

        _geqp3rk_laqp2rk_complex(m, n_sub, nrhs, ioffset, jmax - j + 1, abstol, reltol, kp1, maxc2nrm, *a_sub, lda, *kf,
                                 maxc2nrmk, relmaxc2nrmk, *jpiv_sub, *tau_sub, *vn1, *vn2, *work_apply, *iinfo, RT, MT);

        k.* = j - 1 + kf;
        if iinfo > n_sub && info.* == 0 {
            info.* = 2*ioffset + iinfo;
        } else if iinfo <= n_sub && iinfo > 0 {
            info.* = ioffset + iinfo;
        }
    } else {
        k.* = jmax;
        if k.* < minmn {
            rem_norms := _geqp3rk_subvector(rwork.*, k.* + 1, n - k.*);
            jmaxc2nrm = k.* + _geqp3rk_amax_index(n - k.*, rem_norms, 1, RT);
            maxc2nrmk.* = fget(rwork.*, jmaxc2nrm);
            if k.* == 0 {
                relmaxc2nrmk.* = one_r;
            } else {
                relmaxc2nrmk.* = maxc2nrmk.* / maxc2nrm;
            }

            for j:k.* + 1..minmn {
                fset(tau, j, zero_mt);
            }
        } else {
            maxc2nrmk.* = zero_r;
            relmaxc2nrmk.* = zero_r;
        }
    }

    _geqp3rk_set_work1(work, _geqp3rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT);
}
