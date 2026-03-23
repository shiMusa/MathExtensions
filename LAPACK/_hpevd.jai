// ? LAPACK *hpevd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chpevd/zhpevd
// *HPEVD computes all eigenvalues and, optionally, eigenvectors of a Hermitian packed matrix.

chpevd :: #bake_arguments _hpevd(RT=float32, MT=C32);
zhpevd :: #bake_arguments _hpevd(RT=float64, MT=C64);

#scope_file
_hpevd_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hpevd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hpevd_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hpevd_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(cmach.data); } else { return dlamch(cmach.data); }
}

_hpevd_lanhp :: (
    $$norm: u8, $$uplo: u8, $$n: int, ap: $AP/ArrayType1D, work: *$RWORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { return clanhp(norm, uplo, n, ap, work); } else { return zlanhp(norm, uplo, n, ap, work); }
}

_hpevd_ptrd :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, d: *$D/ArrayType1D, e: *$E/ArrayType1D,
    tau: *$TAU/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && D.data_type == RT && E.data_type == RT && TAU.data_type == MT;
} {
    #if MT == C32 { chptrd(uplo, n, ap, d, e, tau, info); } else { zhptrd(uplo, n, ap, d, e, tau, info); }
}

_hpevd_scal_packed :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_hpevd_stedc :: (
    $$compz: u8, $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int;
} {
    #if MT == C32 {
        cstedc(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
    } else {
        zstedc(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
    }
}

_hpevd_upmtr :: (
    $$side: u8, $$uplo: u8, $$trans: u8, $$m: int, $$n: int, ap: $AP/ArrayType1D, tau: $TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           AP.data_type == MT && TAU.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        cupmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
    } else {
        zupmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
    }
}

_hpevd_sterf :: (
    $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && D.data_type == RT && E.data_type == RT;
} {
    #if RT == float32 { ssterf(n, d, e, info); } else { dsterf(n, d, e, info); }
}

_hpevd_scal_real :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 { sscal(n, alpha, x, incx); } else { dscal(n, alpha, x, incx); }
}

_hpevd :: (
    $$jobz: u8, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, $$lrwork: int, iwork: *$IWORK/ArrayType1D, $$liwork: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    cone   :: #run _hpevd_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHPEVD" } else { "ZHPEVD" };

    lquery, wantz : bool = ---;
    iinfo, imax, inde, inddrwk, indtau, indwrk, iscale, liwmin, llrwk, llwrk, lrwmin, lwmin : int = ---;
    anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum : RT = ---;

    wantz = lsame(jobz, #char "V");
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) {
        info.* = -1;
    } else if !(lsame(uplo, #char "L") || lsame(uplo, #char "U")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ldz < 1 || (wantz && ldz < n) {
        info.* = -7;
    }

    if info.* == 0 {
        if n <= 1 {
            lwmin = 1; liwmin = 1; lrwmin = 1;
        } else if wantz {
            lwmin = 2*n;
            lrwmin = 1 + 5*n + 2*n*n;
            liwmin = 3 + 5*n;
        } else {
            lwmin = n;
            lrwmin = n;
            liwmin = 1;
        }
        fset(work, 1, _hpevd_complex_from_real(_hpevd_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
        fset(rwork, 1, ncast(RT, lrwmin));
        fset(iwork, 1, liwmin);
        if lwork < lwmin && !lquery {
            info.* = -9;
        } else if lrwork < lrwmin && !lquery {
            info.* = -11;
        } else if liwork < liwmin && !lquery {
            info.* = -13;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 { return; }
    if n == 1 {
        fset(w, 1, real(fget(ap.*, 1)));
        if wantz { fset(z, 1, 1, cone); }
        return;
    }

    safmin = _hpevd_lamch("Safe minimum", RT);
    eps = _hpevd_lamch("Precision", RT);
    smlnum = safmin / eps;
    bignum = one_r / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

    anrm = _hpevd_lanhp(#char "M", uplo, n, ap.*, rwork, RT, MT);
    iscale = 0;
    if anrm > zero_r && anrm < rmin {
        iscale = 1;
        sigma = rmin / anrm;
    } else if anrm > rmax {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if iscale == 1 {
        _hpevd_scal_packed((n*(n + 1)) / 2, sigma, ap, 1, RT, MT);
    }

    inde = 1;
    indtau = 1;
    inddrwk = inde + n;
    indwrk = indtau + n;
    llwrk = lwork - indwrk + 1;
    llrwk = lrwork - inddrwk + 1;
    e_view := _hpevd_subvector(rwork.*, inde, n - 1);
    tau_view := _hpevd_subvector(work.*, indtau, n - 1);
    work_dc := _hpevd_subvector(work.*, indwrk, max(1, llwrk));
    rwork_dc := _hpevd_subvector(rwork.*, inddrwk, max(1, llrwk));
    _hpevd_ptrd(uplo, n, ap, w, *e_view, *tau_view, *iinfo, RT, MT);

    if !wantz {
        _hpevd_sterf(n, w, *e_view, info, RT);
    } else {
        _hpevd_stedc(#char "I", n, w, *e_view, z, ldz, *work_dc, llwrk, *rwork_dc, llrwk, iwork, liwork, info, RT, MT);
        _hpevd_upmtr(#char "L", uplo, #char "N", n, n, ap.*, tau_view, z, ldz, *work_dc, *iinfo, MT);
    }

    if iscale == 1 {
        if info.* == 0 { imax = n; } else { imax = info.* - 1; }
        _hpevd_scal_real(imax, one_r / sigma, w, 1, RT);
    }

    fset(work, 1, _hpevd_complex_from_real(_hpevd_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
    fset(rwork, 1, ncast(RT, lrwmin));
    fset(iwork, 1, liwmin);
}
