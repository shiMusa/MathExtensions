// ? LAPACK *heev_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheev_2stage/zheev_2stage
// *HEEV_2STAGE computes all eigenvalues and, optionally, eigenvectors of
// a complex Hermitian matrix using the 2stage reduction.

cheev_2stage :: #bake_arguments _heev_2stage(RT=float32, MT=C32);
zheev_2stage :: #bake_arguments _heev_2stage(RT=float64, MT=C64);

#scope_file
_heev_2stage_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_heev_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_heev_2stage_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_heev_2stage_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(cmach.data); } else { return dlamch(cmach.data); }
}

_heev_2stage_lanhe :: (
    $$norm: u8, $$uplo: u8, $$n: int, a: $A/ArrayType2D, $$lda: int,
    work: *$RWORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { return clanhe(norm, uplo, n, a, lda, work); } else { return zlanhe(norm, uplo, n, a, lda, work); }
}

_heev_2stage_lascl :: (
    $$type: u8, $$kl: int, $$ku: int, cfrom: $RT, cto: $RT, $$m: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && A.data_type == MT;
} {
    #if MT == C32 { clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info); } else { zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info); }
}

_heev_2stage_scal :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 { sscal(n, alpha, x, incx); } else { dscal(n, alpha, x, incx); }
}

_heev_2stage_sterf :: (
    $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && D.data_type == RT && E.data_type == RT;
} {
    #if RT == float32 { ssterf(n, d, e, info); } else { dsterf(n, d, e, info); }
}

_heev_2stage_steqr :: (
    $$compz: u8, $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, z: *$Z/ArrayType2D,
    $$ldz: int, work: *$RWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && Z.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { csteqr(compz, n, d, e, z, ldz, work, info); } else { zsteqr(compz, n, d, e, z, ldz, work, info); }
}

_heev_2stage_ungtr :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 { cungtr(uplo, n, a, lda, tau, work, lwork, info); } else { zungtr(uplo, n, a, lda, tau, work, lwork, info); }
}

_heev_2stage_hetrd_2stage :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    d: *$D/ArrayType1D, e: *$E/ArrayType1D, tau: *$TAU/ArrayType1D,
    hous: *$HOUS/ArrayType1D, $$lhous: int, work: *$WORK/ArrayType1D,
    $$lwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && D.data_type == RT && E.data_type == RT &&
           TAU.data_type == MT && HOUS.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        chetrd_2stage(jobz, uplo, n, a, lda, d, e, tau, hous, lhous, work, lwork, info);
    } else {
        zhetrd_2stage(jobz, uplo, n, a, lda, d, e, tau, hous, lhous, work, lwork, info);
    }
}

_heev_2stage :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && W.data_type == RT && WORK.data_type == MT && RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    one_t  :: #run _heev_2stage_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHEEV_2STAGE " } else { "ZHEEV_2STAGE " };

    lower, lquery, wantz : bool = ---;
    iinfo, imax, inde, indtau, indwrk, iscale, llwork, lwmin, lhtrd, lwtrd, kd, ib, indhous : int = ---;
    anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum : RT = ---;

    wantz = lsame(jobz, #char "V");
    lower = lsame(uplo, #char "L");
    lquery = lwork == -1;

    info.* = 0;
    if !lsame(jobz, #char "N") {
        info.* = -1;
    } else if !(lower || lsame(uplo, #char "U")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    }

    if info.* == 0 {
        kd    = ilaenv2stage(1, "CHETRD_2STAGE", jobz, n, -1, -1, -1);
        ib    = ilaenv2stage(2, "CHETRD_2STAGE", jobz, n, kd, -1, -1);
        lhtrd = ilaenv2stage(3, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
        lwtrd = ilaenv2stage(4, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
        lwmin = n + lhtrd + lwtrd;
        fset(work, 1, _heev_2stage_complex_from_real(_heev_2stage_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
        if lwork < lwmin && !lquery {
            info.* = -8;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 then return;

    if n == 1 {
        fset(w, 1, real(fget(a.*, 1, 1)));
        fset(work, 1, one_t);
        if wantz { fset(a, 1, 1, one_t); }
        return;
    }

    safmin = _heev_2stage_lamch("Safe minimum"s, RT);
    eps    = _heev_2stage_lamch("Precision"s, RT);
    smlnum = safmin / eps;
    bignum = one_r / smlnum;
    rmin   = sqrt(smlnum);
    rmax   = sqrt(bignum);

    anrm = _heev_2stage_lanhe(#char "M", uplo, n, a.*, lda, rwork, RT, MT);
    iscale = 0;
    if anrm > zero_r && anrm < rmin {
        iscale = 1; sigma = rmin / anrm;
    } else if anrm > rmax {
        iscale = 1; sigma = rmax / anrm;
    }
    if iscale == 1 {
        _heev_2stage_lascl(uplo, 0, 0, one_r, sigma, n, n, a, lda, info, RT, MT);
    }

    inde = 1;
    indtau = 1;
    indhous = indtau + n;
    indwrk = indhous + lhtrd;
    llwork = lwork - indwrk + 1;
    tau_view := _heev_2stage_subvector(work.*, indtau, n);
    hous_view := _heev_2stage_subvector(work.*, indhous, lhtrd);
    work_view := _heev_2stage_subvector(work.*, indwrk, max(0, llwork));
    d_view := _heev_2stage_subvector(w.*, 1, n);
    e_view := _heev_2stage_subvector(rwork.*, inde, max(0, n - 1));
    _heev_2stage_hetrd_2stage(jobz, uplo, n, a, lda, *d_view, *e_view, *tau_view, *hous_view,
                              lhtrd, *work_view, llwork, info, RT, MT);

    if !wantz {
        _heev_2stage_sterf(n, *d_view, *e_view, info, RT);
    } else {
        _heev_2stage_ungtr(uplo, n, a, lda, *tau_view, *work_view, llwork, *iinfo, MT);
        rwrk_tail := _heev_2stage_subvector(rwork.*, inde + n, max(0, n - 1));
        _heev_2stage_steqr(jobz, n, *d_view, *e_view, a, lda, *rwrk_tail, info, RT, MT);
    }

    if iscale == 1 {
        if info.* == 0 { imax = n; } else { imax = info.* - 1; }
        _heev_2stage_scal(imax, one_r / sigma, w, 1, RT);
    }

    fset(work, 1, _heev_2stage_complex_from_real(_heev_2stage_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
}
