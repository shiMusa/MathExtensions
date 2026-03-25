// ? LAPACK *heevd_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheevd_2stage/zheevd_2stage
// *HEEVD_2STAGE computes all eigenvalues and, optionally, eigenvectors of
// a complex Hermitian matrix using the 2stage reduction and D&C.

cheevd_2stage :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int
) {
    _heevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float32, MT=C32);
}

zheevd_2stage :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int
) {
    _heevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float64, MT=C64);
}

#scope_file
_heevd_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} { return cpx(x, zero(RT)); }

_heevd_2stage_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_heevd_2stage_matrix_storage :: inline (data: *$T, $$ld: int, $$rows: int, $$cols: int) -> [] T {
    count := ifx rows <= 0 || cols <= 0 then 0 else ld*(cols - 1) + rows;
    return _heevd_2stage_raw_slice(data, count);
}

_heevd_2stage_lanhe :: (
    $$norm: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && A.data_type == MT && WORK.data_type == RT;
} {
    #if MT == C32 {
        return clanhe(norm, uplo, n, a, lda, work);
    } else {
        return zlanhe(norm, uplo, n, a, lda, work);
    }
}

_heevd_2stage_lascl :: (
    $$type: u8, $$kl: int, $$ku: int, $$cfrom: $RT, $$cto: RT,
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && A.data_type == MT;
} {
    #if MT == C32 {
        clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_heevd_2stage_hetrd :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D,
    tau: *$TAU/ArrayType1D, hous2: *$HOUS2/ArrayType1D, $$lhous2: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && D.data_type == RT && E.data_type == RT &&
           TAU.data_type == MT && HOUS2.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        chetrd_2stage(jobz, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);
    } else {
        zhetrd_2stage(jobz, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);
    }
}

_heevd_2stage_stedc :: (
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

_heevd_2stage_unmtr :: (
    $$side: u8, $$uplo: u8, $$trans: u8, $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && TAU.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        cunmtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_heevd_2stage_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_heevd_2stage_sterf :: (
    $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && D.data_type == RT && E.data_type == RT;
} {
    #if RT == float32 {
        ssterf(n, d, e, info);
    } else {
        dsterf(n, d, e, info);
    }
}

_heevd_2stage_scal :: (
    $$n: int, $$alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        sscal(n, alpha, x, incx);
    } else {
        dscal(n, alpha, x, incx);
    }
}

_heevd_2stage :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && W.data_type == RT && WORK.data_type == MT &&
           RWORK.data_type == RT && IWORK.data_type == int;
} {
    one_r :: #run one(RT);
    one_t :: #run _heevd_2stage_complex_from_real(one_r, RT, MT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == C32 { "CHEEVD_2STAGE" } else { "ZHEEVD_2STAGE" };
    wantz := lsame(jobz, #char "V");
    lower := lsame(uplo, #char "L");
    lquery := lwork == -1 || lrwork == -1 || liwork == -1;
    kd, ib, lhtrd, lwtrd, lwmin, lrwmin, liwmin : int = ---;
    anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum : RT = ---;
    iinfo, imax, inde, indrwk, indtau, indwk2, indwrk, iscale, llrwk, llwork, llwrk2, indhous : int = ---;

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
        if n <= 1 {
            lwmin = 1; lrwmin = 1; liwmin = 1;
        } else {
            kd = ilaenv2stage(1, "CHETRD_2STAGE", jobz, n, -1, -1, -1);
            ib = ilaenv2stage(2, "CHETRD_2STAGE", jobz, n, kd, -1, -1);
            lhtrd = ilaenv2stage(3, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
            lwtrd = ilaenv2stage(4, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
            if wantz {
                lwmin = 2*n + n*n;
                lrwmin = 1 + 5*n + 2*n*n;
                liwmin = 3 + 5*n;
            } else {
                lwmin = max(1, lhtrd + lwtrd + n + 1);
                lrwmin = n;
                liwmin = 1;
            }
        }
        fset(work, 1, _heevd_2stage_complex_from_real(ncast(RT, lwmin), RT, MT));
        fset(rwork, 1, ncast(RT, lrwmin));
        fset(iwork, 1, liwmin);
        if lwork < lwmin && !lquery { info.* = -8; }
        else if lrwork < lrwmin && !lquery { info.* = -10; }
        else if liwork < liwmin && !lquery { info.* = -12; }
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        return;
    }
    if n == 1 {
        fset(w, 1, real(fget(a.*, 1, 1)));
        if wantz {
            fset(a, 1, 1, one_t);
        }
        return;
    }

    safmin = ifx RT == float32 then slamch("Safe minimum") else dlamch("Safe minimum");
    eps    = ifx RT == float32 then slamch("Precision") else dlamch("Precision");
    smlnum = safmin / eps;
    bignum = one_r / smlnum;
    rmin   = sqrt(smlnum);
    rmax   = sqrt(bignum);

    anrm = _heevd_2stage_lanhe(#char "M", uplo, n, a, lda, rwork, RT, MT);
    iscale = 0;
    if anrm > zero_r && anrm < rmin {
        iscale = 1;
        sigma = rmin / anrm;
    } else if anrm > rmax {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if iscale == 1 {
        _heevd_2stage_lascl(uplo, 0, 0, one_r, sigma, n, n, a, lda, info, RT, MT);
    }

    inde = 1;
    indrwk = inde + n;
    llrwk = lrwork - indrwk + 1;
    indtau = 1;
    indhous = indtau + n;
    indwrk = indhous + lhtrd;
    llwork = lwork - indwrk + 1;
    indwk2 = indwrk + n*n;
    llwrk2 = lwork - indwk2 + 1;

    e_view := _heevd_2stage_raw_slice(rwork.data + (inde - 1), n);
    tau_view := _heevd_2stage_raw_slice(work.data + (indtau - 1), n);
    hous_view := _heevd_2stage_raw_slice(work.data + (indhous - 1), max(0, lhtrd));
    work1_view := _heevd_2stage_raw_slice(work.data + (indwrk - 1), max(0, llwork));
    _heevd_2stage_hetrd(jobz, uplo, n, a, lda, w, *e_view, *tau_view, *hous_view, lhtrd, *work1_view, llwork, info, RT, MT);

    if !wantz {
        _heevd_2stage_sterf(n, w, *e_view, info, RT);
    } else {
        zmat_raw := _heevd_2stage_matrix_storage(work.data + (indwrk - 1), n, n, n);
        zmat := array2d_view(zmat_raw, n, row_major=false);
        work2_view := _heevd_2stage_raw_slice(work.data + (indwk2 - 1), max(0, llwrk2));
        rwork2_view := _heevd_2stage_raw_slice(rwork.data + (indrwk - 1), max(0, llrwk));
        _heevd_2stage_stedc(#char "I", n, w, *e_view, *zmat, n, *work2_view, llwrk2, *rwork2_view, llrwk, iwork, liwork, info, RT, MT);
        _heevd_2stage_unmtr(#char "L", uplo, #char "N", n, n, a, lda, *tau_view, *zmat, n, *work2_view, llwrk2, info, RT, MT);
        a_full_raw := _heevd_2stage_matrix_storage(a.data, lda, n, n);
        a_full := array2d_view(a_full_raw, lda, row_major=false);
        _heevd_2stage_lacpy(#char "A", n, n, *zmat, n, *a_full, lda, MT);
    }

    if iscale == 1 {
        if info.* == 0 {
            imax = n;
        } else {
            imax = info.* - 1;
        }
        _heevd_2stage_scal(imax, one_r / sigma, w, 1, RT);
    }

    fset(work, 1, _heevd_2stage_complex_from_real(ncast(RT, lwmin), RT, MT));
    fset(rwork, 1, ncast(RT, lrwmin));
    fset(iwork, 1, liwmin);
}
