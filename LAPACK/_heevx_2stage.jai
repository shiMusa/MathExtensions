// ? LAPACK *heevx_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheevx_2stage/zheevx_2stage
// *HEEVX_2STAGE computes selected eigenvalues and, optionally,
// eigenvectors of a complex Hermitian matrix using the 2stage reduction.

cheevx_2stage :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: float32, vu: float32, $$il: int, $$iu: int, abstol: float32, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int
) {
    _heevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info, RT=float32, MT=C32);
}

zheevx_2stage :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: float64, vu: float64, $$il: int, $$iu: int, abstol: float64, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int
) {
    _heevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info, RT=float64, MT=C64);
}

#scope_file
_heevx_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} { return cpx(x, zero(RT)); }

_heevx_2stage_query_value :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        return cpx(sroundup_lwork(x), zero(RT));
    } else {
        return cpx(x, zero(RT));
    }
}

_heevx_2stage_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_heevx_2stage_matrix_storage :: inline (data: *$T, $$ld: int, $$rows: int, $$cols: int) -> [] T {
    count := ifx rows <= 0 || cols <= 0 then 0 else ld*(cols - 1) + rows;
    return _heevx_2stage_raw_slice(data, count);
}

_heevx_2stage_lanhe :: (
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

_heevx_2stage_scal_complex_real :: (
    $$n: int, $$alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_heevx_2stage_hetrd_2stage :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D,
    tau: *$TAU/ArrayType1D, hous: *$HOUS/ArrayType1D, $$lhous: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
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

_heevx_2stage_ungtr :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        cungtr(uplo, n, a, lda, tau, work, lwork, info);
    } else {
        zungtr(uplo, n, a, lda, tau, work, lwork, info);
    }
}

_heevx_2stage_steqr :: (
    $$compz: u8, $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && Z.data_type == MT && WORK.data_type == RT;
} {
    #if MT == C32 {
        csteqr(compz, n, d, e, z, ldz, work, info);
    } else {
        zsteqr(compz, n, d, e, z, ldz, work, info);
    }
}

_heevx_2stage_stebz :: (
    $$range: u8, $$order: u8, $$n: int, $$vl: $RT, $$vu: RT, $$il: int, $$iu: int, $$abstol: RT,
    d: *$D/ArrayType1D, e: *$E/ArrayType1D, m: *int, nsplit: *int, w: *$W/ArrayType1D,
    iblock: *$IBLOCK/ArrayType1D, isplit: *$ISPLIT/ArrayType1D, work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           D.data_type == RT && E.data_type == RT && W.data_type == RT &&
           IBLOCK.data_type == int && ISPLIT.data_type == int &&
           WORK.data_type == RT && IWORK.data_type == int;
} {
    #if RT == float32 {
        sstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info);
    } else {
        dstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info);
    }
}

_heevx_2stage_stein :: (
    $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, $$m: int, w: *$W/ArrayType1D,
    iblock: *$IBLOCK/ArrayType1D, isplit: *$ISPLIT/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && W.data_type == RT &&
           IBLOCK.data_type == int && ISPLIT.data_type == int &&
           Z.data_type == MT && WORK.data_type == RT && IWORK.data_type == int &&
           IFAIL.data_type == int;
} {
    #if MT == C32 {
        cstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info);
    } else {
        zstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info);
    }
}

_heevx_2stage_unmtr :: (
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

_heevx_2stage_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        cswap(n, x, incx, y, incy);
    } else {
        zswap(n, x, incx, y, incy);
    }
}

_heevx_2stage_copy_real :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 {
        scopy(n, x, incx, y, incy);
    } else {
        dcopy(n, x, incx, y, incy);
    }
}

_heevx_2stage_lacpy :: (
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

_heevx_2stage_scal_real :: (
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

_heevx_2stage_sterf :: (
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

_heevx_2stage :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: $RT, vu: RT, $$il: int, $$iu: int, abstol: RT, m: *int, w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int && IFAIL.data_type == int;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    cone :: #run _heevx_2stage_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CHEEVX_2STAGE" } else { "ZHEEVX_2STAGE" };
    wantz := lsame(jobz, #char "V");
    lower := lsame(uplo, #char "L");
    alleig := lsame(range, #char "A");
    valeig := lsame(range, #char "V");
    indeig := lsame(range, #char "I");
    test, use_fallback : bool = ---;
    order : u8 = ---;
    lquery := lwork == -1;
    abstll, anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum, tmp1, vll, vuu : RT = ---;
    i, iinfo, imax, indd, inde, indee, indibl, indisp, indiwk, indrwk, indtau, indwrk, iscale,
    itmp1, j, jj, kd, ib, lhtrd, lwtrd, lwmin, llwork, nsplit, indhous : int = ---;

    info.* = 0;
    if !lsame(jobz, #char "N") { info.* = -1; }
    else if !(alleig || valeig || indeig) { info.* = -2; }
    else if !(lower || lsame(uplo, #char "U")) { info.* = -3; }
    else if n < 0 { info.* = -4; }
    else if lda < max(1, n) { info.* = -6; }
    else {
        if valeig {
            if n > 0 && vu <= vl { info.* = -8; }
        } else if indeig {
            if il < 1 || il > max(1, n) { info.* = -9; }
            else if iu < min(n, il) || iu > n { info.* = -10; }
        }
    }
    if info.* == 0 && (ldz < 1 || (wantz && ldz < n)) { info.* = -15; }
    if info.* == 0 {
        if n <= 1 {
            lwmin = 1;
        } else {
            kd = ilaenv2stage(1, "CHETRD_2STAGE", jobz, n, -1, -1, -1);
            ib = ilaenv2stage(2, "CHETRD_2STAGE", jobz, n, kd, -1, -1);
            lhtrd = ilaenv2stage(3, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
            lwtrd = ilaenv2stage(4, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
            lwmin = n + lhtrd + lwtrd;
        }
        fset(work, 1, _heevx_2stage_query_value(ncast(RT, lwmin), RT, MT));
        if lwork < lwmin && !lquery { info.* = -17; }
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }
    m.* = 0;
    if n == 0 then return;
    if n == 1 {
        if alleig || indeig {
            m.* = 1;
            fset(w, 1, real(fget(a.*, 1, 1)));
        } else if valeig {
            if vl < real(fget(a.*, 1, 1)) && vu >= real(fget(a.*, 1, 1)) {
                m.* = 1;
                fset(w, 1, real(fget(a.*, 1, 1)));
            }
        }
        if wantz {
            fset(z, 1, 1, cone);
        }
        return;
    }

    safmin = ifx RT == float32 then slamch("Safe minimum") else dlamch("Safe minimum");
    eps = ifx RT == float32 then slamch("Precision") else dlamch("Precision");
    smlnum = safmin / eps;
    bignum = one_r / smlnum;
    rmin = sqrt(smlnum);
    rmax = min(sqrt(bignum), one_r / sqrt(sqrt(safmin)));

    iscale = 0;
    abstll = abstol;
    if valeig {
        vll = vl;
        vuu = vu;
    }
    anrm = _heevx_2stage_lanhe(#char "M", uplo, n, a, lda, rwork, RT, MT);
    if anrm > zero_r && anrm < rmin {
        iscale = 1;
        sigma = rmin / anrm;
    } else if anrm > rmax {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if iscale == 1 {
        if lower {
            for j:1..n {
                col := _heevx_2stage_raw_slice(a.data + (j - 1)*lda + (j - 1), n - j + 1);
                _heevx_2stage_scal_complex_real(n - j + 1, sigma, *col, 1, RT, MT);
            }
        } else {
            for j:1..n {
                col := _heevx_2stage_raw_slice(a.data + (j - 1)*lda, j);
                _heevx_2stage_scal_complex_real(j, sigma, *col, 1, RT, MT);
            }
        }
        if abstol > zero_r {
            abstll = abstol * sigma;
        }
        if valeig {
            vll = vl * sigma;
            vuu = vu * sigma;
        }
    }

    indd = 1;
    inde = indd + n;
    indrwk = inde + n;
    indtau = 1;
    indhous = indtau + n;
    indwrk = indhous + lhtrd;
    llwork = lwork - indwrk + 1;
    d_view := _heevx_2stage_raw_slice(rwork.data + (indd - 1), n);
    e_view := _heevx_2stage_raw_slice(rwork.data + (inde - 1), n);
    tau_view := _heevx_2stage_raw_slice(work.data + (indtau - 1), n);
    hous_view := _heevx_2stage_raw_slice(work.data + (indhous - 1), lhtrd);
    work_view := _heevx_2stage_raw_slice(work.data + (indwrk - 1), max(0, llwork));
    _heevx_2stage_hetrd_2stage(jobz, uplo, n, a, lda, *d_view, *e_view, *tau_view, *hous_view, lhtrd, *work_view, llwork, *iinfo, RT, MT);

    test = false;
    if indeig {
        if il == 1 && iu == n {
            test = true;
        }
    }
    if (alleig || test) && (abstol <= zero_r) {
        _heevx_2stage_copy_real(n, *d_view, 1, w, 1, RT);
        indee = indrwk + 2*n;
        ee_view := _heevx_2stage_raw_slice(rwork.data + (indee - 1), max(0, n - 1));
        if !wantz {
            _heevx_2stage_copy_real(n - 1, *e_view, 1, *ee_view, 1, RT);
            _heevx_2stage_sterf(n, w, *ee_view, info, RT);
        } else {
            z_full_raw := _heevx_2stage_matrix_storage(a.data, lda, n, n);
            z_full := array2d_view(z_full_raw, lda, row_major=false);
            z_dst_raw := _heevx_2stage_matrix_storage(z.data, ldz, n, n);
            z_dst := array2d_view(z_dst_raw, ldz, row_major=false);
            _heevx_2stage_lacpy(#char "A", n, n, *z_full, lda, *z_dst, ldz, MT);
            _heevx_2stage_ungtr(uplo, n, z, ldz, *tau_view, *work_view, llwork, *iinfo, RT, MT);
            _heevx_2stage_copy_real(n - 1, *e_view, 1, *ee_view, 1, RT);
            _heevx_2stage_steqr(jobz, n, w, *ee_view, z, ldz, rwork, info, RT, MT);
            if info.* == 0 {
                for i:1..n {
                    fset(ifail, i, 0);
                }
            }
        }
        use_fallback = true;
        if info.* == 0 {
            m.* = n;
            use_fallback = false;
        } else {
            info.* = 0;
        }
    } else {
        use_fallback = true;
    }

    if use_fallback {
        if wantz { order = #char "B"; } else { order = #char "E"; }
        indibl = 1;
        indisp = indibl + n;
        indiwk = indisp + n;
        rw_view := _heevx_2stage_raw_slice(rwork.data + (indrwk - 1), max(0, 7*n - indrwk + 1));
        iw_view := _heevx_2stage_raw_slice(iwork.data + (indiwk - 1), max(0, 5*n - indiwk + 1));
        iblock_view := _heevx_2stage_raw_slice(iwork.data + (indibl - 1), n);
        isplit_view := _heevx_2stage_raw_slice(iwork.data + (indisp - 1), n);
        _heevx_2stage_stebz(range, order, n, vll, vuu, il, iu, abstll, *d_view, *e_view, m, *nsplit, w, *iblock_view, *isplit_view, *rw_view, *iw_view, info, RT);

        if wantz {
            _heevx_2stage_stein(n, *d_view, *e_view, m.*, w, *iblock_view, *isplit_view, z, ldz, *rw_view, *iw_view, ifail, info, RT, MT);
            _heevx_2stage_unmtr(#char "L", uplo, #char "N", n, m.*, a, lda, *tau_view, z, ldz, *work_view, llwork, *iinfo, RT, MT);
        }
    }

    if iscale == 1 {
        if info.* == 0 { imax = m.*; } else { imax = info.* - 1; }
        _heevx_2stage_scal_real(imax, one_r / sigma, w, 1, RT);
    }

    if wantz {
        for j:1..m.* - 1 {
            i = 0;
            tmp1 = fget(w, j);
            for jj:j + 1..m.* {
                if fget(w, jj) < tmp1 {
                    i = jj;
                    tmp1 = fget(w, jj);
                }
            }
            if i != 0 {
                itmp1 = fget(iwork, indibl + i - 1);
                fset(w, i, fget(w, j));
                fset(iwork, indibl + i - 1, fget(iwork, indibl + j - 1));
                fset(w, j, tmp1);
                fset(iwork, indibl + j - 1, itmp1);
                zcoli := _heevx_2stage_raw_slice(z.data + (i - 1)*ldz, n);
                zcolj := _heevx_2stage_raw_slice(z.data + (j - 1)*ldz, n);
                _heevx_2stage_swap(n, *zcoli, 1, *zcolj, 1, MT);
                if info.* != 0 {
                    itmp1 = fget(ifail, i);
                    fset(ifail, i, fget(ifail, j));
                    fset(ifail, j, itmp1);
                }
            }
        }
    }

    fset(work, 1, _heevx_2stage_query_value(ncast(RT, lwmin), RT, MT));
}
