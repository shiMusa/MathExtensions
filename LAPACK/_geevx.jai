// ? LAPACK *geevx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeevx/cgeevx/dgeevx/zgeevx
// *GEEVX computes the eigenvalues and, optionally, the left and/or
// right eigenvectors for a general matrix, plus balancing and
// condition estimates.

sgeevx :: (
    $$balanc: u8,
    $$jobvl: u8,
    $$jobvr: u8,
    $$sense: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    ilo: *int,
    ihi: *int,
    scale: *$SCALE/ArrayType1D,
    abnrm: *float32,
    rconde: *$RCONDE/ArrayType1D,
    rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _geevx_real(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info, RT=float32);
}

dgeevx :: (
    $$balanc: u8,
    $$jobvl: u8,
    $$jobvr: u8,
    $$sense: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    ilo: *int,
    ihi: *int,
    scale: *$SCALE/ArrayType1D,
    abnrm: *float64,
    rconde: *$RCONDE/ArrayType1D,
    rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _geevx_real(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info, RT=float64);
}

cgeevx :: (
    $$balanc: u8,
    $$jobvl: u8,
    $$jobvr: u8,
    $$sense: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    w: *$W/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    ilo: *int,
    ihi: *int,
    scale: *$SCALE/ArrayType1D,
    abnrm: *float32,
    rconde: *$RCONDE/ArrayType1D,
    rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _geevx_complex(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info, RT=float32, MT=C32);
}

zgeevx :: (
    $$balanc: u8,
    $$jobvl: u8,
    $$jobvr: u8,
    $$sense: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    w: *$W/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    ilo: *int,
    ihi: *int,
    scale: *$SCALE/ArrayType1D,
    abnrm: *float64,
    rconde: *$RCONDE/ArrayType1D,
    rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _geevx_complex(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info, RT=float64, MT=C64);
}

#scope_file
_geevx_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geevx_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_geevx_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _geevx_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geevx_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_geevx_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_geevx_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_geevx_lange :: (
    $$norm: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           WORK.data_type == RT;
} {
    #if MT == float32 {
        return slange(norm, m, n, a, lda, work);
    } else if MT == float64 {
        return dlange(norm, m, n, a, lda, work);
    } else if MT == C32 {
        return clange(norm, m, n, a, lda, work);
    } else {
        return zlange(norm, m, n, a, lda, work);
    }
}

_geevx_lacpy :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == float32 {
        slacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == float64 {
        dlacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_geevx_lascl :: (
    $$type: u8,
    $$kl: int,
    $$ku: int,
    cfrom: $RT,
    cto: RT,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT;
} {
    #if MT == float32 {
        slascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else if MT == float64 {
        dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else if MT == C32 {
        clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    } else {
        zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
}

_geevx_nrm2 :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return snrm2(n, x, incx);
    } else if MT == float64 {
        return dnrm2(n, x, incx);
    } else if MT == C32 {
        return scnrm2(n, x, incx);
    } else {
        return dznrm2(n, x, incx);
    }
}

_geevx_scal_real :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sscal(n, alpha, x, incx);
    } else if MT == float64 {
        dscal(n, alpha, x, incx);
    } else if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_geevx_scale_complex :: (
    $$n: int,
    alpha: $MT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        cscal(n, alpha, x, incx);
    } else {
        zscal(n, alpha, x, incx);
    }
}

_geevx_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $RT: Type
) -> int #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return isamax(n, x, incx) + 1;
    } else {
        return idamax(n, x, incx) + 1;
    }
}

_geevx_rotg :: (a: $RT, b: RT, c: *RT, s: *RT, r: *RT, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        slartg(a, b, c, s, r);
    } else {
        dlartg(a, b, c, s, r);
    }
}

_geevx_rot :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    c: $RT,
    s: RT,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           Y.data_type == RT;
} {
    #if RT == float32 {
        srot(n, x, incx, y, incy, c, s);
    } else {
        drot(n, x, incx, y, incy, c, s);
    }
}

_geevx_lapy2 :: (x: $RT, y: RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slapy2(x, y);
    } else {
        return dlapy2(x, y);
    }
}

_geevx_scale_real_segment :: (
    x: *$X/ArrayType1D,
    $$offset: int,
    $$count: int,
    factor: $RT,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    if count <= 0 then return;
    for i:0..count-1 {
        fset(x, offset + i, factor * fget(x.*, offset + i));
    }
}

_geevx_scale_complex_segment :: (
    x: *$X/ArrayType1D,
    $$offset: int,
    $$count: int,
    factor: $RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    if count <= 0 then return;
    factor_t := cpx(factor, zero(RT));
    for i:0..count-1 {
        fset(x, offset + i, factor_t * fget(x.*, offset + i));
    }
}

_geevx_trsna :: (
    $$sense: u8,
    $$job: u8,
    select: $SELECT/ArrayType1D,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    vl: $VL/ArrayType2D,
    $$ldvl: int,
    vr: $VR/ArrayType2D,
    $$ldvr: int,
    rconde: *$RCONDE/ArrayType1D,
    rcondv: *$RCONDV/ArrayType1D,
    $$mm: int,
    m: *int,
    work: *$WORK/ArrayType1D,
    $$ldwork: int,
    aux: *$AUX/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           VL.data_type == MT &&
           VR.data_type == MT &&
           RCONDE.data_type == RT &&
           RCONDV.data_type == RT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) ||
            (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        strsna(sense, job, select, n, a, lda, vl, ldvl, vr, ldvr, rconde, rcondv, mm, m, work, ldwork, aux, info);
    } else if MT == float64 {
        dtrsna(sense, job, select, n, a, lda, vl, ldvl, vr, ldvr, rconde, rcondv, mm, m, work, ldwork, aux, info);
    } else if MT == C32 {
        ctrsna(sense, job, select, n, a, lda, vl, ldvl, vr, ldvr, rconde, rcondv, mm, m, work, ldwork, aux, info);
    } else {
        ztrsna(sense, job, select, n, a, lda, vl, ldvl, vr, ldvr, rconde, rcondv, mm, m, work, ldwork, aux, info);
    }
}

_geevx_real :: (
    $$balanc: u8,
    $$jobvl: u8,
    $$jobvr: u8,
    $$sense: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    ilo: *int,
    ihi: *int,
    scale: *$SCALE/ArrayType1D,
    abnrm: *RT,
    rconde: *$RCONDE/ArrayType1D,
    rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           WR.data_type == RT &&
           WI.data_type == RT &&
           VL.data_type == RT &&
           VR.data_type == RT &&
           SCALE.data_type == RT &&
           RCONDE.data_type == RT &&
           RCONDV.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if RT == float32 { "SGEEVX" } else { "DGEEVX" };
    gehrd_name :: #run if RT == float32 { "SGEHRD" } else { "DGEHRD" };
    orghr_name :: #run if RT == float32 { "SORGHR" } else { "DORGHR" };

    wantvl, wantvr, wntsnn, wntsne, wntsnv, wntsnb, lquery, scalea : bool = ---;
    job, side : u8 = ---;
    hswork, i, icond, ierr, itau, iwrk, k, lwork_trevc, maxwrk, minwrk, nout : int = ---;
    anrm, bignum, cs, cscale, eps, r, scl, smlnum, sn : RT = ---;
    select_dummy : [1] bool;
    dum : [1] RT;
    hseqr_query : [1] RT;
    trevc_query : [1] RT;

    info.* = 0;
    lquery = lwork == -1;
    wantvl = lsame(jobvl, #char "V");
    wantvr = lsame(jobvr, #char "V");
    wntsnn = lsame(sense, #char "N");
    wntsne = lsame(sense, #char "E");
    wntsnv = lsame(sense, #char "V");
    wntsnb = lsame(sense, #char "B");

    if !(lsame(balanc, #char "N") || lsame(balanc, #char "S") || lsame(balanc, #char "P") || lsame(balanc, #char "B")) {
        info.* = -1;
    } else if !wantvl && !lsame(jobvl, #char "N") {
        info.* = -2;
    } else if !wantvr && !lsame(jobvr, #char "N") {
        info.* = -3;
    } else if !(wntsnn || wntsne || wntsnb || wntsnv) || ((wntsne || wntsnb) && !(wantvl && wantvr)) {
        info.* = -4;
    } else if n < 0 {
        info.* = -5;
    } else if lda < max(1, n) {
        info.* = -7;
    } else if ldvl < 1 || (wantvl && ldvl < n) {
        info.* = -11;
    } else if ldvr < 1 || (wantvr && ldvr < n) {
        info.* = -13;
    }

    if info.* == 0 {
        if n == 0 {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = n + n*ilaenv(1, gehrd_name, " "s, n, 1, n, 0);

            if wantvl {
                #if RT == float32 {
                    strevc3(#char "L", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, *ierr);
                    shseqr(#char "S", #char "V", n, 1, n, a, lda, wr, wi, vl, ldvl, hseqr_query, -1, *info);
                } else {
                    dtrevc3(#char "L", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, *ierr);
                    dhseqr(#char "S", #char "V", n, 1, n, a, lda, wr, wi, vl, ldvl, hseqr_query, -1, *info);
                }
                lwork_trevc = cast(int) trevc_query[0];
                maxwrk = max(maxwrk, n + lwork_trevc);
            } else if wantvr {
                #if RT == float32 {
                    strevc3(#char "R", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, *ierr);
                    shseqr(#char "S", #char "V", n, 1, n, a, lda, wr, wi, vr, ldvr, hseqr_query, -1, *info);
                } else {
                    dtrevc3(#char "R", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, *ierr);
                    dhseqr(#char "S", #char "V", n, 1, n, a, lda, wr, wi, vr, ldvr, hseqr_query, -1, *info);
                }
                lwork_trevc = cast(int) trevc_query[0];
                maxwrk = max(maxwrk, n + lwork_trevc);
            } else {
                job = ifx wntsnn then #char "E" else #char "S";
                #if RT == float32 {
                    shseqr(job, #char "N", n, 1, n, a, lda, wr, wi, vr, ldvr, hseqr_query, -1, *info);
                } else {
                    dhseqr(job, #char "N", n, 1, n, a, lda, wr, wi, vr, ldvr, hseqr_query, -1, *info);
                }
            }
            hswork = cast(int) hseqr_query[0];

            if !wantvl && !wantvr {
                minwrk = 2*n;
                if !wntsnn {
                    minwrk = max(minwrk, n*n + 6*n);
                }
                maxwrk = max(maxwrk, hswork);
                if !wntsnn {
                    maxwrk = max(maxwrk, n*n + 6*n);
                }
            } else {
                minwrk = 3*n;
                if !wntsnn && !wntsne {
                    minwrk = max(minwrk, n*n + 6*n);
                }
                maxwrk = max(maxwrk, hswork);
                maxwrk = max(maxwrk, n + (n - 1)*ilaenv(1, orghr_name, " "s, n, 1, n, -1));
                if !wntsnn && !wntsne {
                    maxwrk = max(maxwrk, n*n + 6*n);
                }
                maxwrk = max(maxwrk, 3*n);
            }

            maxwrk = max(maxwrk, minwrk);
        }

        fset(work, 1, ncast(RT, maxwrk));
        if lwork < minwrk && !lquery {
            info.* = -21;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 then return;

    eps = _geevx_lamch_precision(RT);
    smlnum = _geevx_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    icond = 0;
    anrm = _geevx_lange(#char "M", n, n, a.*, lda, *dum, RT, RT);
    scalea = false;
    if anrm > zero_r && anrm < smlnum {
        scalea = true;
        cscale = smlnum;
    } else if anrm > bignum {
        scalea = true;
        cscale = bignum;
    }
    if scalea {
        _geevx_lascl(#char "G", 0, 0, anrm, cscale, n, n, a, lda, *ierr, RT, RT);
    }

    #if RT == float32 {
        sgebal(balanc, n, a, lda, ilo, ihi, scale, *ierr);
    } else {
        dgebal(balanc, n, a, lda, ilo, ihi, scale, *ierr);
    }

    abnrm.* = _geevx_lange(#char "1", n, n, a.*, lda, *dum, RT, RT);
    if scalea {
        abnrm.* *= anrm / cscale;
    }

    itau = 1;
    iwrk = itau + n;
    tau := _geevx_subvector(work.*, itau, n);
    wrk := _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if RT == float32 {
        sgehrd(n, ilo.*, ihi.*, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    } else {
        dgehrd(n, ilo.*, ihi.*, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    }

    if wantvl {
        side = #char "L";
        _geevx_lacpy(#char "L", n, n, a.*, lda, vl, ldvl, RT);
        #if RT == float32 {
            sorghr(n, ilo.*, ihi.*, vl, ldvl, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            dorghr(n, ilo.*, ihi.*, vl, ldvl, tau, wrk, lwork - iwrk + 1, *ierr);
        }
        iwrk = itau;
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if RT == float32 {
            shseqr(#char "S", #char "V", n, ilo.*, ihi.*, a, lda, wr, wi, vl, ldvl, wrk, lwork - iwrk + 1, info);
        } else {
            dhseqr(#char "S", #char "V", n, ilo.*, ihi.*, a, lda, wr, wi, vl, ldvl, wrk, lwork - iwrk + 1, info);
        }
        if wantvr {
            side = #char "B";
            _geevx_lacpy(#char "F", n, n, vl.*, ldvl, vr, ldvr, RT);
        }
    } else if wantvr {
        side = #char "R";
        _geevx_lacpy(#char "L", n, n, a.*, lda, vr, ldvr, RT);
        #if RT == float32 {
            sorghr(n, ilo.*, ihi.*, vr, ldvr, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            dorghr(n, ilo.*, ihi.*, vr, ldvr, tau, wrk, lwork - iwrk + 1, *ierr);
        }
        iwrk = itau;
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if RT == float32 {
            shseqr(#char "S", #char "V", n, ilo.*, ihi.*, a, lda, wr, wi, vr, ldvr, wrk, lwork - iwrk + 1, info);
        } else {
            dhseqr(#char "S", #char "V", n, ilo.*, ihi.*, a, lda, wr, wi, vr, ldvr, wrk, lwork - iwrk + 1, info);
        }
    } else {
        job = ifx wntsnn then #char "E" else #char "S";
        iwrk = itau;
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if RT == float32 {
            shseqr(job, #char "N", n, ilo.*, ihi.*, a, lda, wr, wi, vr, ldvr, wrk, lwork - iwrk + 1, info);
        } else {
            dhseqr(job, #char "N", n, ilo.*, ihi.*, a, lda, wr, wi, vr, ldvr, wrk, lwork - iwrk + 1, info);
        }
    }

    if info.* != 0 goto real_end_scaling;

    if wantvl || wantvr {
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if RT == float32 {
            strevc3(side, #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, wrk, lwork - iwrk + 1, *ierr);
        } else {
            dtrevc3(side, #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, wrk, lwork - iwrk + 1, *ierr);
        }
    }

    if !wntsnn {
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        _geevx_trsna(sense, #char "A", select_dummy, n, a.*, lda, vl.*, ldvl, vr.*, ldvr, rconde, rcondv, n, *nout, wrk, n, iwork, *icond, RT, RT);
    }

    if wantvl {
        #if RT == float32 {
            sgebak(balanc, #char "L", n, ilo.*, ihi.*, scale.*, n, vl, ldvl, *ierr);
        } else {
            dgebak(balanc, #char "L", n, ilo.*, ihi.*, scale.*, n, vl, ldvl, *ierr);
        }

        i = 1;
        while i <= n {
            if fget(wi.*, i) == zero_r {
                col := _geevx_matrix_column(vl.*, ldvl, 1, i, n);
                scl = one_r / _geevx_nrm2(n, col, 1, RT, RT);
                _geevx_scal_real(n, scl, *col, 1, RT, RT);
                i += 1;
            } else if fget(wi.*, i) > zero_r {
                col1 := _geevx_matrix_column(vl.*, ldvl, 1, i, n);
                col2 := _geevx_matrix_column(vl.*, ldvl, 1, i + 1, n);
                scl = one_r / _geevx_lapy2(_geevx_nrm2(n, col1, 1, RT, RT), _geevx_nrm2(n, col2, 1, RT, RT), RT);
                _geevx_scal_real(n, scl, *col1, 1, RT, RT);
                _geevx_scal_real(n, scl, *col2, 1, RT, RT);
                for k:1..n {
                    fset(work, k, fget(vl.*, k, i)*fget(vl.*, k, i) + fget(vl.*, k, i + 1)*fget(vl.*, k, i + 1));
                }
                scratch := _geevx_subvector(work.*, 1, n);
                k = _geevx_iamax_1based(n, scratch, 1, RT);
                _geevx_rotg(fget(vl.*, k, i), fget(vl.*, k, i + 1), *cs, *sn, *r, RT);
                _geevx_rot(n, *col1, 1, *col2, 1, cs, sn, RT);
                fset(vl, k, i + 1, zero_r);
                i += 2;
            } else {
                i += 1;
            }
        }
    }

    if wantvr {
        #if RT == float32 {
            sgebak(balanc, #char "R", n, ilo.*, ihi.*, scale.*, n, vr, ldvr, *ierr);
        } else {
            dgebak(balanc, #char "R", n, ilo.*, ihi.*, scale.*, n, vr, ldvr, *ierr);
        }

        i = 1;
        while i <= n {
            if fget(wi.*, i) == zero_r {
                col := _geevx_matrix_column(vr.*, ldvr, 1, i, n);
                scl = one_r / _geevx_nrm2(n, col, 1, RT, RT);
                _geevx_scal_real(n, scl, *col, 1, RT, RT);
                i += 1;
            } else if fget(wi.*, i) > zero_r {
                col1 := _geevx_matrix_column(vr.*, ldvr, 1, i, n);
                col2 := _geevx_matrix_column(vr.*, ldvr, 1, i + 1, n);
                scl = one_r / _geevx_lapy2(_geevx_nrm2(n, col1, 1, RT, RT), _geevx_nrm2(n, col2, 1, RT, RT), RT);
                _geevx_scal_real(n, scl, *col1, 1, RT, RT);
                _geevx_scal_real(n, scl, *col2, 1, RT, RT);
                for k:1..n {
                    fset(work, k, fget(vr.*, k, i)*fget(vr.*, k, i) + fget(vr.*, k, i + 1)*fget(vr.*, k, i + 1));
                }
                scratch := _geevx_subvector(work.*, 1, n);
                k = _geevx_iamax_1based(n, scratch, 1, RT);
                _geevx_rotg(fget(vr.*, k, i), fget(vr.*, k, i + 1), *cs, *sn, *r, RT);
                _geevx_rot(n, *col1, 1, *col2, 1, cs, sn, RT);
                fset(vr, k, i + 1, zero_r);
                i += 2;
            } else {
                i += 1;
            }
        }
    }

real_end_scaling:
    if scalea {
        factor := anrm / cscale;
        _geevx_scale_real_segment(wr, info.* + 1, n - info.*, factor, RT);
        _geevx_scale_real_segment(wi, info.* + 1, n - info.*, factor, RT);
        if info.* == 0 {
            if (wntsnv || wntsnb) && icond == 0 {
                _geevx_scale_real_segment(rcondv, 1, n, factor, RT);
            }
        } else {
            _geevx_scale_real_segment(wr, 1, ilo.* - 1, factor, RT);
            _geevx_scale_real_segment(wi, 1, ilo.* - 1, factor, RT);
        }
    }

    fset(work, 1, ncast(RT, maxwrk));
}

_geevx_complex :: (
    $$balanc: u8,
    $$jobvl: u8,
    $$jobvr: u8,
    $$sense: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    w: *$W/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    ilo: *int,
    ihi: *int,
    scale: *$SCALE/ArrayType1D,
    abnrm: *RT,
    rconde: *$RCONDE/ArrayType1D,
    rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           W.data_type == MT &&
           VL.data_type == MT &&
           VR.data_type == MT &&
           SCALE.data_type == RT &&
           RCONDE.data_type == RT &&
           RCONDV.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == C32 { "CGEEVX" } else { "ZGEEVX" };
    gehrd_name :: #run if MT == C32 { "CGEHRD" } else { "ZGEHRD" };
    unghr_name :: #run if MT == C32 { "CUNGHR" } else { "ZUNGHR" };

    wantvl, wantvr, wntsnn, wntsne, wntsnv, wntsnb, lquery, scalea : bool = ---;
    job, side : u8 = ---;
    hswork, i, icond, ierr, itau, iwrk, k, lwork_trevc, maxwrk, minwrk, nout : int = ---;
    anrm, bignum, cscale, eps, scl, smlnum : RT = ---;
    tmp : MT = ---;
    select_dummy : [1] bool;
    dum : [1] RT;
    hseqr_query : [1] MT;
    trevc_query : [1] MT;

    info.* = 0;
    lquery = lwork == -1;
    wantvl = lsame(jobvl, #char "V");
    wantvr = lsame(jobvr, #char "V");
    wntsnn = lsame(sense, #char "N");
    wntsne = lsame(sense, #char "E");
    wntsnv = lsame(sense, #char "V");
    wntsnb = lsame(sense, #char "B");

    if !(lsame(balanc, #char "N") || lsame(balanc, #char "S") || lsame(balanc, #char "P") || lsame(balanc, #char "B")) {
        info.* = -1;
    } else if !wantvl && !lsame(jobvl, #char "N") {
        info.* = -2;
    } else if !wantvr && !lsame(jobvr, #char "N") {
        info.* = -3;
    } else if !(wntsnn || wntsne || wntsnb || wntsnv) || ((wntsne || wntsnb) && !(wantvl && wantvr)) {
        info.* = -4;
    } else if n < 0 {
        info.* = -5;
    } else if lda < max(1, n) {
        info.* = -7;
    } else if ldvl < 1 || (wantvl && ldvl < n) {
        info.* = -10;
    } else if ldvr < 1 || (wantvr && ldvr < n) {
        info.* = -12;
    }

    if info.* == 0 {
        if n == 0 {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = n + n*ilaenv(1, gehrd_name, " "s, n, 1, n, 0);

            if wantvl {
                #if MT == C32 {
                    ctrevc3(#char "L", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, rwork, -1, *ierr);
                    chseqr(#char "S", #char "V", n, 1, n, a, lda, w, vl, ldvl, hseqr_query, -1, *info);
                } else {
                    ztrevc3(#char "L", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, rwork, -1, *ierr);
                    zhseqr(#char "S", #char "V", n, 1, n, a, lda, w, vl, ldvl, hseqr_query, -1, *info);
                }
                lwork_trevc = cast(int) real(trevc_query[0]);
                maxwrk = max(maxwrk, lwork_trevc);
            } else if wantvr {
                #if MT == C32 {
                    ctrevc3(#char "R", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, rwork, -1, *ierr);
                    chseqr(#char "S", #char "V", n, 1, n, a, lda, w, vr, ldvr, hseqr_query, -1, *info);
                } else {
                    ztrevc3(#char "R", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, rwork, -1, *ierr);
                    zhseqr(#char "S", #char "V", n, 1, n, a, lda, w, vr, ldvr, hseqr_query, -1, *info);
                }
                lwork_trevc = cast(int) real(trevc_query[0]);
                maxwrk = max(maxwrk, lwork_trevc);
            } else {
                job = ifx wntsnn then #char "E" else #char "S";
                #if MT == C32 {
                    chseqr(job, #char "N", n, 1, n, a, lda, w, vr, ldvr, hseqr_query, -1, *info);
                } else {
                    zhseqr(job, #char "N", n, 1, n, a, lda, w, vr, ldvr, hseqr_query, -1, *info);
                }
            }
            hswork = cast(int) real(hseqr_query[0]);

            if !wantvl && !wantvr {
                minwrk = 2*n;
                if !(wntsnn || wntsne) {
                    minwrk = max(minwrk, n*n + 2*n);
                }
                maxwrk = max(maxwrk, hswork);
                if !(wntsnn || wntsne) {
                    maxwrk = max(maxwrk, n*n + 2*n);
                }
            } else {
                minwrk = 2*n;
                if !(wntsnn || wntsne) {
                    minwrk = max(minwrk, n*n + 2*n);
                }
                maxwrk = max(maxwrk, hswork);
                maxwrk = max(maxwrk, n + (n - 1)*ilaenv(1, unghr_name, " "s, n, 1, n, -1));
                if !(wntsnn || wntsne) {
                    maxwrk = max(maxwrk, n*n + 2*n);
                }
                maxwrk = max(maxwrk, 2*n);
            }

            maxwrk = max(maxwrk, minwrk);
        }

        fset(work, 1, _geevx_complex_from_real(ncast(RT, maxwrk), RT, MT));
        if lwork < minwrk && !lquery {
            info.* = -20;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 then return;

    eps = _geevx_lamch_precision(RT);
    smlnum = _geevx_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    icond = 0;
    anrm = _geevx_lange(#char "M", n, n, a.*, lda, *dum, RT, MT);
    scalea = false;
    if anrm > zero_r && anrm < smlnum {
        scalea = true;
        cscale = smlnum;
    } else if anrm > bignum {
        scalea = true;
        cscale = bignum;
    }
    if scalea {
        _geevx_lascl(#char "G", 0, 0, anrm, cscale, n, n, a, lda, *ierr, RT, MT);
    }

    #if MT == C32 {
        cgebal(balanc, n, a, lda, ilo, ihi, scale, *ierr);
    } else {
        zgebal(balanc, n, a, lda, ilo, ihi, scale, *ierr);
    }

    abnrm.* = _geevx_lange(#char "1", n, n, a.*, lda, *dum, RT, MT);
    if scalea {
        abnrm.* *= anrm / cscale;
    }

    itau = 1;
    iwrk = itau + n;
    tau := _geevx_subvector(work.*, itau, n);
    wrk := _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if MT == C32 {
        cgehrd(n, ilo.*, ihi.*, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    } else {
        zgehrd(n, ilo.*, ihi.*, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    }

    if wantvl {
        side = #char "L";
        _geevx_lacpy(#char "L", n, n, a.*, lda, vl, ldvl, MT);
        #if MT == C32 {
            cunghr(n, ilo.*, ihi.*, vl, ldvl, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            zunghr(n, ilo.*, ihi.*, vl, ldvl, tau, wrk, lwork - iwrk + 1, *ierr);
        }
        iwrk = itau;
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if MT == C32 {
            chseqr(#char "S", #char "V", n, ilo.*, ihi.*, a, lda, w, vl, ldvl, wrk, lwork - iwrk + 1, info);
        } else {
            zhseqr(#char "S", #char "V", n, ilo.*, ihi.*, a, lda, w, vl, ldvl, wrk, lwork - iwrk + 1, info);
        }
        if wantvr {
            side = #char "B";
            _geevx_lacpy(#char "F", n, n, vl.*, ldvl, vr, ldvr, MT);
        }
    } else if wantvr {
        side = #char "R";
        _geevx_lacpy(#char "L", n, n, a.*, lda, vr, ldvr, MT);
        #if MT == C32 {
            cunghr(n, ilo.*, ihi.*, vr, ldvr, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            zunghr(n, ilo.*, ihi.*, vr, ldvr, tau, wrk, lwork - iwrk + 1, *ierr);
        }
        iwrk = itau;
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if MT == C32 {
            chseqr(#char "S", #char "V", n, ilo.*, ihi.*, a, lda, w, vr, ldvr, wrk, lwork - iwrk + 1, info);
        } else {
            zhseqr(#char "S", #char "V", n, ilo.*, ihi.*, a, lda, w, vr, ldvr, wrk, lwork - iwrk + 1, info);
        }
    } else {
        job = ifx wntsnn then #char "E" else #char "S";
        iwrk = itau;
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if MT == C32 {
            chseqr(job, #char "N", n, ilo.*, ihi.*, a, lda, w, vr, ldvr, wrk, lwork - iwrk + 1, info);
        } else {
            zhseqr(job, #char "N", n, ilo.*, ihi.*, a, lda, w, vr, ldvr, wrk, lwork - iwrk + 1, info);
        }
    }

    if info.* != 0 goto complex_end_scaling;

    if wantvl || wantvr {
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if MT == C32 {
            ctrevc3(side, #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, wrk, lwork - iwrk + 1, rwork, n, *ierr);
        } else {
            ztrevc3(side, #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, wrk, lwork - iwrk + 1, rwork, n, *ierr);
        }
    }

    if !wntsnn {
        wrk = _geevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        _geevx_trsna(sense, #char "A", select_dummy, n, a.*, lda, vl.*, ldvl, vr.*, ldvr, rconde, rcondv, n, *nout, wrk, n, rwork, *icond, RT, MT);
    }

    if wantvl {
        #if MT == C32 {
            cgebak(balanc, #char "L", n, ilo.*, ihi.*, scale.*, n, vl, ldvl, *ierr);
        } else {
            zgebak(balanc, #char "L", n, ilo.*, ihi.*, scale.*, n, vl, ldvl, *ierr);
        }

        for i:1..n {
            col := _geevx_matrix_column(vl.*, ldvl, 1, i, n);
            scl = one_r / _geevx_nrm2(n, col, 1, RT, MT);
            _geevx_scal_real(n, scl, *col, 1, RT, MT);
            for k:1..n {
                fset(rwork, k, real(fget(vl.*, k, i))*real(fget(vl.*, k, i)) + imag(fget(vl.*, k, i))*imag(fget(vl.*, k, i)));
            }
            scratch := _geevx_subvector(rwork.*, 1, n);
            k = _geevx_iamax_1based(n, scratch, 1, RT);
            tmp = conjugate(fget(vl.*, k, i)) / sqrt(fget(rwork.*, k));
            _geevx_scale_complex(n, tmp, *col, 1, MT);
            fset(vl, k, i, cpx(real(fget(vl.*, k, i)), zero_r));
        }
    }

    if wantvr {
        #if MT == C32 {
            cgebak(balanc, #char "R", n, ilo.*, ihi.*, scale.*, n, vr, ldvr, *ierr);
        } else {
            zgebak(balanc, #char "R", n, ilo.*, ihi.*, scale.*, n, vr, ldvr, *ierr);
        }

        for i:1..n {
            col := _geevx_matrix_column(vr.*, ldvr, 1, i, n);
            scl = one_r / _geevx_nrm2(n, col, 1, RT, MT);
            _geevx_scal_real(n, scl, *col, 1, RT, MT);
            for k:1..n {
                fset(rwork, k, real(fget(vr.*, k, i))*real(fget(vr.*, k, i)) + imag(fget(vr.*, k, i))*imag(fget(vr.*, k, i)));
            }
            scratch := _geevx_subvector(rwork.*, 1, n);
            k = _geevx_iamax_1based(n, scratch, 1, RT);
            tmp = conjugate(fget(vr.*, k, i)) / sqrt(fget(rwork.*, k));
            _geevx_scale_complex(n, tmp, *col, 1, MT);
            fset(vr, k, i, cpx(real(fget(vr.*, k, i)), zero_r));
        }
    }

complex_end_scaling:
    if scalea {
        factor := anrm / cscale;
        _geevx_scale_complex_segment(w, info.* + 1, n - info.*, factor, RT, MT);
        if info.* == 0 {
            if (wntsnv || wntsnb) && icond == 0 {
                _geevx_scale_real_segment(rcondv, 1, n, factor, RT);
            }
        } else {
            _geevx_scale_complex_segment(w, 1, ilo.* - 1, factor, RT, MT);
        }
    }

    fset(work, 1, _geevx_complex_from_real(ncast(RT, maxwrk), RT, MT));
}
