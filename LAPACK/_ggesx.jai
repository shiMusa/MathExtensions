// ? LAPACK *ggesx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggesx/cggesx/dggesx/zggesx
// *GGESX computes the generalized eigenvalues, generalized Schur form,
// optional generalized Schur vectors, and reciprocal condition numbers
// for selected eigenvalues and/or deflating subspaces.

sggesx :: (
    $$jobvsl: u8, $$jobvsr: u8, $$sort: u8,
    selctg: (float32, float32, float32) -> bool,
    $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    sdim: *int,
    alphar: *$ALPHAR/ArrayType1D, alphai: *$ALPHAI/ArrayType1D, beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D, $$ldvsl: int,
    vsr: *$VSR/ArrayType2D, $$ldvsr: int,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int,
    bwork: *$BWORK/ArrayType1D, info: *int
) {
    _ggesx_real(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, bwork, info, RT=float32);
}

dggesx :: (
    $$jobvsl: u8, $$jobvsr: u8, $$sort: u8,
    selctg: (float64, float64, float64) -> bool,
    $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    sdim: *int,
    alphar: *$ALPHAR/ArrayType1D, alphai: *$ALPHAI/ArrayType1D, beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D, $$ldvsl: int,
    vsr: *$VSR/ArrayType2D, $$ldvsr: int,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int,
    bwork: *$BWORK/ArrayType1D, info: *int
) {
    _ggesx_real(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, bwork, info, RT=float64);
}

cggesx :: (
    $$jobvsl: u8, $$jobvsr: u8, $$sort: u8,
    selctg: (C32, C32) -> bool,
    $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    sdim: *int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D, $$ldvsl: int,
    vsr: *$VSR/ArrayType2D, $$ldvsr: int,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, $$liwork: int,
    bwork: *$BWORK/ArrayType1D, info: *int
) {
    _ggesx_complex(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork, bwork, info, RT=float32, MT=C32);
}

zggesx :: (
    $$jobvsl: u8, $$jobvsr: u8, $$sort: u8,
    selctg: (C64, C64) -> bool,
    $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    sdim: *int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D, $$ldvsl: int,
    vsr: *$VSR/ArrayType2D, $$ldvsr: int,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, $$liwork: int,
    bwork: *$BWORK/ArrayType1D, info: *int
) {
    _ggesx_complex(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork, bwork, info, RT=float64, MT=C64);
}

#scope_file
_ggesx_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_ggesx_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_ggesx_matrix_storage :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _ggesx_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_ggesx_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_ggesx_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_ggesx_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_ggesx_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_ggesx_lange :: (
    $$norm: u8, $$m: int, $$n: int,
    a: $A/ArrayType2D, $$lda: int, work: *$WORK/ArrayType1D,
    $RT: Type, $MT: Type
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

_ggesx_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int,
    a: $A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
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

_ggesx_laset :: (
    $$uplo: u8, $$m: int, $$n: int,
    alpha: $MT, beta: MT,
    a: *$A/ArrayType2D, $$lda: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT;
} {
    #if MT == float32 {
        slaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == float64 {
        dlaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == C32 {
        claset(uplo, m, n, alpha, beta, a, lda);
    } else {
        zlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_ggesx_lascl :: (
    $$type: u8, $$kl: int, $$ku: int,
    cfrom: $RT, cto: RT,
    $$m: int, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    info: *int, $RT: Type, $MT: Type
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

_ggesx_scale_real_segment :: (x: *$X/ArrayType1D, $$offset: int, $$count: int, factor: $RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    if count <= 0 then return;
    for i:offset..offset + count - 1 {
        fset(x, i, factor * fget(x.*, i));
    }
}

_ggesx_scale_complex_segment :: (x: *$X/ArrayType1D, $$offset: int, $$count: int, factor: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    if count <= 0 then return;
    factor_t := _ggesx_complex_from_real(factor, RT, MT);
    for i:offset..offset + count - 1 {
        fset(x, i, factor_t * fget(x.*, i));
    }
}

_ggesx_real :: (
    $$jobvsl: u8, $$jobvsr: u8, $$sort: u8,
    selctg: ($RT, RT, RT) -> bool,
    $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    sdim: *int,
    alphar: *$ALPHAR/ArrayType1D, alphai: *$ALPHAI/ArrayType1D, beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D, $$ldvsl: int,
    vsr: *$VSR/ArrayType2D, $$ldvsr: int,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int,
    bwork: *$BWORK/ArrayType1D, info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT &&
           ALPHAR.data_type == RT && ALPHAI.data_type == RT && BETA.data_type == RT &&
           VSL.data_type == RT && VSR.data_type == RT &&
           RCONDE.data_type == RT && RCONDV.data_type == RT &&
           WORK.data_type == RT && IWORK.data_type == int && BWORK.data_type == bool;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if RT == float32 { "SGGESX" } else { "DGGESX" };
    geqrf_name :: #run if RT == float32 { "SGEQRF" } else { "DGEQRF" };
    ormqr_name :: #run if RT == float32 { "SORMQR" } else { "DORMQR" };
    orgqr_name :: #run if RT == float32 { "SORGQR" } else { "DORGQR" };

    cursl, ilascl, ilbscl, ilvsl, ilvsr, lastsl, lquery, lst2sl, wantsb, wantse, wantsn, wantst, wantsv : bool = ---;
    i, icols, ierr, ijob, ijobvl, ijobvr, ihi, ileft, ilo, ip, iright, irows, itau, iwrk, liwmin, lwrk, maxwrk, minwrk : int = ---;
    anrm, anrmto, bignum, bnrm, bnrmto, eps, pl, pr, safmax, safmin, scale, smlnum : RT = ---;
    dif: [2] RT;

    info.* = 0;
    if lsame(jobvsl, #char "N") {
        ijobvl = 1; ilvsl = false;
    } else if lsame(jobvsl, #char "V") {
        ijobvl = 2; ilvsl = true;
    } else {
        ijobvl = -1; ilvsl = false;
    }

    if lsame(jobvsr, #char "N") {
        ijobvr = 1; ilvsr = false;
    } else if lsame(jobvsr, #char "V") {
        ijobvr = 2; ilvsr = true;
    } else {
        ijobvr = -1; ilvsr = false;
    }

    wantst = lsame(sort, #char "S");
    wantsn = lsame(sense, #char "N");
    wantse = lsame(sense, #char "E");
    wantsv = lsame(sense, #char "V");
    wantsb = lsame(sense, #char "B");
    lquery = lwork == -1 || liwork == -1;

    if wantsn {
        ijob = 0;
    } else if wantse {
        ijob = 1;
    } else if wantsv {
        ijob = 2;
    } else if wantsb {
        ijob = 4;
    } else {
        ijob = -1;
    }

    if ijobvl <= 0 {
        info.* = -1;
    } else if ijobvr <= 0 {
        info.* = -2;
    } else if !wantst && !lsame(sort, #char "N") {
        info.* = -3;
    } else if !(wantsn || wantse || wantsv || wantsb) || (!wantst && !wantsn) {
        info.* = -5;
    } else if n < 0 {
        info.* = -6;
    } else if lda < max(1, n) {
        info.* = -8;
    } else if ldb < max(1, n) {
        info.* = -10;
    } else if ldvsl < 1 || (ilvsl && ldvsl < n) {
        info.* = -16;
    } else if ldvsr < 1 || (ilvsr && ldvsr < n) {
        info.* = -18;
    }

    if info.* == 0 {
        if n > 0 {
            minwrk = max(8*n, 6*n + 16);
            maxwrk = minwrk - n + n*ilaenv(1, geqrf_name, " "s, n, 1, n, 0);
            maxwrk = max(maxwrk, minwrk - n + n*ilaenv(1, ormqr_name, " "s, n, 1, n, -1));
            if ilvsl {
                maxwrk = max(maxwrk, minwrk - n + n*ilaenv(1, orgqr_name, " "s, n, 1, n, -1));
            }
            lwrk = maxwrk;
            if ijob >= 1 {
                lwrk = max(lwrk, (n*n)/2);
            }
        } else {
            minwrk = 1;
            maxwrk = 1;
            lwrk = 1;
        }

        liwmin = ifx wantsn || n == 0 then 1 else n + 6;
        fset(work, 1, _ggesx_roundup_lwork(lwrk, RT));
        fset(iwork, 1, liwmin);
        if lwork < minwrk && !lquery {
            info.* = -22;
        } else if liwork < liwmin && !lquery {
            info.* = -24;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        sdim.* = 0;
        return;
    }

    eps = _ggesx_lamch_precision(RT);
    safmin = _ggesx_lamch_safe_minimum(RT);
    safmax = one_r / safmin;
    smlnum = sqrt(safmin) / eps;
    bignum = one_r / smlnum;

    anrm = _ggesx_lange(#char "M", n, n, a.*, lda, work, RT, RT);
    ilascl = false;
    if anrm > zero_r && anrm < smlnum {
        anrmto = smlnum;
        ilascl = true;
    } else if anrm > bignum {
        anrmto = bignum;
        ilascl = true;
    }
    if ilascl {
        _ggesx_lascl(#char "G", 0, 0, anrm, anrmto, n, n, a, lda, *ierr, RT, RT);
    }

    bnrm = _ggesx_lange(#char "M", n, n, b.*, ldb, work, RT, RT);
    ilbscl = false;
    if bnrm > zero_r && bnrm < smlnum {
        bnrmto = smlnum;
        ilbscl = true;
    } else if bnrm > bignum {
        bnrmto = bignum;
        ilbscl = true;
    }
    if ilbscl {
        _ggesx_lascl(#char "G", 0, 0, bnrm, bnrmto, n, n, b, ldb, *ierr, RT, RT);
    }

    ileft = 1;
    iright = n + 1;
    iwrk = iright + n;
    lscale := _ggesx_subvector(work.*, ileft, n);
    rscale := _ggesx_subvector(work.*, iright, n);
    bal_work := _ggesx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if RT == float32 {
        sggbal(#char "P", n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, bal_work, *ierr);
    } else {
        dggbal(#char "P", n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, bal_work, *ierr);
    }

    irows = ihi - ilo + 1;
    icols = n - ilo + 1;
    itau = iwrk;
    iwrk = itau + irows;
    tau := _ggesx_subvector(work.*, itau, max(0, irows));
    wrk := _ggesx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));

    b_sub_raw := _ggesx_matrix_storage(b.*, ldb, ilo, ilo, irows, icols);
    b_sub := array2d_view(b_sub_raw, ldb, row_major=false);
    #if RT == float32 {
        sgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr);
    } else {
        dgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr);
    }

    a_sub_raw := _ggesx_matrix_storage(a.*, lda, ilo, ilo, irows, icols);
    a_sub := array2d_view(a_sub_raw, lda, row_major=false);
    #if RT == float32 {
        sormqr(#char "L", #char "T", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr);
    } else {
        dormqr(#char "L", #char "T", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr);
    }

    if ilvsl {
        _ggesx_laset(#char "F", n, n, zero_r, one_r, vsl, ldvsl, RT);
        if irows > 1 {
            b_lower_raw := _ggesx_matrix_storage(b.*, ldb, ilo + 1, ilo, irows - 1, irows - 1);
            vsl_lower_raw := _ggesx_matrix_storage(vsl.*, ldvsl, ilo + 1, ilo, irows - 1, irows - 1);
            b_lower := array2d_view(b_lower_raw, ldb, row_major=false);
            vsl_lower := array2d_view(vsl_lower_raw, ldvsl, row_major=false);
            _ggesx_lacpy(#char "L", irows - 1, irows - 1, *b_lower, ldb, vsl_lower, ldvsl, RT);
        }

        vsl_sub_raw := _ggesx_matrix_storage(vsl.*, ldvsl, ilo, ilo, irows, irows);
        vsl_sub := array2d_view(vsl_sub_raw, ldvsl, row_major=false);
        #if RT == float32 {
            sorgqr(irows, irows, irows, *vsl_sub, ldvsl, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            dorgqr(irows, irows, irows, *vsl_sub, ldvsl, tau, wrk, lwork - iwrk + 1, *ierr);
        }
    }

    if ilvsr {
        _ggesx_laset(#char "F", n, n, zero_r, one_r, vsr, ldvsr, RT);
    }

    #if RT == float32 {
        sgghrd(jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr, *ierr);
    } else {
        dgghrd(jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr, *ierr);
    }

    sdim.* = 0;
    iwrk = itau;
    wrk = _ggesx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if RT == float32 {
        shgeqz(#char "S", jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, wrk, lwork - iwrk + 1, *ierr);
    } else {
        dhgeqz(#char "S", jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, wrk, lwork - iwrk + 1, *ierr);
    }
    if ierr != 0 {
        if ierr > 0 && ierr <= n {
            info.* = ierr;
        } else if ierr > n && ierr <= 2*n {
            info.* = ierr - n;
        } else {
            info.* = n + 1;
        }
        fset(work, 1, _ggesx_roundup_lwork(maxwrk, RT));
        fset(iwork, 1, liwmin);
        return;
    }

    if wantst {
        if ilascl {
            _ggesx_scale_real_segment(alphar, 1, n, anrm / anrmto, RT);
            _ggesx_scale_real_segment(alphai, 1, n, anrm / anrmto, RT);
        }
        if ilbscl {
            _ggesx_scale_real_segment(beta, 1, n, bnrm / bnrmto, RT);
        }

        for i:1..n {
            fset(bwork, i, selctg(fget(alphar.*, i), fget(alphai.*, i), fget(beta.*, i)));
        }

        #if RT == float32 {
            stgsen(ijob, ilvsl, ilvsr, bwork, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, sdim, *pl, *pr, dif, wrk, lwork - iwrk + 1, iwork, liwork, *ierr);
        } else {
            dtgsen(ijob, ilvsl, ilvsr, bwork, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, sdim, *pl, *pr, dif, wrk, lwork - iwrk + 1, iwork, liwork, *ierr);
        }

        if ijob >= 1 {
            maxwrk = max(maxwrk, 2*sdim.* * (n - sdim.*));
        }
        if ierr == -22 {
            info.* = -22;
        } else {
            if ijob == 1 || ijob == 4 {
                fset(rconde, 1, pl);
                fset(rconde, 2, pr);
            }
            if ijob == 2 || ijob == 4 {
                fset(rcondv, 1, dif[0]);
                fset(rcondv, 2, dif[1]);
            }
            if ierr == 1 {
                info.* = n + 3;
            }
        }
    }

    if ilvsl {
        #if RT == float32 {
            sggbak(#char "P", #char "L", n, ilo, ihi, lscale, rscale, n, vsl, ldvsl, *ierr);
        } else {
            dggbak(#char "P", #char "L", n, ilo, ihi, lscale, rscale, n, vsl, ldvsl, *ierr);
        }
    }
    if ilvsr {
        #if RT == float32 {
            sggbak(#char "P", #char "R", n, ilo, ihi, lscale, rscale, n, vsr, ldvsr, *ierr);
        } else {
            dggbak(#char "P", #char "R", n, ilo, ihi, lscale, rscale, n, vsr, ldvsr, *ierr);
        }
    }

    if ilascl {
        for i:1..n {
            if fget(alphai.*, i) != zero_r {
                if (fget(alphar.*, i) / safmax) > (anrmto / anrm) || (safmin / fget(alphar.*, i)) > (anrm / anrmto) {
                    scale = abs(fget(a.*, i, i) / fget(alphar.*, i));
                    fset(beta, i, fget(beta.*, i) * scale);
                    fset(alphar, i, fget(alphar.*, i) * scale);
                    fset(alphai, i, fget(alphai.*, i) * scale);
                } else if (fget(alphai.*, i) / safmax) > (anrmto / anrm) || (safmin / fget(alphai.*, i)) > (anrm / anrmto) {
                    scale = abs(fget(a.*, i, i + 1) / fget(alphai.*, i));
                    fset(beta, i, fget(beta.*, i) * scale);
                    fset(alphar, i, fget(alphar.*, i) * scale);
                    fset(alphai, i, fget(alphai.*, i) * scale);
                }
            }
        }
    }

    if ilbscl {
        for i:1..n {
            if fget(alphai.*, i) != zero_r {
                if (fget(beta.*, i) / safmax) > (bnrmto / bnrm) || (safmin / fget(beta.*, i)) > (bnrm / bnrmto) {
                    scale = abs(fget(b.*, i, i) / fget(beta.*, i));
                    fset(beta, i, fget(beta.*, i) * scale);
                    fset(alphar, i, fget(alphar.*, i) * scale);
                    fset(alphai, i, fget(alphai.*, i) * scale);
                }
            }
        }
    }

    if ilascl {
        _ggesx_lascl(#char "H", 0, 0, anrmto, anrm, n, n, a, lda, *ierr, RT, RT);
        _ggesx_scale_real_segment(alphar, 1, n, anrm / anrmto, RT);
        _ggesx_scale_real_segment(alphai, 1, n, anrm / anrmto, RT);
    }

    if ilbscl {
        _ggesx_lascl(#char "U", 0, 0, bnrmto, bnrm, n, n, b, ldb, *ierr, RT, RT);
        _ggesx_scale_real_segment(beta, 1, n, bnrm / bnrmto, RT);
    }

    if wantst {
        lastsl = true;
        lst2sl = true;
        sdim.* = 0;
        ip = 0;
        for i:1..n {
            cursl = selctg(fget(alphar.*, i), fget(alphai.*, i), fget(beta.*, i));
            if fget(alphai.*, i) == zero_r {
                if cursl then sdim.* += 1;
                ip = 0;
                if cursl && !lastsl {
                    info.* = n + 2;
                }
            } else {
                if ip == 1 {
                    cursl = cursl || lastsl;
                    lastsl = cursl;
                    if cursl then sdim.* += 2;
                    ip = -1;
                    if cursl && !lst2sl {
                        info.* = n + 2;
                    }
                } else {
                    ip = 1;
                }
            }
            lst2sl = lastsl;
            lastsl = cursl;
        }
    }

    fset(work, 1, _ggesx_roundup_lwork(maxwrk, RT));
    fset(iwork, 1, liwmin);
}

_ggesx_complex :: (
    $$jobvsl: u8, $$jobvsr: u8, $$sort: u8,
    selctg: ($MT, MT) -> bool,
    $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    sdim: *int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D, $$ldvsl: int,
    vsr: *$VSR/ArrayType2D, $$ldvsr: int,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, $$liwork: int,
    bwork: *$BWORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && B.data_type == MT &&
           ALPHA.data_type == MT && BETA.data_type == MT &&
           VSL.data_type == MT && VSR.data_type == MT &&
           RCONDE.data_type == RT && RCONDV.data_type == RT &&
           WORK.data_type == MT && RWORK.data_type == RT &&
           IWORK.data_type == int && BWORK.data_type == bool;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    zero_t : MT = _ggesx_complex_from_real(zero_r, RT, MT);
    one_t  : MT = _ggesx_complex_from_real(one_r, RT, MT);
    routine_name :: #run if MT == C32 { "CGGESX" } else { "ZGGESX" };
    geqrf_name :: #run if MT == C32 { "CGEQRF" } else { "ZGEQRF" };
    unmqr_name :: #run if MT == C32 { "CUNMQR" } else { "ZUNMQR" };
    ungqr_name :: #run if MT == C32 { "CUNGQR" } else { "ZUNGQR" };

    cursl, ilascl, ilbscl, ilvsl, ilvsr, lastsl, lquery, wantsb, wantse, wantsn, wantst, wantsv : bool = ---;
    i, icols, ierr, ijob, ijobvl, ijobvr, ihi, ileft, ilo, iright, irows, irwrk, itau, iwrk, liwmin, lwrk, maxwrk, minwrk : int = ---;
    anrm, anrmto, bignum, bnrm, bnrmto, eps, pl, pr, smlnum : RT = ---;
    dif: [2] RT;

    info.* = 0;
    if lsame(jobvsl, #char "N") {
        ijobvl = 1; ilvsl = false;
    } else if lsame(jobvsl, #char "V") {
        ijobvl = 2; ilvsl = true;
    } else {
        ijobvl = -1; ilvsl = false;
    }

    if lsame(jobvsr, #char "N") {
        ijobvr = 1; ilvsr = false;
    } else if lsame(jobvsr, #char "V") {
        ijobvr = 2; ilvsr = true;
    } else {
        ijobvr = -1; ilvsr = false;
    }

    wantst = lsame(sort, #char "S");
    wantsn = lsame(sense, #char "N");
    wantse = lsame(sense, #char "E");
    wantsv = lsame(sense, #char "V");
    wantsb = lsame(sense, #char "B");
    lquery = lwork == -1 || liwork == -1;

    if wantsn {
        ijob = 0;
    } else if wantse {
        ijob = 1;
    } else if wantsv {
        ijob = 2;
    } else if wantsb {
        ijob = 4;
    } else {
        ijob = -1;
    }

    if ijobvl <= 0 {
        info.* = -1;
    } else if ijobvr <= 0 {
        info.* = -2;
    } else if !wantst && !lsame(sort, #char "N") {
        info.* = -3;
    } else if !(wantsn || wantse || wantsv || wantsb) || (!wantst && !wantsn) {
        info.* = -5;
    } else if n < 0 {
        info.* = -6;
    } else if lda < max(1, n) {
        info.* = -8;
    } else if ldb < max(1, n) {
        info.* = -10;
    } else if ldvsl < 1 || (ilvsl && ldvsl < n) {
        info.* = -15;
    } else if ldvsr < 1 || (ilvsr && ldvsr < n) {
        info.* = -17;
    }

    if info.* == 0 {
        if n > 0 {
            minwrk = 2*n;
            maxwrk = n*(1 + ilaenv(1, geqrf_name, " "s, n, 1, n, 0));
            maxwrk = max(maxwrk, n*(1 + ilaenv(1, unmqr_name, " "s, n, 1, n, -1)));
            if ilvsl {
                maxwrk = max(maxwrk, n*(1 + ilaenv(1, ungqr_name, " "s, n, 1, n, -1)));
            }
            lwrk = maxwrk;
            if ijob >= 1 {
                lwrk = max(lwrk, (n*n)/2);
            }
        } else {
            minwrk = 1;
            maxwrk = 1;
            lwrk = 1;
        }

        liwmin = ifx wantsn || n == 0 then 1 else n + 2;
        if n == 0 {
            fset(work, 1, _ggesx_complex_from_real(one_r, RT, MT));
        } else {
            fset(work, 1, _ggesx_complex_from_real(_ggesx_roundup_lwork(lwrk, RT), RT, MT));
        }
        fset(iwork, 1, liwmin);
        if lwork < minwrk && !lquery {
            info.* = -21;
        } else if liwork < liwmin && !lquery {
            info.* = -24;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        sdim.* = 0;
        return;
    }

    eps = _ggesx_lamch_precision(RT);
    smlnum = _ggesx_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    anrm = _ggesx_lange(#char "M", n, n, a.*, lda, rwork, RT, MT);
    ilascl = false;
    if anrm > zero_r && anrm < smlnum {
        anrmto = smlnum;
        ilascl = true;
    } else if anrm > bignum {
        anrmto = bignum;
        ilascl = true;
    }
    if ilascl {
        _ggesx_lascl(#char "G", 0, 0, anrm, anrmto, n, n, a, lda, *ierr, RT, MT);
    }

    bnrm = _ggesx_lange(#char "M", n, n, b.*, ldb, rwork, RT, MT);
    ilbscl = false;
    if bnrm > zero_r && bnrm < smlnum {
        bnrmto = smlnum;
        ilbscl = true;
    } else if bnrm > bignum {
        bnrmto = bignum;
        ilbscl = true;
    }
    if ilbscl {
        _ggesx_lascl(#char "G", 0, 0, bnrm, bnrmto, n, n, b, ldb, *ierr, RT, MT);
    }

    ileft = 1;
    iright = n + 1;
    irwrk = iright + n;
    lscale := _ggesx_subvector(rwork.*, ileft, n);
    rscale := _ggesx_subvector(rwork.*, iright, n);
    bal_work := _ggesx_subvector(rwork.*, irwrk, max(0, 8*n - irwrk + 1));
    #if MT == C32 {
        cggbal(#char "P", n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, bal_work, *ierr);
    } else {
        zggbal(#char "P", n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, bal_work, *ierr);
    }

    irows = ihi - ilo + 1;
    icols = n - ilo + 1;
    itau = 1;
    iwrk = itau + irows;
    tau := _ggesx_subvector(work.*, itau, max(0, irows));
    wrk := _ggesx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));

    b_sub_raw := _ggesx_matrix_storage(b.*, ldb, ilo, ilo, irows, icols);
    b_sub := array2d_view(b_sub_raw, ldb, row_major=false);
    #if MT == C32 {
        cgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr);
    } else {
        zgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr);
    }

    a_sub_raw := _ggesx_matrix_storage(a.*, lda, ilo, ilo, irows, icols);
    a_sub := array2d_view(a_sub_raw, lda, row_major=false);
    #if MT == C32 {
        cunmqr(#char "L", #char "C", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr);
    } else {
        zunmqr(#char "L", #char "C", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr);
    }

    if ilvsl {
        _ggesx_laset(#char "F", n, n, zero_t, one_t, vsl, ldvsl, MT);
        if irows > 1 {
            b_lower_raw := _ggesx_matrix_storage(b.*, ldb, ilo + 1, ilo, irows - 1, irows - 1);
            vsl_lower_raw := _ggesx_matrix_storage(vsl.*, ldvsl, ilo + 1, ilo, irows - 1, irows - 1);
            b_lower := array2d_view(b_lower_raw, ldb, row_major=false);
            vsl_lower := array2d_view(vsl_lower_raw, ldvsl, row_major=false);
            _ggesx_lacpy(#char "L", irows - 1, irows - 1, *b_lower, ldb, vsl_lower, ldvsl, MT);
        }

        vsl_sub_raw := _ggesx_matrix_storage(vsl.*, ldvsl, ilo, ilo, irows, irows);
        vsl_sub := array2d_view(vsl_sub_raw, ldvsl, row_major=false);
        #if MT == C32 {
            cungqr(irows, irows, irows, *vsl_sub, ldvsl, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            zungqr(irows, irows, irows, *vsl_sub, ldvsl, tau, wrk, lwork - iwrk + 1, *ierr);
        }
    }

    if ilvsr {
        _ggesx_laset(#char "F", n, n, zero_t, one_t, vsr, ldvsr, MT);
    }

    #if MT == C32 {
        cgghrd(jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr, *ierr);
    } else {
        zgghrd(jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr, *ierr);
    }

    sdim.* = 0;
    iwrk = itau;
    wrk = _ggesx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    qz_rwork := _ggesx_subvector(rwork.*, irwrk, max(0, n));
    #if MT == C32 {
        chgeqz(#char "S", jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, wrk, lwork - iwrk + 1, qz_rwork, *ierr);
    } else {
        zhgeqz(#char "S", jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, wrk, lwork - iwrk + 1, qz_rwork, *ierr);
    }
    if ierr != 0 {
        if ierr > 0 && ierr <= n {
            info.* = ierr;
        } else if ierr > n && ierr <= 2*n {
            info.* = ierr - n;
        } else {
            info.* = n + 1;
        }
        fset(work, 1, _ggesx_complex_from_real(_ggesx_roundup_lwork(maxwrk, RT), RT, MT));
        fset(iwork, 1, liwmin);
        return;
    }

    if wantst {
        if ilascl {
            _ggesx_scale_complex_segment(alpha, 1, n, anrm / anrmto, RT, MT);
        }
        if ilbscl {
            _ggesx_scale_complex_segment(beta, 1, n, bnrm / bnrmto, RT, MT);
        }

        for i:1..n {
            fset(bwork, i, selctg(fget(alpha.*, i), fget(beta.*, i)));
        }

        #if MT == C32 {
            ctgsen(ijob, ilvsl, ilvsr, bwork, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, sdim, *pl, *pr, dif, wrk, lwork - iwrk + 1, iwork, liwork, *ierr);
        } else {
            ztgsen(ijob, ilvsl, ilvsr, bwork, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, sdim, *pl, *pr, dif, wrk, lwork - iwrk + 1, iwork, liwork, *ierr);
        }

        if ijob >= 1 {
            maxwrk = max(maxwrk, 2*sdim.* * (n - sdim.*));
        }
        if ierr == -21 {
            info.* = -21;
        } else {
            if ijob == 1 || ijob == 4 {
                fset(rconde, 1, pl);
                fset(rconde, 2, pr);
            }
            if ijob == 2 || ijob == 4 {
                fset(rcondv, 1, dif[0]);
                fset(rcondv, 2, dif[1]);
            }
            if ierr == 1 {
                info.* = n + 3;
            }
        }
    }

    if ilvsl {
        #if MT == C32 {
            cggbak(#char "P", #char "L", n, ilo, ihi, lscale, rscale, n, vsl, ldvsl, *ierr);
        } else {
            zggbak(#char "P", #char "L", n, ilo, ihi, lscale, rscale, n, vsl, ldvsl, *ierr);
        }
    }
    if ilvsr {
        #if MT == C32 {
            cggbak(#char "P", #char "R", n, ilo, ihi, lscale, rscale, n, vsr, ldvsr, *ierr);
        } else {
            zggbak(#char "P", #char "R", n, ilo, ihi, lscale, rscale, n, vsr, ldvsr, *ierr);
        }
    }

    if ilascl {
        _ggesx_lascl(#char "U", 0, 0, anrmto, anrm, n, n, a, lda, *ierr, RT, MT);
        _ggesx_scale_complex_segment(alpha, 1, n, anrm / anrmto, RT, MT);
    }
    if ilbscl {
        _ggesx_lascl(#char "U", 0, 0, bnrmto, bnrm, n, n, b, ldb, *ierr, RT, MT);
        _ggesx_scale_complex_segment(beta, 1, n, bnrm / bnrmto, RT, MT);
    }

    if wantst {
        lastsl = true;
        sdim.* = 0;
        for i:1..n {
            cursl = selctg(fget(alpha.*, i), fget(beta.*, i));
            if cursl then sdim.* += 1;
            if cursl && !lastsl {
                info.* = n + 2;
            }
            lastsl = cursl;
        }
    }

    fset(work, 1, _ggesx_complex_from_real(_ggesx_roundup_lwork(maxwrk, RT), RT, MT));
    fset(iwork, 1, liwmin);
}
