// ? LAPACK *gees()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgees/cgees/dgees/zgees
// *GEES computes the eigenvalues, Schur form, and optionally the
// Schur vectors for a general matrix, with optional sorting.

sgees :: (
    $$jobvs: u8,
    $$sort: u8,
    select: (float32, float32) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    sdim: *int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    vs: *$VS/ArrayType2D,
    $$ldvs: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _gees_real(jobvs, sort, select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info, RT=float32);
}

dgees :: (
    $$jobvs: u8,
    $$sort: u8,
    select: (float64, float64) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    sdim: *int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    vs: *$VS/ArrayType2D,
    $$ldvs: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _gees_real(jobvs, sort, select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info, RT=float64);
}

cgees :: (
    $$jobvs: u8,
    $$sort: u8,
    select: (C32) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    sdim: *int,
    w: *$W/ArrayType1D,
    vs: *$VS/ArrayType2D,
    $$ldvs: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _gees_complex(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info, RT=float32, MT=C32);
}

zgees :: (
    $$jobvs: u8,
    $$sort: u8,
    select: (C64) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    sdim: *int,
    w: *$W/ArrayType1D,
    vs: *$VS/ArrayType2D,
    $$ldvs: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _gees_complex(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info, RT=float64, MT=C64);
}

#scope_file
_gees_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gees_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gees_matrix_column :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    return _gees_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gees_matrix_row_suffix :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*lda + (row - 1);
    view.count = if count <= 0 then 0 else 1 + (count - 1)*lda;
    return view;
}

_gees_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_gees_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_gees_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gees_lange :: (
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

_gees_lacpy :: (
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

_gees_lascl :: (
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

_gees_swap :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sswap(n, x, incx, y, incy);
    } else if MT == float64 {
        dswap(n, x, incx, y, incy);
    } else if MT == C32 {
        cswap(n, x, incx, y, incy);
    } else {
        zswap(n, x, incx, y, incy);
    }
}

_gees_copy_diag :: (a: $A/ArrayType2D, $$lda: int, $$n: int, d: *$D/ArrayType1D) #modify {
    return A.data_type == D.data_type;
} {
    for i:1..n {
        fset(d, i, fget(a, i, i));
    }
}

_gees_swap_real_schur_pair :: (
    $$n: int,
    $$i: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    vs: *$VS/ArrayType2D,
    $$ldvs: int,
    $$wantvs: bool,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           VS.data_type == RT;
} {
    if i > 1 {
        col1 := _gees_matrix_column(a.*, lda, 1, i, i - 1);
        col2 := _gees_matrix_column(a.*, lda, 1, i + 1, i - 1);
        _gees_swap(i - 1, *col1, 1, *col2, 1, RT);
    }

    if n > i + 1 {
        row1 := _gees_matrix_row_suffix(a.*, lda, i, i + 2, n - i - 1);
        row2 := _gees_matrix_row_suffix(a.*, lda, i + 1, i + 2, n - i - 1);
        _gees_swap(n - i - 1, *row1, lda, *row2, lda, RT);
    }

    if wantvs {
        col1 := _gees_matrix_column(vs.*, ldvs, 1, i, n);
        col2 := _gees_matrix_column(vs.*, ldvs, 1, i + 1, n);
        _gees_swap(n, *col1, 1, *col2, 1, RT);
    }
}

_gees_trsen_real :: (
    $$job: u8,
    $$compq: u8,
    select: $SELECT/ArrayType1D,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    vs: *$VS/ArrayType2D,
    $$ldvs: int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    m: *int,
    s: *$RT,
    sep: *RT,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    iwork: *$IWORK/ArrayType1D,
    $$liwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           VS.data_type == RT &&
           WR.data_type == RT &&
           WI.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int &&
           SELECT.data_type == bool;
} {
    #if RT == float32 {
        strsen(job, compq, select, n, a, lda, vs, ldvs, wr, wi, m, s, sep, work, lwork, iwork, liwork, info);
    } else {
        dtrsen(job, compq, select, n, a, lda, vs, ldvs, wr, wi, m, s, sep, work, lwork, iwork, liwork, info);
    }
}

_gees_trsen_complex :: (
    $$job: u8,
    $$compq: u8,
    select: $SELECT/ArrayType1D,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    vs: *$VS/ArrayType2D,
    $$ldvs: int,
    w: *$W/ArrayType1D,
    m: *int,
    s: *$RT,
    sep: *RT,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           VS.data_type == MT &&
           W.data_type == MT &&
           WORK.data_type == MT &&
           SELECT.data_type == bool;
} {
    #if MT == C32 {
        ctrsen(job, compq, select, n, a, lda, vs, ldvs, w, m, s, sep, work, lwork, info);
    } else {
        ztrsen(job, compq, select, n, a, lda, vs, ldvs, w, m, s, sep, work, lwork, info);
    }
}

_gees_real :: (
    $$jobvs: u8,
    $$sort: u8,
    select: ($RT, RT) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    sdim: *int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    vs: *$VS/ArrayType2D,
    $$ldvs: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    bwork: *$BWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           WR.data_type == RT &&
           WI.data_type == RT &&
           VS.data_type == RT &&
           WORK.data_type == RT &&
           BWORK.data_type == bool;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if RT == float32 { "SGEES " } else { "DGEES " };
    gehrd_name :: #run if RT == float32 { "SGEHRD" } else { "DGEHRD" };
    orghr_name :: #run if RT == float32 { "SORGHR" } else { "DORGHR" };

    cursl, lastsl, lquery, lst2sl, scalea, wantst, wantvs : bool = ---;
    hswork, i, i1, i2, ibal, icond, ierr, ieval, ihi, ilo, inxt, ip, itau, iwrk, maxwrk, minwrk : int = ---;
    anrm, bignum, cscale, eps, s, sep, smlnum : RT = ---;
    dum : [1] RT;
    idum : [1] int;
    hseqr_query : [1] RT;

    info.* = 0;
    lquery = lwork == -1;
    wantvs = lsame(jobvs, #char "V");
    wantst = lsame(sort, #char "S");

    if !wantvs && !lsame(jobvs, #char "N") {
        info.* = -1;
    } else if !wantst && !lsame(sort, #char "N") {
        info.* = -2;
    } else if n < 0 {
        info.* = -4;
    } else if lda < max(1, n) {
        info.* = -6;
    } else if ldvs < 1 || (wantvs && ldvs < n) {
        info.* = -11;
    }

    if info.* == 0 {
        if n == 0 {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = 2*n + n*ilaenv(1, gehrd_name, " "s, n, 1, n, 0);
            minwrk = 3*n;

            #if RT == float32 {
                shseqr(#char "S", jobvs, n, 1, n, a, lda, wr, wi, vs, ldvs, hseqr_query, -1, *ieval);
            } else {
                dhseqr(#char "S", jobvs, n, 1, n, a, lda, wr, wi, vs, ldvs, hseqr_query, -1, *ieval);
            }
            hswork = cast(int) hseqr_query[0];

            if !wantvs {
                maxwrk = max(maxwrk, n + hswork);
            } else {
                maxwrk = max(maxwrk, 2*n + (n - 1)*ilaenv(1, orghr_name, " "s, n, 1, n, -1));
                maxwrk = max(maxwrk, n + hswork);
            }
        }

        fset(work, 1, ncast(RT, maxwrk));
        if lwork < minwrk && !lquery {
            info.* = -13;
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

    eps = _gees_lamch_precision(RT);
    smlnum = _gees_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    anrm = _gees_lange(#char "M", n, n, a.*, lda, *dum, RT, RT);
    scalea = false;
    if anrm > zero_r && anrm < smlnum {
        scalea = true;
        cscale = smlnum;
    } else if anrm > bignum {
        scalea = true;
        cscale = bignum;
    }
    if scalea {
        _gees_lascl(#char "G", 0, 0, anrm, cscale, n, n, a, lda, *ierr, RT, RT);
    }

    ibal = 1;
    bal_scale := _gees_subvector(work.*, ibal, n);
    #if RT == float32 {
        sgebal(#char "P", n, a, lda, *ilo, *ihi, bal_scale, *ierr);
    } else {
        dgebal(#char "P", n, a, lda, *ilo, *ihi, bal_scale, *ierr);
    }

    itau = n + ibal;
    iwrk = n + itau;
    tau := _gees_subvector(work.*, itau, n);
    wrk := _gees_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if RT == float32 {
        sgehrd(n, ilo, ihi, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    } else {
        dgehrd(n, ilo, ihi, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    }

    if wantvs {
        _gees_lacpy(#char "L", n, n, a.*, lda, vs, ldvs, RT);
        #if RT == float32 {
            sorghr(n, ilo, ihi, vs, ldvs, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            dorghr(n, ilo, ihi, vs, ldvs, tau, wrk, lwork - iwrk + 1, *ierr);
        }
    }

    sdim.* = 0;
    iwrk = itau;
    wrk = _gees_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if RT == float32 {
        shseqr(#char "S", jobvs, n, ilo, ihi, a, lda, wr, wi, vs, ldvs, wrk, lwork - iwrk + 1, *ieval);
    } else {
        dhseqr(#char "S", jobvs, n, ilo, ihi, a, lda, wr, wi, vs, ldvs, wrk, lwork - iwrk + 1, *ieval);
    }
    if ieval > 0 {
        info.* = ieval;
    }

    if wantst && info.* == 0 {
        if scalea {
            _gees_scale_real_segment(wr, 1, n, anrm / cscale, RT);
            _gees_scale_real_segment(wi, 1, n, anrm / cscale, RT);
        }

        for i:1..n {
            fset(bwork, i, select(fget(wr.*, i), fget(wi.*, i)));
        }

        wrk = _gees_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        _gees_trsen_real(#char "N", jobvs, bwork, n, a, lda, vs, ldvs, wr, wi, sdim, *s, *sep, wrk, lwork - iwrk + 1, *idum, 1, *icond, RT);
        if icond > 0 {
            info.* = n + icond;
        }
    }

    if wantvs {
        #if RT == float32 {
            sgebak(#char "P", #char "R", n, ilo, ihi, bal_scale, n, vs, ldvs, *ierr);
        } else {
            dgebak(#char "P", #char "R", n, ilo, ihi, bal_scale, n, vs, ldvs, *ierr);
        }
    }

    if scalea {
        _gees_lascl(#char "H", 0, 0, cscale, anrm, n, n, a, lda, *ierr, RT, RT);
        _gees_copy_diag(a.*, lda, n, wr);

        if cscale == smlnum {
            if ieval > 0 {
                i1 = ieval + 1;
                i2 = ihi - 1;
                _gees_scale_real_segment(wi, 1, ilo - 1, anrm / cscale, RT);
            } else if wantst {
                i1 = 1;
                i2 = n - 1;
            } else {
                i1 = ilo;
                i2 = ihi - 1;
            }

            inxt = i1 - 1;
            for i:i1..i2 {
                if i < inxt continue;
                if fget(wi.*, i) == zero_r {
                    inxt = i + 1;
                } else {
                    if fget(a.*, i + 1, i) == zero_r {
                        fset(wi, i, zero_r);
                        fset(wi, i + 1, zero_r);
                    } else if fget(a.*, i + 1, i) != zero_r && fget(a.*, i, i + 1) == zero_r {
                        fset(wi, i, zero_r);
                        fset(wi, i + 1, zero_r);
                        _gees_swap_real_schur_pair(n, i, a, lda, vs, ldvs, wantvs, RT);
                        fset(a, i, i + 1, fget(a.*, i + 1, i));
                        fset(a, i + 1, i, zero_r);
                    }
                    inxt = i + 2;
                }
            }
        }

        _gees_scale_real_segment(wi, ieval + 1, n - ieval, anrm / cscale, RT);
    }

    if wantst && info.* == 0 {
        lastsl = true;
        lst2sl = true;
        sdim.* = 0;
        ip = 0;
        for i:1..n {
            cursl = select(fget(wr.*, i), fget(wi.*, i));
            if fget(wi.*, i) == zero_r {
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

    fset(work, 1, ncast(RT, maxwrk));
}

_gees_complex :: (
    $$jobvs: u8,
    $$sort: u8,
    select: ($MT) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    sdim: *int,
    w: *$W/ArrayType1D,
    vs: *$VS/ArrayType2D,
    $$ldvs: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    bwork: *$BWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           W.data_type == MT &&
           VS.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           BWORK.data_type == bool;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == C32 { "CGEES " } else { "ZGEES " };
    gehrd_name :: #run if MT == C32 { "CGEHRD" } else { "ZGEHRD" };
    unghr_name :: #run if MT == C32 { "CUNGHR" } else { "ZUNGHR" };

    lquery, scalea, wantst, wantvs : bool = ---;
    hswork, i, ibal, icond, ierr, ieval, ihi, ilo, itau, iwrk, maxwrk, minwrk : int = ---;
    anrm, bignum, cscale, eps, s, sep, smlnum : RT = ---;
    dum : [1] RT;
    hseqr_query : [1] MT;

    info.* = 0;
    lquery = lwork == -1;
    wantvs = lsame(jobvs, #char "V");
    wantst = lsame(sort, #char "S");

    if !wantvs && !lsame(jobvs, #char "N") {
        info.* = -1;
    } else if !wantst && !lsame(sort, #char "N") {
        info.* = -2;
    } else if n < 0 {
        info.* = -4;
    } else if lda < max(1, n) {
        info.* = -6;
    } else if ldvs < 1 || (wantvs && ldvs < n) {
        info.* = -10;
    }

    if info.* == 0 {
        if n == 0 {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = n + n*ilaenv(1, gehrd_name, " "s, n, 1, n, 0);
            minwrk = 2*n;

            #if MT == C32 {
                chseqr(#char "S", jobvs, n, 1, n, a, lda, w, vs, ldvs, hseqr_query, -1, *ieval);
            } else {
                zhseqr(#char "S", jobvs, n, 1, n, a, lda, w, vs, ldvs, hseqr_query, -1, *ieval);
            }
            hswork = cast(int) real(hseqr_query[0]);

            if !wantvs {
                maxwrk = max(maxwrk, hswork);
            } else {
                maxwrk = max(maxwrk, n + (n - 1)*ilaenv(1, unghr_name, " "s, n, 1, n, -1));
                maxwrk = max(maxwrk, hswork);
            }
        }

        fset(work, 1, _gees_complex_from_real(ncast(RT, maxwrk), RT, MT));
        if lwork < minwrk && !lquery {
            info.* = -12;
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

    eps = _gees_lamch_precision(RT);
    smlnum = _gees_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    anrm = _gees_lange(#char "M", n, n, a.*, lda, *dum, RT, MT);
    scalea = false;
    if anrm > zero_r && anrm < smlnum {
        scalea = true;
        cscale = smlnum;
    } else if anrm > bignum {
        scalea = true;
        cscale = bignum;
    }
    if scalea {
        _gees_lascl(#char "G", 0, 0, anrm, cscale, n, n, a, lda, *ierr, RT, MT);
    }

    ibal = 1;
    bal_scale := _gees_subvector(rwork.*, ibal, n);
    #if MT == C32 {
        cgebal(#char "P", n, a, lda, *ilo, *ihi, bal_scale, *ierr);
    } else {
        zgebal(#char "P", n, a, lda, *ilo, *ihi, bal_scale, *ierr);
    }

    itau = 1;
    iwrk = n + itau;
    tau := _gees_subvector(work.*, itau, n);
    wrk := _gees_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if MT == C32 {
        cgehrd(n, ilo, ihi, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    } else {
        zgehrd(n, ilo, ihi, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    }

    if wantvs {
        _gees_lacpy(#char "L", n, n, a.*, lda, vs, ldvs, MT);
        #if MT == C32 {
            cunghr(n, ilo, ihi, vs, ldvs, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            zunghr(n, ilo, ihi, vs, ldvs, tau, wrk, lwork - iwrk + 1, *ierr);
        }
    }

    sdim.* = 0;
    iwrk = itau;
    wrk = _gees_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if MT == C32 {
        chseqr(#char "S", jobvs, n, ilo, ihi, a, lda, w, vs, ldvs, wrk, lwork - iwrk + 1, *ieval);
    } else {
        zhseqr(#char "S", jobvs, n, ilo, ihi, a, lda, w, vs, ldvs, wrk, lwork - iwrk + 1, *ieval);
    }
    if ieval > 0 {
        info.* = ieval;
    }

    if wantst && info.* == 0 {
        if scalea {
            _gees_scale_complex_segment(w, 1, n, anrm / cscale, RT, MT);
        }

        for i:1..n {
            fset(bwork, i, select(fget(w.*, i)));
        }

        wrk = _gees_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        _gees_trsen_complex(#char "N", jobvs, bwork, n, a, lda, vs, ldvs, w, sdim, *s, *sep, wrk, lwork - iwrk + 1, *icond, RT, MT);
    }

    if wantvs {
        #if MT == C32 {
            cgebak(#char "P", #char "R", n, ilo, ihi, bal_scale, n, vs, ldvs, *ierr);
        } else {
            zgebak(#char "P", #char "R", n, ilo, ihi, bal_scale, n, vs, ldvs, *ierr);
        }
    }

    if scalea {
        _gees_lascl(#char "U", 0, 0, cscale, anrm, n, n, a, lda, *ierr, RT, MT);
        _gees_copy_diag(a.*, lda, n, w);
    }

    fset(work, 1, _gees_complex_from_real(ncast(RT, maxwrk), RT, MT));
}

_gees_scale_real_segment :: (
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

_gees_scale_complex_segment :: (
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
