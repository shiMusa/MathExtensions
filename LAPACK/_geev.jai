// ? LAPACK *geev()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeev/cgeev/dgeev/zgeev
// *GEEV computes the eigenvalues and, optionally, the left and/or
// right eigenvectors for a general matrix.

sgeev :: #bake_arguments _geev_real(RT=float32);
dgeev :: #bake_arguments _geev_real(RT=float64);
cgeev :: #bake_arguments _geev_complex(RT=float32, MT=C32);
zgeev :: #bake_arguments _geev_complex(RT=float64, MT=C64);

#scope_file
_geev_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geev_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_geev_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _geev_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geev_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_geev_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_geev_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_geev_lange :: (
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

_geev_lacpy :: (
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

_geev_lascl :: (
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

_geev_nrm2 :: (
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

_geev_scal_real :: (
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

_geev_scale_complex :: (
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

_geev_iamax_1based :: (
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

_geev_rotg :: (
    a: $RT,
    b: RT,
    c: *RT,
    s: *RT,
    r: *RT,
    $RT: Type
) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        slartg(a, b, c, s, r);
    } else {
        dlartg(a, b, c, s, r);
    }
}

_geev_rot :: (
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
           X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 {
        srot(n, x, incx, y, incy, c, s);
    } else {
        drot(n, x, incx, y, incy, c, s);
    }
}

_geev_lapy2 :: (x: $RT, y: RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slapy2(x, y);
    } else {
        return dlapy2(x, y);
    }
}

_geev_scale_real_segment :: (
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

_geev_scale_complex_segment :: (
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
    scale_t := cpx(factor, zero(RT));
    for i:0..count-1 {
        fset(x, offset + i, scale_t * fget(x.*, offset + i));
    }
}

_geev_real :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           WR.data_type == RT &&
           WI.data_type == RT &&
           VL.data_type == RT &&
           VR.data_type == RT &&
           WORK.data_type == RT;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if RT == float32 { "SGEEV " } else { "DGEEV " };
    gehrd_name :: #run if RT == float32 { "SGEHRD" } else { "DGEHRD" };
    orghr_name :: #run if RT == float32 { "SORGHR" } else { "DORGHR" };

    wantvl, wantvr, lquery, scalea : bool = ---;
    side : u8 = ---;
    hswork, i, ibal, ierr, ihi, ilo, itau, iwrk, k, lwork_trevc, maxwrk, minwrk, nout : int = ---;
    anrm, bignum, cs, cscale, eps, r, scl, smlnum, sn : RT = ---;
    select_dummy : [1] bool;
    dum : [1] RT;
    hseqr_query : [1] RT;
    trevc_query : [1] RT;

    info.* = 0;
    lquery = lwork == -1;
    wantvl = lsame(jobvl, #char "V");
    wantvr = lsame(jobvr, #char "V");

    if !wantvl && !lsame(jobvl, #char "N") {
        info.* = -1;
    } else if !wantvr && !lsame(jobvr, #char "N") {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ldvl < 1 || (wantvl && ldvl < n) {
        info.* = -9;
    } else if ldvr < 1 || (wantvr && ldvr < n) {
        info.* = -11;
    }

    if info.* == 0 {
        if n == 0 {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = 2*n + n*ilaenv(1, gehrd_name, " "s, n, 1, n, 0);
            if wantvl {
                minwrk = 4*n;
                maxwrk = max(maxwrk, 2*n + (n - 1)*ilaenv(1, orghr_name, " "s, n, 1, n, -1));
                #if RT == float32 {
                    strevc3(#char "L", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, *ierr);
                    shseqr(#char "S", #char "V", n, 1, n, a, lda, wr, wi, vl, ldvl, hseqr_query, -1, *info);
                } else {
                    dtrevc3(#char "L", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, *ierr);
                    dhseqr(#char "S", #char "V", n, 1, n, a, lda, wr, wi, vl, ldvl, hseqr_query, -1, *info);
                }
                lwork_trevc = cast(int) trevc_query[0];
                hswork = cast(int) hseqr_query[0];
                maxwrk = max(maxwrk, n + hswork, n + lwork_trevc, 4*n);
            } else if wantvr {
                minwrk = 4*n;
                maxwrk = max(maxwrk, 2*n + (n - 1)*ilaenv(1, orghr_name, " "s, n, 1, n, -1));
                #if RT == float32 {
                    strevc3(#char "R", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, *ierr);
                    shseqr(#char "S", #char "V", n, 1, n, a, lda, wr, wi, vr, ldvr, hseqr_query, -1, *info);
                } else {
                    dtrevc3(#char "R", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, *ierr);
                    dhseqr(#char "S", #char "V", n, 1, n, a, lda, wr, wi, vr, ldvr, hseqr_query, -1, *info);
                }
                lwork_trevc = cast(int) trevc_query[0];
                hswork = cast(int) hseqr_query[0];
                maxwrk = max(maxwrk, n + hswork, n + lwork_trevc, 4*n);
            } else {
                minwrk = 3*n;
                #if RT == float32 {
                    shseqr(#char "E", #char "N", n, 1, n, a, lda, wr, wi, vr, ldvr, hseqr_query, -1, *info);
                } else {
                    dhseqr(#char "E", #char "N", n, 1, n, a, lda, wr, wi, vr, ldvr, hseqr_query, -1, *info);
                }
                hswork = cast(int) hseqr_query[0];
                maxwrk = max(maxwrk, n + hswork);
            }
            maxwrk = max(maxwrk, minwrk);
        }
        fset(work, 1, ncast(RT, maxwrk));
        if (lwork < minwrk) && !lquery {
            info.* = -13;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 then return;

    eps = _geev_lamch_precision(RT);
    smlnum = _geev_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    anrm = _geev_lange(#char "M", n, n, a.*, lda, *dum, RT, RT);
    scalea = false;
    if anrm > zero_r && anrm < smlnum {
        scalea = true;
        cscale = smlnum;
    } else if anrm > bignum {
        scalea = true;
        cscale = bignum;
    }
    if scalea {
        _geev_lascl(#char "G", 0, 0, anrm, cscale, n, n, a, lda, *ierr, RT, RT);
    }

    ibal = 1;
    #if RT == float32 {
        sgebal(#char "B", n, a, lda, *ilo, *ihi, _geev_subvector(work.*, ibal, n), *ierr);
    } else {
        dgebal(#char "B", n, a, lda, *ilo, *ihi, _geev_subvector(work.*, ibal, n), *ierr);
    }

    itau = ibal + n;
    iwrk = itau + n;
    tau := _geev_subvector(work.*, itau, n);
    wrk := _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if RT == float32 {
        sgehrd(n, ilo, ihi, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    } else {
        dgehrd(n, ilo, ihi, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    }

    if wantvl {
        side = #char "L";
        _geev_lacpy(#char "L", n, n, a.*, lda, vl, ldvl, RT);
        #if RT == float32 {
            sorghr(n, ilo, ihi, vl, ldvl, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            dorghr(n, ilo, ihi, vl, ldvl, tau, wrk, lwork - iwrk + 1, *ierr);
        }
        iwrk = itau;
        wrk = _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if RT == float32 {
            shseqr(#char "S", #char "V", n, ilo, ihi, a, lda, wr, wi, vl, ldvl, wrk, lwork - iwrk + 1, info);
        } else {
            dhseqr(#char "S", #char "V", n, ilo, ihi, a, lda, wr, wi, vl, ldvl, wrk, lwork - iwrk + 1, info);
        }
        if wantvr {
            side = #char "B";
            _geev_lacpy(#char "F", n, n, vl.*, ldvl, vr, ldvr, RT);
        }
    } else if wantvr {
        side = #char "R";
        _geev_lacpy(#char "L", n, n, a.*, lda, vr, ldvr, RT);
        #if RT == float32 {
            sorghr(n, ilo, ihi, vr, ldvr, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            dorghr(n, ilo, ihi, vr, ldvr, tau, wrk, lwork - iwrk + 1, *ierr);
        }
        iwrk = itau;
        wrk = _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if RT == float32 {
            shseqr(#char "S", #char "V", n, ilo, ihi, a, lda, wr, wi, vr, ldvr, wrk, lwork - iwrk + 1, info);
        } else {
            dhseqr(#char "S", #char "V", n, ilo, ihi, a, lda, wr, wi, vr, ldvr, wrk, lwork - iwrk + 1, info);
        }
    } else {
        iwrk = itau;
        wrk = _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if RT == float32 {
            shseqr(#char "E", #char "N", n, ilo, ihi, a, lda, wr, wi, vr, ldvr, wrk, lwork - iwrk + 1, info);
        } else {
            dhseqr(#char "E", #char "N", n, ilo, ihi, a, lda, wr, wi, vr, ldvr, wrk, lwork - iwrk + 1, info);
        }
    }

    if info.* != 0 goto end_scaling;

    if wantvl || wantvr {
        wrk = _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if RT == float32 {
            strevc3(side, #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, wrk, lwork - iwrk + 1, *ierr);
        } else {
            dtrevc3(side, #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, wrk, lwork - iwrk + 1, *ierr);
        }
    }

    if wantvl {
        bal_scale := _geev_subvector(work.*, ibal, n);
        #if RT == float32 {
            sgebak(#char "B", #char "L", n, ilo, ihi, bal_scale, n, vl, ldvl, *ierr);
        } else {
            dgebak(#char "B", #char "L", n, ilo, ihi, bal_scale, n, vl, ldvl, *ierr);
        }
        i = 1;
        while i <= n {
            if fget(wi.*, i) == zero_r {
                col := _geev_matrix_column(vl.*, ldvl, 1, i, n);
                scl = one_r / _geev_nrm2(n, col, 1, RT, RT);
                _geev_scal_real(n, scl, *col, 1, RT, RT);
                i += 1;
            } else if fget(wi.*, i) > zero_r {
                col1 := _geev_matrix_column(vl.*, ldvl, 1, i, n);
                col2 := _geev_matrix_column(vl.*, ldvl, 1, i + 1, n);
                scl = one_r / _geev_lapy2(_geev_nrm2(n, col1, 1, RT, RT), _geev_nrm2(n, col2, 1, RT, RT), RT);
                _geev_scal_real(n, scl, *col1, 1, RT, RT);
                _geev_scal_real(n, scl, *col2, 1, RT, RT);
                for k:1..n {
                    fset(work, iwrk + k - 1, fget(vl.*, k, i)*fget(vl.*, k, i) + fget(vl.*, k, i + 1)*fget(vl.*, k, i + 1));
                }
                scratch := _geev_subvector(work.*, iwrk, n);
                k = _geev_iamax_1based(n, scratch, 1, RT);
                _geev_rotg(fget(vl.*, k, i), fget(vl.*, k, i + 1), *cs, *sn, *r, RT);
                _geev_rot(n, *col1, 1, *col2, 1, cs, sn, RT);
                fset(vl, k, i + 1, zero_r);
                i += 2;
            } else {
                i += 1;
            }
        }
    }

    if wantvr {
        bal_scale := _geev_subvector(work.*, ibal, n);
        #if RT == float32 {
            sgebak(#char "B", #char "R", n, ilo, ihi, bal_scale, n, vr, ldvr, *ierr);
        } else {
            dgebak(#char "B", #char "R", n, ilo, ihi, bal_scale, n, vr, ldvr, *ierr);
        }
        i = 1;
        while i <= n {
            if fget(wi.*, i) == zero_r {
                col := _geev_matrix_column(vr.*, ldvr, 1, i, n);
                scl = one_r / _geev_nrm2(n, col, 1, RT, RT);
                _geev_scal_real(n, scl, *col, 1, RT, RT);
                i += 1;
            } else if fget(wi.*, i) > zero_r {
                col1 := _geev_matrix_column(vr.*, ldvr, 1, i, n);
                col2 := _geev_matrix_column(vr.*, ldvr, 1, i + 1, n);
                scl = one_r / _geev_lapy2(_geev_nrm2(n, col1, 1, RT, RT), _geev_nrm2(n, col2, 1, RT, RT), RT);
                _geev_scal_real(n, scl, *col1, 1, RT, RT);
                _geev_scal_real(n, scl, *col2, 1, RT, RT);
                for k:1..n {
                    fset(work, iwrk + k - 1, fget(vr.*, k, i)*fget(vr.*, k, i) + fget(vr.*, k, i + 1)*fget(vr.*, k, i + 1));
                }
                scratch := _geev_subvector(work.*, iwrk, n);
                k = _geev_iamax_1based(n, scratch, 1, RT);
                _geev_rotg(fget(vr.*, k, i), fget(vr.*, k, i + 1), *cs, *sn, *r, RT);
                _geev_rot(n, *col1, 1, *col2, 1, cs, sn, RT);
                fset(vr, k, i + 1, zero_r);
                i += 2;
            } else {
                i += 1;
            }
        }
    }

end_scaling:
    if scalea {
        factor := anrm / cscale;
        _geev_scale_real_segment(wr, info.* + 1, n - info.*, factor, RT);
        _geev_scale_real_segment(wi, info.* + 1, n - info.*, factor, RT);
        if info.* > 0 {
            _geev_scale_real_segment(wr, 1, ilo - 1, factor, RT);
            _geev_scale_real_segment(wi, 1, ilo - 1, factor, RT);
        }
    }

    fset(work, 1, ncast(RT, maxwrk));
}

_geev_complex :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    w: *$W/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
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
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == C32 { "CGEEV " } else { "ZGEEV " };
    gehrd_name :: #run if MT == C32 { "CGEHRD" } else { "ZGEHRD" };
    unghr_name :: #run if MT == C32 { "CUNGHR" } else { "ZUNGHR" };

    wantvl, wantvr, lquery, scalea : bool = ---;
    side : u8 = ---;
    hswork, i, ibal, ierr, ihi, ilo, irwork, itau, iwrk, k, lwork_trevc, maxwrk, minwrk, nout : int = ---;
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

    if !wantvl && !lsame(jobvl, #char "N") {
        info.* = -1;
    } else if !wantvr && !lsame(jobvr, #char "N") {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ldvl < 1 || (wantvl && ldvl < n) {
        info.* = -8;
    } else if ldvr < 1 || (wantvr && ldvr < n) {
        info.* = -10;
    }

    if info.* == 0 {
        if n == 0 {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = n + n*ilaenv(1, gehrd_name, " "s, n, 1, n, 0);
            minwrk = 2*n;
            if wantvl {
                maxwrk = max(maxwrk, n + (n - 1)*ilaenv(1, unghr_name, " "s, n, 1, n, -1));
                #if MT == C32 {
                    ctrevc3(#char "L", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, rwork, -1, *ierr);
                    chseqr(#char "S", #char "V", n, 1, n, a, lda, w, vl, ldvl, hseqr_query, -1, *info);
                } else {
                    ztrevc3(#char "L", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, rwork, -1, *ierr);
                    zhseqr(#char "S", #char "V", n, 1, n, a, lda, w, vl, ldvl, hseqr_query, -1, *info);
                }
                lwork_trevc = cast(int) real(trevc_query[0]);
                hswork = cast(int) real(hseqr_query[0]);
                maxwrk = max(maxwrk, hswork, lwork_trevc, minwrk);
            } else if wantvr {
                maxwrk = max(maxwrk, n + (n - 1)*ilaenv(1, unghr_name, " "s, n, 1, n, -1));
                #if MT == C32 {
                    ctrevc3(#char "R", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, rwork, -1, *ierr);
                    chseqr(#char "S", #char "V", n, 1, n, a, lda, w, vr, ldvr, hseqr_query, -1, *info);
                } else {
                    ztrevc3(#char "R", #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, trevc_query, -1, rwork, -1, *ierr);
                    zhseqr(#char "S", #char "V", n, 1, n, a, lda, w, vr, ldvr, hseqr_query, -1, *info);
                }
                lwork_trevc = cast(int) real(trevc_query[0]);
                hswork = cast(int) real(hseqr_query[0]);
                maxwrk = max(maxwrk, hswork, lwork_trevc, minwrk);
            } else {
                #if MT == C32 {
                    chseqr(#char "E", #char "N", n, 1, n, a, lda, w, vr, ldvr, hseqr_query, -1, *info);
                } else {
                    zhseqr(#char "E", #char "N", n, 1, n, a, lda, w, vr, ldvr, hseqr_query, -1, *info);
                }
                hswork = cast(int) real(hseqr_query[0]);
                maxwrk = max(maxwrk, hswork, minwrk);
            }
        }
        fset(work, 1, _geev_complex_from_real(ncast(RT, maxwrk), RT, MT));
        if (lwork < minwrk) && !lquery {
            info.* = -12;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 then return;

    eps = _geev_lamch_precision(RT);
    smlnum = _geev_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    anrm = _geev_lange(#char "M", n, n, a.*, lda, *dum, RT, MT);
    scalea = false;
    if anrm > zero_r && anrm < smlnum {
        scalea = true;
        cscale = smlnum;
    } else if anrm > bignum {
        scalea = true;
        cscale = bignum;
    }
    if scalea {
        _geev_lascl(#char "G", 0, 0, anrm, cscale, n, n, a, lda, *ierr, RT, MT);
    }

    ibal = 1;
    bal_scale := _geev_subvector(rwork.*, ibal, n);
    #if MT == C32 {
        cgebal(#char "B", n, a, lda, *ilo, *ihi, bal_scale, *ierr);
    } else {
        zgebal(#char "B", n, a, lda, *ilo, *ihi, bal_scale, *ierr);
    }

    itau = 1;
    iwrk = itau + n;
    tau := _geev_subvector(work.*, itau, n);
    wrk := _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if MT == C32 {
        cgehrd(n, ilo, ihi, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    } else {
        zgehrd(n, ilo, ihi, a, lda, tau, wrk, lwork - iwrk + 1, *ierr);
    }

    if wantvl {
        side = #char "L";
        _geev_lacpy(#char "L", n, n, a.*, lda, vl, ldvl, MT);
        #if MT == C32 {
            cunghr(n, ilo, ihi, vl, ldvl, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            zunghr(n, ilo, ihi, vl, ldvl, tau, wrk, lwork - iwrk + 1, *ierr);
        }
        iwrk = itau;
        wrk = _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if MT == C32 {
            chseqr(#char "S", #char "V", n, ilo, ihi, a, lda, w, vl, ldvl, wrk, lwork - iwrk + 1, info);
        } else {
            zhseqr(#char "S", #char "V", n, ilo, ihi, a, lda, w, vl, ldvl, wrk, lwork - iwrk + 1, info);
        }
        if wantvr {
            side = #char "B";
            _geev_lacpy(#char "F", n, n, vl.*, ldvl, vr, ldvr, MT);
        }
    } else if wantvr {
        side = #char "R";
        _geev_lacpy(#char "L", n, n, a.*, lda, vr, ldvr, MT);
        #if MT == C32 {
            cunghr(n, ilo, ihi, vr, ldvr, tau, wrk, lwork - iwrk + 1, *ierr);
        } else {
            zunghr(n, ilo, ihi, vr, ldvr, tau, wrk, lwork - iwrk + 1, *ierr);
        }
        iwrk = itau;
        wrk = _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if MT == C32 {
            chseqr(#char "S", #char "V", n, ilo, ihi, a, lda, w, vr, ldvr, wrk, lwork - iwrk + 1, info);
        } else {
            zhseqr(#char "S", #char "V", n, ilo, ihi, a, lda, w, vr, ldvr, wrk, lwork - iwrk + 1, info);
        }
    } else {
        iwrk = itau;
        wrk = _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if MT == C32 {
            chseqr(#char "E", #char "N", n, ilo, ihi, a, lda, w, vr, ldvr, wrk, lwork - iwrk + 1, info);
        } else {
            zhseqr(#char "E", #char "N", n, ilo, ihi, a, lda, w, vr, ldvr, wrk, lwork - iwrk + 1, info);
        }
    }

    if info.* != 0 goto complex_end_scaling;

    if wantvl || wantvr {
        irwork = ibal + n;
        rwrk := _geev_subvector(rwork.*, irwork, n);
        wrk = _geev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
        #if MT == C32 {
            ctrevc3(side, #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, wrk, lwork - iwrk + 1, rwrk, n, *ierr);
        } else {
            ztrevc3(side, #char "B", select_dummy, n, a, lda, vl, ldvl, vr, ldvr, n, *nout, wrk, lwork - iwrk + 1, rwrk, n, *ierr);
        }
    }

    if wantvl {
        #if MT == C32 {
            cgebak(#char "B", #char "L", n, ilo, ihi, bal_scale, n, vl, ldvl, *ierr);
        } else {
            zgebak(#char "B", #char "L", n, ilo, ihi, bal_scale, n, vl, ldvl, *ierr);
        }
        irwork = ibal + n;
        for i:1..n {
            col := _geev_matrix_column(vl.*, ldvl, 1, i, n);
            scl = one_r / _geev_nrm2(n, col, 1, RT, MT);
            _geev_scal_real(n, scl, *col, 1, RT, MT);
            for k:1..n {
                fset(rwork, irwork + k - 1, real(fget(vl.*, k, i))*real(fget(vl.*, k, i)) + imag(fget(vl.*, k, i))*imag(fget(vl.*, k, i)));
            }
            scratch := _geev_subvector(rwork.*, irwork, n);
            k = _geev_iamax_1based(n, scratch, 1, RT);
            tmp = conjugate(fget(vl.*, k, i)) / sqrt(fget(rwork.*, irwork + k - 1));
            _geev_scale_complex(n, tmp, *col, 1, MT);
            fset(vl, k, i, cpx(real(fget(vl.*, k, i)), zero_r));
        }
    }

    if wantvr {
        #if MT == C32 {
            cgebak(#char "B", #char "R", n, ilo, ihi, bal_scale, n, vr, ldvr, *ierr);
        } else {
            zgebak(#char "B", #char "R", n, ilo, ihi, bal_scale, n, vr, ldvr, *ierr);
        }
        irwork = ibal + n;
        for i:1..n {
            col := _geev_matrix_column(vr.*, ldvr, 1, i, n);
            scl = one_r / _geev_nrm2(n, col, 1, RT, MT);
            _geev_scal_real(n, scl, *col, 1, RT, MT);
            for k:1..n {
                fset(rwork, irwork + k - 1, real(fget(vr.*, k, i))*real(fget(vr.*, k, i)) + imag(fget(vr.*, k, i))*imag(fget(vr.*, k, i)));
            }
            scratch := _geev_subvector(rwork.*, irwork, n);
            k = _geev_iamax_1based(n, scratch, 1, RT);
            tmp = conjugate(fget(vr.*, k, i)) / sqrt(fget(rwork.*, irwork + k - 1));
            _geev_scale_complex(n, tmp, *col, 1, MT);
            fset(vr, k, i, cpx(real(fget(vr.*, k, i)), zero_r));
        }
    }

complex_end_scaling:
    if scalea {
        factor := anrm / cscale;
        _geev_scale_complex_segment(w, info.* + 1, n - info.*, factor, RT, MT);
        if info.* > 0 {
            _geev_scale_complex_segment(w, 1, ilo - 1, factor, RT, MT);
        }
    }

    fset(work, 1, _geev_complex_from_real(ncast(RT, maxwrk), RT, MT));
}
