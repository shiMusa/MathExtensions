// ? LAPACK *ggev()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggev/cggev/dggev/zggev
// *GGEV computes the generalized eigenvalues and, optionally, the
// left and/or right generalized eigenvectors of a matrix pair.

sggev :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D,
    alphai: *$ALPHAI/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int
) {
    _ggev_real(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info, blocked=false, RT=float32);
}

dggev :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D,
    alphai: *$ALPHAI/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int
) {
    _ggev_real(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info, blocked=false, RT=float64);
}

cggev :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alpha: *$ALPHA/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _ggev_complex(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info, blocked=false, RT=float32, MT=C32);
}

zggev :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alpha: *$ALPHA/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _ggev_complex(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info, blocked=false, RT=float64, MT=C64);
}

#scope_file
_ggev_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_ggev_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_ggev_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _ggev_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_ggev_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_ggev_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return ncast(MT, x);
    } else {
        return cpx(x, zero(RT));
    }
}

_ggev_query_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_ggev_lange :: (
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

_ggev_lacpy :: (
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

_ggev_laset :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    beta: MT,
    a: *$A/ArrayType2D,
    $$lda: int,
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

_ggev_lascl_matrix :: (
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

_ggev_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_ggev_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "E");
    } else {
        return dlamch(#char "E");
    }
}

_ggev_lamch_base :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "B");
    } else {
        return dlamch(#char "B");
    }
}

_ggev_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_ggev_scale_real_segment :: (
    x: *$X/ArrayType1D,
    $$offset: int,
    $$count: int,
    scale: $RT,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    if count <= 0 then return;
    for i:offset..offset + count - 1 {
        fset(x, i, fget(x.*, i)*scale);
    }
}

_ggev_scale_complex_segment :: (
    x: *$X/ArrayType1D,
    $$offset: int,
    $$count: int,
    scale: $RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    if count <= 0 then return;
    scale_t := cpx(scale, zero(RT));
    for i:offset..offset + count - 1 {
        fset(x, i, fget(x.*, i)*scale_t);
    }
}

_ggev_abs1 :: (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return abs(real(x)) + abs(imag(x));
}

_ggev_real :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D,
    alphai: *$ALPHAI/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int,
    $$blocked: bool,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           B.data_type == RT &&
           ALPHAR.data_type == RT &&
           ALPHAI.data_type == RT &&
           BETA.data_type == RT &&
           VL.data_type == RT &&
           VR.data_type == RT &&
           WORK.data_type == RT;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if RT == float32 { if blocked then "SGGEV3 " else "SGGEV " } else { if blocked then "DGGEV3 " else "DGGEV " };
    geqrf_name :: #run if RT == float32 { "SGEQRF" } else { "DGEQRF" };
    ormqr_name :: #run if RT == float32 { "SORMQR" } else { "DORMQR" };
    orgqr_name :: #run if RT == float32 { "SORGQR" } else { "DORGQR" };

    ilascl, ilbscl, ilv, ilvl, ilvr, lquery : bool = ---;
    chtemp : u8 = ---;
    anrm, anrmto, bignum, bnrm, bnrmto, eps, smlnum, temp : RT = ---;
    ierr, ihi, ijobvl, ijobvr, ileft, ilo, in, iright, irows, icols, itau, iwrk, jc, jr, lwkmin, lwkopt, maxwrk, minwrk : int = ---;
    ldumma: [1] bool;
    query: [1] RT;

    info.* = 0;
    lquery = lwork == -1;
    if lsame(jobvl, #char "N") { ijobvl = 1; ilvl = false; } else if lsame(jobvl, #char "V") { ijobvl = 2; ilvl = true; } else { ijobvl = -1; ilvl = false; }
    if lsame(jobvr, #char "N") { ijobvr = 1; ilvr = false; } else if lsame(jobvr, #char "V") { ijobvr = 2; ilvr = true; } else { ijobvr = -1; ilvr = false; }
    ilv = ilvl || ilvr;

    if ijobvl <= 0 { info.* = -1; }
    else if ijobvr <= 0 { info.* = -2; }
    else if n < 0 { info.* = -3; }
    else if lda < max(1, n) { info.* = -5; }
    else if ldb < max(1, n) { info.* = -7; }
    else if ldvl < 1 || (ilvl && ldvl < n) { info.* = -12; }
    else if ldvr < 1 || (ilvr && ldvr < n) { info.* = -14; }

    if info.* == 0 {
        if blocked {
            lwkmin = max(1, 8*n);
            #if RT == float32 { sgeqrf(n, n, b, ldb, query, query, -1, *ierr); } else { dgeqrf(n, n, b, ldb, query, query, -1, *ierr); }
            lwkopt = max(lwkmin, 3*n + cast(int) query[0]);
            #if RT == float32 { sormqr(#char "L", #char "T", n, n, n, b, ldb, query, a, lda, query, -1, *ierr); } else { dormqr(#char "L", #char "T", n, n, n, b, ldb, query, a, lda, query, -1, *ierr); }
            lwkopt = max(lwkopt, 3*n + cast(int) query[0]);
            #if RT == float32 { sgghd3(jobvl, jobvr, n, 1, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, query, -1, *ierr); } else { dgghd3(jobvl, jobvr, n, 1, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, query, -1, *ierr); }
            lwkopt = max(lwkopt, 3*n + cast(int) query[0]);
            if ilvl {
                #if RT == float32 { sorgqr(n, n, n, vl, ldvl, query, query, -1, *ierr); } else { dorgqr(n, n, n, vl, ldvl, query, query, -1, *ierr); }
                lwkopt = max(lwkopt, 3*n + cast(int) query[0]);
                #if RT == float32 { slaqz0(#char "S", jobvl, jobvr, n, 1, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, query, -1, 0, *ierr); } else { dlaqz0(#char "S", jobvl, jobvr, n, 1, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, query, -1, 0, *ierr); }
                lwkopt = max(lwkopt, 2*n + cast(int) query[0]);
            } else {
                #if RT == float32 { slaqz0(#char "E", jobvl, jobvr, n, 1, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, query, -1, 0, *ierr); } else { dlaqz0(#char "E", jobvl, jobvr, n, 1, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, query, -1, 0, *ierr); }
                lwkopt = max(lwkopt, 2*n + cast(int) query[0]);
            }
            fset(work, 1, ifx n == 0 then one_r else _ggev_roundup_lwork(lwkopt, RT));
            if lwork < lwkmin && !lquery info.* = -16;
        } else {
            minwrk = max(1, 8*n);
            maxwrk = max(1, n*(7 + ilaenv(1, geqrf_name, " "s, n, 1, n, 0)));
            maxwrk = max(maxwrk, n*(7 + ilaenv(1, ormqr_name, " "s, n, 1, n, 0)));
            if ilvl maxwrk = max(maxwrk, n*(7 + ilaenv(1, orgqr_name, " "s, n, 1, n, -1)));
            fset(work, 1, _ggev_roundup_lwork(maxwrk, RT));
            if lwork < minwrk && !lquery info.* = -16;
        }
    }

    if info.* != 0 { xerbla(routine_name, -info.*); return; }
    if lquery || n == 0 then return;

    eps = _ggev_lamch_precision(RT);
    smlnum = _ggev_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    anrm = _ggev_lange(#char "M", n, n, a.*, lda, work, RT, RT);
    ilascl = false;
    if anrm > zero_r && anrm < smlnum { anrmto = smlnum; ilascl = true; } else if anrm > bignum { anrmto = bignum; ilascl = true; }
    if ilascl _ggev_lascl_matrix(#char "G", 0, 0, anrm, anrmto, n, n, a, lda, *ierr, RT, RT);

    bnrm = _ggev_lange(#char "M", n, n, b.*, ldb, work, RT, RT);
    ilbscl = false;
    if bnrm > zero_r && bnrm < smlnum { bnrmto = smlnum; ilbscl = true; } else if bnrm > bignum { bnrmto = bignum; ilbscl = true; }
    if ilbscl _ggev_lascl_matrix(#char "G", 0, 0, bnrm, bnrmto, n, n, b, ldb, *ierr, RT, RT);

    ileft = 1; iright = n + 1; iwrk = iright + n;
    lscale := _ggev_subvector(work.*, ileft, n);
    rscale := _ggev_subvector(work.*, iright, n);
    bal_work := _ggev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    #if RT == float32 { sggbal(#char "P", n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, bal_work, *ierr); } else { dggbal(#char "P", n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, bal_work, *ierr); }

    irows = ihi + 1 - ilo;
    icols = ifx ilv then n + 1 - ilo else irows;
    itau = iwrk;
    iwrk = itau + irows;
    tau := _ggev_subvector(work.*, itau, max(0, irows));
    wrk := _ggev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    b_raw := _ggev_matrix_storage(b.*, ldb, ilo, ilo, irows, icols);
    b_sub := array2d_view(b_raw, ldb, row_major=false);
    #if RT == float32 { sgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr); } else { dgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr); }
    a_raw := _ggev_matrix_storage(a.*, lda, ilo, ilo, irows, icols);
    a_sub := array2d_view(a_raw, lda, row_major=false);
    #if RT == float32 { sormqr(#char "L", #char "T", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr); } else { dormqr(#char "L", #char "T", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr); }

    if ilvl {
        _ggev_laset(#char "F", n, n, zero_r, one_r, vl, ldvl, RT);
        if irows > 1 {
            b_lower_raw := _ggev_matrix_storage(b.*, ldb, ilo + 1, ilo, irows - 1, irows - 1);
            vl_lower_raw := _ggev_matrix_storage(vl.*, ldvl, ilo + 1, ilo, irows - 1, irows - 1);
            b_lower := array2d_view(b_lower_raw, ldb, row_major=false);
            vl_lower := array2d_view(vl_lower_raw, ldvl, row_major=false);
            _ggev_lacpy(#char "L", irows - 1, irows - 1, *b_lower, ldb, vl_lower, ldvl, RT);
        }
        vl_raw := _ggev_matrix_storage(vl.*, ldvl, ilo, ilo, irows, irows);
        vl_sub := array2d_view(vl_raw, ldvl, row_major=false);
        #if RT == float32 { sorgqr(irows, irows, irows, *vl_sub, ldvl, tau, wrk, lwork - iwrk + 1, *ierr); } else { dorgqr(irows, irows, irows, *vl_sub, ldvl, tau, wrk, lwork - iwrk + 1, *ierr); }
    }
    if ilvr _ggev_laset(#char "F", n, n, zero_r, one_r, vr, ldvr, RT);

    if blocked {
        if ilv {
            #if RT == float32 { sgghd3(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); } else { dgghd3(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); }
        } else {
            sub_a_raw := _ggev_matrix_storage(a.*, lda, ilo, ilo, irows, irows);
            sub_b_raw := _ggev_matrix_storage(b.*, ldb, ilo, ilo, irows, irows);
            sub_a := array2d_view(sub_a_raw, lda, row_major=false);
            sub_b := array2d_view(sub_b_raw, ldb, row_major=false);
            #if RT == float32 { sgghd3(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); } else { dgghd3(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); }
        }
    } else {
        if ilv {
            #if RT == float32 { sgghrd(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, *ierr); } else { dgghrd(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, *ierr); }
        } else {
            sub_a_raw := _ggev_matrix_storage(a.*, lda, ilo, ilo, irows, irows);
            sub_b_raw := _ggev_matrix_storage(b.*, ldb, ilo, ilo, irows, irows);
            sub_a := array2d_view(sub_a_raw, lda, row_major=false);
            sub_b := array2d_view(sub_b_raw, ldb, row_major=false);
            #if RT == float32 { sgghrd(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, *ierr); } else { dgghrd(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, *ierr); }
        }
    }

    iwrk = itau;
    wrk = _ggev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    chtemp = ifx ilv then #char "S" else #char "E";
    if blocked {
        #if RT == float32 { slaqz0(chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, 0, *ierr); } else { dlaqz0(chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, 0, *ierr); }
    } else {
        #if RT == float32 { shgeqz(chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); } else { dhgeqz(chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); }
    }
    if ierr != 0 {
        if ierr > 0 && ierr <= n { info.* = ierr; } else if ierr > n && ierr <= 2*n { info.* = ierr - n; } else { info.* = n + 1; }
        fset(work, 1, ifx blocked then _ggev_roundup_lwork(lwkopt, RT) else _ggev_roundup_lwork(maxwrk, RT));
        return;
    }

    if ilv {
        if ilvl { if ilvr then chtemp = #char "B" else chtemp = #char "L"; } else { chtemp = #char "R"; }
        #if RT == float32 { stgevc(chtemp, #char "B", ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, *in, wrk, *ierr); } else { dtgevc(chtemp, #char "B", ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, *in, wrk, *ierr); }
        if ierr != 0 { info.* = n + 2; }
        else {
            if ilvl {
                #if RT == float32 { sggbak(#char "P", #char "L", n, ilo, ihi, lscale, rscale, n, vl, ldvl, *ierr); } else { dggbak(#char "P", #char "L", n, ilo, ihi, lscale, rscale, n, vl, ldvl, *ierr); }
                for jc:1..n {
                    if fget(alphai.*, jc) < zero_r then continue;
                    temp = zero_r;
                    if fget(alphai.*, jc) == zero_r { for jr:1..n temp = max(temp, abs(fget(vl.*, jr, jc))); }
                    else { for jr:1..n temp = max(temp, abs(fget(vl.*, jr, jc)) + abs(fget(vl.*, jr, jc + 1))); }
                    if temp < smlnum then continue;
                    temp = one_r / temp;
                    if fget(alphai.*, jc) == zero_r { for jr:1..n fset(vl, jr, jc, fget(vl.*, jr, jc)*temp); }
                    else { for jr:1..n { fset(vl, jr, jc, fget(vl.*, jr, jc)*temp); fset(vl, jr, jc + 1, fget(vl.*, jr, jc + 1)*temp); } }
                }
            }
            if ilvr {
                #if RT == float32 { sggbak(#char "P", #char "R", n, ilo, ihi, lscale, rscale, n, vr, ldvr, *ierr); } else { dggbak(#char "P", #char "R", n, ilo, ihi, lscale, rscale, n, vr, ldvr, *ierr); }
                for jc:1..n {
                    if fget(alphai.*, jc) < zero_r then continue;
                    temp = zero_r;
                    if fget(alphai.*, jc) == zero_r { for jr:1..n temp = max(temp, abs(fget(vr.*, jr, jc))); }
                    else { for jr:1..n temp = max(temp, abs(fget(vr.*, jr, jc)) + abs(fget(vr.*, jr, jc + 1))); }
                    if temp < smlnum then continue;
                    temp = one_r / temp;
                    if fget(alphai.*, jc) == zero_r { for jr:1..n fset(vr, jr, jc, fget(vr.*, jr, jc)*temp); }
                    else { for jr:1..n { fset(vr, jr, jc, fget(vr.*, jr, jc)*temp); fset(vr, jr, jc + 1, fget(vr.*, jr, jc + 1)*temp); } }
                }
            }
        }
    }

    if ilascl { _ggev_scale_real_segment(alphar, 1, n, anrm / anrmto, RT); _ggev_scale_real_segment(alphai, 1, n, anrm / anrmto, RT); }
    if ilbscl _ggev_scale_real_segment(beta, 1, n, bnrm / bnrmto, RT);
    fset(work, 1, ifx blocked then _ggev_roundup_lwork(lwkopt, RT) else _ggev_roundup_lwork(maxwrk, RT));
}

_ggev_complex :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alpha: *$ALPHA/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $$blocked: bool,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT &&
           ALPHA.data_type == MT &&
           BETA.data_type == MT &&
           VL.data_type == MT &&
           VR.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    zero_t :: #run cpx(zero(RT), zero(RT));
    one_t :: #run cpx(one(RT), zero(RT));
    routine_name :: #run if RT == float32 { if blocked then "CGGEV3 " else "CGGEV " } else { if blocked then "ZGGEV3 " else "ZGGEV " };
    geqrf_name :: #run if RT == float32 { "CGEQRF" } else { "ZGEQRF" };
    unmqr_name :: #run if RT == float32 { "CUNMQR" } else { "ZUNMQR" };
    ungqr_name :: #run if RT == float32 { "CUNGQR" } else { "ZUNGQR" };

    ilascl, ilbscl, ilv, ilvl, ilvr, lquery : bool = ---;
    chtemp : u8 = ---;
    anrm, anrmto, bignum, bnrm, bnrmto, eps, smlnum, temp : RT = ---;
    ierr, ihi, ijobvl, ijobvr, ileft, ilo, in, iright, irwrk, irows, icols, itau, iwrk, jc, jr, lwkmin, lwkopt : int = ---;
    ldumma: [1] bool;
    query: [1] MT;

    info.* = 0;
    lquery = lwork == -1;
    if lsame(jobvl, #char "N") { ijobvl = 1; ilvl = false; } else if lsame(jobvl, #char "V") { ijobvl = 2; ilvl = true; } else { ijobvl = -1; ilvl = false; }
    if lsame(jobvr, #char "N") { ijobvr = 1; ilvr = false; } else if lsame(jobvr, #char "V") { ijobvr = 2; ilvr = true; } else { ijobvr = -1; ilvr = false; }
    ilv = ilvl || ilvr;

    if ijobvl <= 0 { info.* = -1; }
    else if ijobvr <= 0 { info.* = -2; }
    else if n < 0 { info.* = -3; }
    else if lda < max(1, n) { info.* = -5; }
    else if ldb < max(1, n) { info.* = -7; }
    else if ldvl < 1 || (ilvl && ldvl < n) { info.* = -11; }
    else if ldvr < 1 || (ilvr && ldvr < n) { info.* = -13; }

    if info.* == 0 {
        lwkmin = max(1, 2*n);
        if blocked {
            #if RT == float32 { cgeqrf(n, n, b, ldb, query, query, -1, *ierr); } else { zgeqrf(n, n, b, ldb, query, query, -1, *ierr); }
            lwkopt = max(lwkmin, n + _ggev_query_to_int(query[0], RT, MT));
            #if RT == float32 { cunmqr(#char "L", #char "C", n, n, n, b, ldb, query, a, lda, query, -1, *ierr); } else { zunmqr(#char "L", #char "C", n, n, n, b, ldb, query, a, lda, query, -1, *ierr); }
            lwkopt = max(lwkopt, n + _ggev_query_to_int(query[0], RT, MT));
            if ilvl {
                #if RT == float32 { cungqr(n, n, n, vl, ldvl, query, query, -1, *ierr); } else { zungqr(n, n, n, vl, ldvl, query, query, -1, *ierr); }
                lwkopt = max(lwkopt, n + _ggev_query_to_int(query[0], RT, MT));
            }
            if ilv {
                #if RT == float32 { cgghd3(jobvl, jobvr, n, 1, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, query, -1, *ierr); } else { zgghd3(jobvl, jobvr, n, 1, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, query, -1, *ierr); }
                lwkopt = max(lwkopt, n + _ggev_query_to_int(query[0], RT, MT));
                #if RT == float32 { claqz0(#char "S", jobvl, jobvr, n, 1, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, query, -1, rwork, 0, *ierr); } else { zlaqz0(#char "S", jobvl, jobvr, n, 1, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, query, -1, rwork, 0, *ierr); }
                lwkopt = max(lwkopt, n + _ggev_query_to_int(query[0], RT, MT));
            } else {
                #if RT == float32 { cgghd3(#char "N", #char "N", n, 1, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, query, -1, *ierr); } else { zgghd3(#char "N", #char "N", n, 1, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, query, -1, *ierr); }
                lwkopt = max(lwkopt, n + _ggev_query_to_int(query[0], RT, MT));
                #if RT == float32 { claqz0(#char "E", jobvl, jobvr, n, 1, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, query, -1, rwork, 0, *ierr); } else { zlaqz0(#char "E", jobvl, jobvr, n, 1, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, query, -1, rwork, 0, *ierr); }
                lwkopt = max(lwkopt, n + _ggev_query_to_int(query[0], RT, MT));
            }
        } else {
            lwkopt = max(1, n + n*ilaenv(1, geqrf_name, " "s, n, 1, n, 0));
            lwkopt = max(lwkopt, n + n*ilaenv(1, unmqr_name, " "s, n, 1, n, 0));
            if ilvl lwkopt = max(lwkopt, n + n*ilaenv(1, ungqr_name, " "s, n, 1, n, -1));
        }
        fset(work, 1, ifx n == 0 then one_t else _ggev_complex_from_real(_ggev_roundup_lwork(lwkopt, RT), RT, MT));
        if lwork < lwkmin && !lquery info.* = -15;
    }

    if info.* != 0 { xerbla(routine_name, -info.*); return; }
    if lquery || n == 0 then return;

    eps = _ggev_lamch_epsilon(RT) * _ggev_lamch_base(RT);
    smlnum = _ggev_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    anrm = _ggev_lange(#char "M", n, n, a.*, lda, rwork, RT, MT);
    ilascl = false;
    if anrm > zero_r && anrm < smlnum { anrmto = smlnum; ilascl = true; } else if anrm > bignum { anrmto = bignum; ilascl = true; }
    if ilascl _ggev_lascl_matrix(#char "G", 0, 0, anrm, anrmto, n, n, a, lda, *ierr, RT, MT);
    bnrm = _ggev_lange(#char "M", n, n, b.*, ldb, rwork, RT, MT);
    ilbscl = false;
    if bnrm > zero_r && bnrm < smlnum { bnrmto = smlnum; ilbscl = true; } else if bnrm > bignum { bnrmto = bignum; ilbscl = true; }
    if ilbscl _ggev_lascl_matrix(#char "G", 0, 0, bnrm, bnrmto, n, n, b, ldb, *ierr, RT, MT);

    ileft = 1; iright = n + 1; irwrk = iright + n;
    lscale := _ggev_subvector(rwork.*, ileft, n);
    rscale := _ggev_subvector(rwork.*, iright, n);
    bal_work := _ggev_subvector(rwork.*, irwrk, max(0, rwork.count - irwrk + 1));
    #if RT == float32 { cggbal(#char "P", n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, bal_work, *ierr); } else { zggbal(#char "P", n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, bal_work, *ierr); }

    irows = ihi + 1 - ilo; icols = ifx ilv then n + 1 - ilo else irows;
    itau = 1; iwrk = itau + irows;
    tau := _ggev_subvector(work.*, itau, max(0, irows));
    wrk := _ggev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    rwrk := _ggev_subvector(rwork.*, irwrk, max(0, rwork.count - irwrk + 1));
    b_raw := _ggev_matrix_storage(b.*, ldb, ilo, ilo, irows, icols);
    b_sub := array2d_view(b_raw, ldb, row_major=false);
    #if RT == float32 { cgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr); } else { zgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr); }
    a_raw := _ggev_matrix_storage(a.*, lda, ilo, ilo, irows, icols);
    a_sub := array2d_view(a_raw, lda, row_major=false);
    #if RT == float32 { cunmqr(#char "L", #char "C", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr); } else { zunmqr(#char "L", #char "C", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr); }

    if ilvl {
        _ggev_laset(#char "F", n, n, zero_t, one_t, vl, ldvl, MT);
        if irows > 1 {
            b_lower_raw := _ggev_matrix_storage(b.*, ldb, ilo + 1, ilo, irows - 1, irows - 1);
            vl_lower_raw := _ggev_matrix_storage(vl.*, ldvl, ilo + 1, ilo, irows - 1, irows - 1);
            b_lower := array2d_view(b_lower_raw, ldb, row_major=false);
            vl_lower := array2d_view(vl_lower_raw, ldvl, row_major=false);
            _ggev_lacpy(#char "L", irows - 1, irows - 1, *b_lower, ldb, vl_lower, ldvl, MT);
        }
        vl_raw := _ggev_matrix_storage(vl.*, ldvl, ilo, ilo, irows, irows);
        vl_sub := array2d_view(vl_raw, ldvl, row_major=false);
        #if RT == float32 { cungqr(irows, irows, irows, *vl_sub, ldvl, tau, wrk, lwork - iwrk + 1, *ierr); } else { zungqr(irows, irows, irows, *vl_sub, ldvl, tau, wrk, lwork - iwrk + 1, *ierr); }
    }
    if ilvr _ggev_laset(#char "F", n, n, zero_t, one_t, vr, ldvr, MT);

    if blocked {
        if ilv {
            #if RT == float32 { cgghd3(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); } else { zgghd3(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); }
        } else {
            sub_a_raw := _ggev_matrix_storage(a.*, lda, ilo, ilo, irows, irows);
            sub_b_raw := _ggev_matrix_storage(b.*, ldb, ilo, ilo, irows, irows);
            sub_a := array2d_view(sub_a_raw, lda, row_major=false);
            sub_b := array2d_view(sub_b_raw, ldb, row_major=false);
            #if RT == float32 { cgghd3(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); } else { zgghd3(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, *ierr); }
        }
    } else {
        if ilv {
            #if RT == float32 { cgghrd(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, *ierr); } else { zgghrd(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, *ierr); }
        } else {
            sub_a_raw := _ggev_matrix_storage(a.*, lda, ilo, ilo, irows, irows);
            sub_b_raw := _ggev_matrix_storage(b.*, ldb, ilo, ilo, irows, irows);
            sub_a := array2d_view(sub_a_raw, lda, row_major=false);
            sub_b := array2d_view(sub_b_raw, ldb, row_major=false);
            #if RT == float32 { cgghrd(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, *ierr); } else { zgghrd(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, *ierr); }
        }
    }

    iwrk = itau;
    wrk = _ggev_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    rwrk = _ggev_subvector(rwork.*, irwrk, max(0, rwork.count - irwrk + 1));
    chtemp = ifx ilv then #char "S" else #char "E";
    if blocked {
        #if RT == float32 { claqz0(chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, rwrk, 0, *ierr); } else { zlaqz0(chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, rwrk, 0, *ierr); }
    } else {
        #if RT == float32 { chgeqz(chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, rwrk, *ierr); } else { zhgeqz(chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, rwrk, *ierr); }
    }
    if ierr != 0 {
        if ierr > 0 && ierr <= n { info.* = ierr; } else if ierr > n && ierr <= 2*n { info.* = ierr - n; } else { info.* = n + 1; }
        fset(work, 1, _ggev_complex_from_real(_ggev_roundup_lwork(lwkopt, RT), RT, MT));
        return;
    }

    if ilv {
        if ilvl { if ilvr then chtemp = #char "B" else chtemp = #char "L"; } else { chtemp = #char "R"; }
        #if RT == float32 { ctgevc(chtemp, #char "B", ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, *in, wrk, rwrk, *ierr); } else { ztgevc(chtemp, #char "B", ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, *in, wrk, rwrk, *ierr); }
        if ierr != 0 { info.* = n + 2; }
        else {
            if ilvl {
                #if RT == float32 { cggbak(#char "P", #char "L", n, ilo, ihi, lscale, rscale, n, vl, ldvl, *ierr); } else { zggbak(#char "P", #char "L", n, ilo, ihi, lscale, rscale, n, vl, ldvl, *ierr); }
                for jc:1..n {
                    temp = zero_r;
                    for jr:1..n temp = max(temp, _ggev_abs1(fget(vl.*, jr, jc), RT, MT));
                    if temp < smlnum then continue;
                    scale_t := cpx(one_r / temp, zero_r);
                    for jr:1..n fset(vl, jr, jc, fget(vl.*, jr, jc)*scale_t);
                }
            }
            if ilvr {
                #if RT == float32 { cggbak(#char "P", #char "R", n, ilo, ihi, lscale, rscale, n, vr, ldvr, *ierr); } else { zggbak(#char "P", #char "R", n, ilo, ihi, lscale, rscale, n, vr, ldvr, *ierr); }
                for jc:1..n {
                    temp = zero_r;
                    for jr:1..n temp = max(temp, _ggev_abs1(fget(vr.*, jr, jc), RT, MT));
                    if temp < smlnum then continue;
                    scale_t := cpx(one_r / temp, zero_r);
                    for jr:1..n fset(vr, jr, jc, fget(vr.*, jr, jc)*scale_t);
                }
            }
        }
    }

    if ilascl _ggev_scale_complex_segment(alpha, 1, n, anrm / anrmto, RT, MT);
    if ilbscl _ggev_scale_complex_segment(beta, 1, n, bnrm / bnrmto, RT, MT);
    fset(work, 1, _ggev_complex_from_real(_ggev_roundup_lwork(lwkopt, RT), RT, MT));
}
