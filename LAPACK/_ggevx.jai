// ? LAPACK *ggevx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggevx/cggevx/dggevx/zggevx
// *GGEVX computes the generalized eigenvalues, optionally the left
// and/or right generalized eigenvectors, balancing data, and
// reciprocal condition numbers.

sggevx :: (
    $$balanc: u8, $$jobvl: u8, $$jobvr: u8, $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D, alphai: *$ALPHAI/ArrayType1D, beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D, $$ldvl: int,
    vr: *$VR/ArrayType2D, $$ldvr: int,
    ilo: *int, ihi: *int,
    lscale: *$LSCALE/ArrayType1D, rscale: *$RSCALE/ArrayType1D,
    abnrm: *float32, bbnrm: *float32,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _ggevx_real(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info, RT=float32);
}

dggevx :: (
    $$balanc: u8, $$jobvl: u8, $$jobvr: u8, $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D, alphai: *$ALPHAI/ArrayType1D, beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D, $$ldvl: int,
    vr: *$VR/ArrayType2D, $$ldvr: int,
    ilo: *int, ihi: *int,
    lscale: *$LSCALE/ArrayType1D, rscale: *$RSCALE/ArrayType1D,
    abnrm: *float64, bbnrm: *float64,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _ggevx_real(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info, RT=float64);
}

cggevx :: (
    $$balanc: u8, $$jobvl: u8, $$jobvr: u8, $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D, $$ldvl: int,
    vr: *$VR/ArrayType2D, $$ldvr: int,
    ilo: *int, ihi: *int,
    lscale: *$LSCALE/ArrayType1D, rscale: *$RSCALE/ArrayType1D,
    abnrm: *float32, bbnrm: *float32,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _ggevx_complex(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork, info, RT=float32, MT=C32);
}

zggevx :: (
    $$balanc: u8, $$jobvl: u8, $$jobvr: u8, $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D, $$ldvl: int,
    vr: *$VR/ArrayType2D, $$ldvr: int,
    ilo: *int, ihi: *int,
    lscale: *$LSCALE/ArrayType1D, rscale: *$RSCALE/ArrayType1D,
    abnrm: *float64, bbnrm: *float64,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _ggevx_complex(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork, info, RT=float64, MT=C64);
}

#scope_file
_ggevx_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_ggevx_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_ggevx_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _ggevx_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_ggevx_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(lwork); } else { return droundup_lwork(lwork); }
}

_ggevx_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 { return ncast(MT, x); } else { return cpx(x, zero(RT)); }
}

_ggevx_lange :: (
    $$norm: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           WORK.data_type == RT;
} {
    #if MT == float32 { return slange(norm, m, n, a, lda, work); }
    else if MT == float64 { return dlange(norm, m, n, a, lda, work); }
    else if MT == C32 { return clange(norm, m, n, a, lda, work); }
    else { return zlange(norm, m, n, a, lda, work); }
}

_ggevx_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == float32 { slacpy(uplo, m, n, a, lda, b, ldb); }
    else if MT == float64 { dlacpy(uplo, m, n, a, lda, b, ldb); }
    else if MT == C32 { clacpy(uplo, m, n, a, lda, b, ldb); }
    else { zlacpy(uplo, m, n, a, lda, b, ldb); }
}

_ggevx_laset :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $MT, beta: MT, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == float32 { slaset(uplo, m, n, alpha, beta, a, lda); }
    else if MT == float64 { dlaset(uplo, m, n, alpha, beta, a, lda); }
    else if MT == C32 { claset(uplo, m, n, alpha, beta, a, lda); }
    else { zlaset(uplo, m, n, alpha, beta, a, lda); }
}

_ggevx_lascl_matrix :: (
    $$type: u8, $$kl: int, $$ku: int, cfrom: $RT, cto: RT, $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT;
} {
    #if MT == float32 { slascl(type, kl, ku, cfrom, cto, m, n, a, lda, info); }
    else if MT == float64 { dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info); }
    else if MT == C32 { clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info); }
    else { zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info); }
}

_ggevx_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "P"); } else { return dlamch(#char "P"); }
}

_ggevx_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "E"); } else { return dlamch(#char "E"); }
}

_ggevx_lamch_base :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "B"); } else { return dlamch(#char "B"); }
}

_ggevx_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "S"); } else { return dlamch(#char "S"); }
}

_ggevx_scale_real_segment :: (
    x: *$X/ArrayType1D, $$offset: int, $$count: int, scale: $RT, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    if count <= 0 then return;
    for i:offset..offset + count - 1 {
        fset(x, i, fget(x.*, i)*scale);
    }
}

_ggevx_scale_complex_segment :: (
    x: *$X/ArrayType1D, $$offset: int, $$count: int, scale: $RT, $RT: Type, $MT: Type
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

_ggevx_abs1 :: (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return abs(real(x)) + abs(imag(x));
}

_ggevx_real :: (
    $$balanc: u8, $$jobvl: u8, $$jobvr: u8, $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D, alphai: *$ALPHAI/ArrayType1D, beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D, $$ldvl: int, vr: *$VR/ArrayType2D, $$ldvr: int,
    ilo: *int, ihi: *int,
    lscale: *$LSCALE/ArrayType1D, rscale: *$RSCALE/ArrayType1D,
    abnrm: *RT, bbnrm: *RT,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    iwork: *$IWORK/ArrayType1D, bwork: *$BWORK/ArrayType1D,
    info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && B.data_type == RT &&
           ALPHAR.data_type == RT && ALPHAI.data_type == RT && BETA.data_type == RT &&
           VL.data_type == RT && VR.data_type == RT &&
           LSCALE.data_type == RT && RSCALE.data_type == RT &&
           RCONDE.data_type == RT && RCONDV.data_type == RT &&
           WORK.data_type == RT && IWORK.data_type == int && BWORK.data_type == bool;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if RT == float32 { "SGGEVX" } else { "DGGEVX" };
    geqrf_name :: #run if RT == float32 { "SGEQRF" } else { "DGEQRF" };
    ormqr_name :: #run if RT == float32 { "SORMQR" } else { "DORMQR" };
    orgqr_name :: #run if RT == float32 { "SORGQR" } else { "DORGQR" };

    ilascl, ilbscl, ilv, ilvl, ilvr, lquery, noscl, pair, wantsb, wantse, wantsn, wantsv : bool = ---;
    chtemp : u8 = ---;
    anrm, anrmto, bignum, bnrm, bnrmto, eps, smlnum, temp : RT = ---;
    i, icols, ierr, ijobvl, ijobvr, irows, itau, iwrk, iwrk1, j, jc, jr, m, maxwrk, minwrk, mm : int = ---;
    ldumma: [1] bool;
    query: [1] RT;

    info.* = 0;
    lquery = lwork == -1;
    if lsame(jobvl, #char "N") { ijobvl = 1; ilvl = false; } else if lsame(jobvl, #char "V") { ijobvl = 2; ilvl = true; } else { ijobvl = -1; ilvl = false; }
    if lsame(jobvr, #char "N") { ijobvr = 1; ilvr = false; } else if lsame(jobvr, #char "V") { ijobvr = 2; ilvr = true; } else { ijobvr = -1; ilvr = false; }
    ilv = ilvl || ilvr;
    noscl = lsame(balanc, #char "N") || lsame(balanc, #char "P");
    wantsn = lsame(sense, #char "N");
    wantse = lsame(sense, #char "E");
    wantsv = lsame(sense, #char "V");
    wantsb = lsame(sense, #char "B");

    if !(noscl || lsame(balanc, #char "S") || lsame(balanc, #char "B")) { info.* = -1; }
    else if ijobvl <= 0 { info.* = -2; }
    else if ijobvr <= 0 { info.* = -3; }
    else if !(wantsn || wantse || wantsv || wantsb) { info.* = -4; }
    else if n < 0 { info.* = -5; }
    else if lda < max(1, n) { info.* = -7; }
    else if ldb < max(1, n) { info.* = -9; }
    else if ldvl < 1 || (ilvl && ldvl < n) { info.* = -14; }
    else if ldvr < 1 || (ilvr && ldvr < n) { info.* = -16; }

    if info.* == 0 {
        if n == 0 {
            minwrk = 1; maxwrk = 1;
        } else {
            if noscl && !ilv { minwrk = 2*n; } else { minwrk = 6*n; }
            if wantse { minwrk = 10*n; }
            else if wantsv || wantsb { minwrk = 2*n*(n + 4) + 16; }
            maxwrk = minwrk;
            maxwrk = max(maxwrk, n + n*ilaenv(1, geqrf_name, " "s, n, 1, n, 0));
            maxwrk = max(maxwrk, n + n*ilaenv(1, ormqr_name, " "s, n, 1, n, 0));
            if ilvl {
                maxwrk = max(maxwrk, n + n*ilaenv(1, orgqr_name, " "s, n, 1, n, 0));
            }
        }
        fset(work, 1, _ggevx_roundup_lwork(maxwrk, RT));
            if lwork < minwrk && !lquery {
                info.* = -26;
            }
    }

    if info.* != 0 { xerbla(routine_name, -info.*); return; }
    if lquery then return;
    if n == 0 { abnrm.* = zero_r; bbnrm.* = zero_r; return; }

    eps = _ggevx_lamch_precision(RT);
    smlnum = _ggevx_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    anrm = _ggevx_lange(#char "M", n, n, a.*, lda, work, RT, RT);
    ilascl = false;
    if anrm > zero_r && anrm < smlnum { anrmto = smlnum; ilascl = true; } else if anrm > bignum { anrmto = bignum; ilascl = true; }
    if ilascl {
        _ggevx_lascl_matrix(#char "G", 0, 0, anrm, anrmto, n, n, a, lda, *ierr, RT, RT);
    }

    bnrm = _ggevx_lange(#char "M", n, n, b.*, ldb, work, RT, RT);
    ilbscl = false;
    if bnrm > zero_r && bnrm < smlnum { bnrmto = smlnum; ilbscl = true; } else if bnrm > bignum { bnrmto = bignum; ilbscl = true; }
    if ilbscl {
        _ggevx_lascl_matrix(#char "G", 0, 0, bnrm, bnrmto, n, n, b, ldb, *ierr, RT, RT);
    }

    #if RT == float32 { sggbal(balanc, n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, work, *ierr); }
    else { dggbal(balanc, n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, work, *ierr); }

    abnrm.* = _ggevx_lange(#char "1", n, n, a.*, lda, work, RT, RT);
    if ilascl {
        abnrm.* = abnrm.* * anrm / anrmto;
    }
    bbnrm.* = _ggevx_lange(#char "1", n, n, b.*, ldb, work, RT, RT);
    if ilbscl {
        bbnrm.* = bbnrm.* * bnrm / bnrmto;
    }

    irows = ihi.* + 1 - ilo.*;
    if ilv || !wantsn { icols = n + 1 - ilo.*; } else { icols = irows; }
    itau = 1;
    iwrk = itau + irows;
    tau := _ggevx_subvector(work.*, itau, max(0, irows));
    wrk := _ggevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    b_raw := _ggevx_matrix_storage(b.*, ldb, ilo.*, ilo.*, irows, icols);
    b_sub := array2d_view(b_raw, ldb, row_major=false);
    #if RT == float32 { sgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr); }
    else { dgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr); }

    a_raw := _ggevx_matrix_storage(a.*, lda, ilo.*, ilo.*, irows, icols);
    a_sub := array2d_view(a_raw, lda, row_major=false);
    #if RT == float32 { sormqr(#char "L", #char "T", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr); }
    else { dormqr(#char "L", #char "T", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr); }

    if ilvl {
        _ggevx_laset(#char "F", n, n, zero_r, one_r, vl, ldvl, RT);
        if irows > 1 {
            b_lower_raw := _ggevx_matrix_storage(b.*, ldb, ilo.* + 1, ilo.*, irows - 1, irows - 1);
            vl_lower_raw := _ggevx_matrix_storage(vl.*, ldvl, ilo.* + 1, ilo.*, irows - 1, irows - 1);
            b_lower := array2d_view(b_lower_raw, ldb, row_major=false);
            vl_lower := array2d_view(vl_lower_raw, ldvl, row_major=false);
            _ggevx_lacpy(#char "L", irows - 1, irows - 1, *b_lower, ldb, vl_lower, ldvl, RT);
        }
        vl_raw := _ggevx_matrix_storage(vl.*, ldvl, ilo.*, ilo.*, irows, irows);
        vl_sub := array2d_view(vl_raw, ldvl, row_major=false);
        #if RT == float32 { sorgqr(irows, irows, irows, *vl_sub, ldvl, tau, wrk, lwork - iwrk + 1, *ierr); }
        else { dorgqr(irows, irows, irows, *vl_sub, ldvl, tau, wrk, lwork - iwrk + 1, *ierr); }
    }
    if ilvr {
        _ggevx_laset(#char "F", n, n, zero_r, one_r, vr, ldvr, RT);
    }

    if ilv || !wantsn {
        #if RT == float32 { sgghrd(jobvl, jobvr, n, ilo.*, ihi.*, a, lda, b, ldb, vl, ldvl, vr, ldvr, *ierr); }
        else { dgghrd(jobvl, jobvr, n, ilo.*, ihi.*, a, lda, b, ldb, vl, ldvl, vr, ldvr, *ierr); }
    } else {
        sub_a_raw := _ggevx_matrix_storage(a.*, lda, ilo.*, ilo.*, irows, irows);
        sub_b_raw := _ggevx_matrix_storage(b.*, ldb, ilo.*, ilo.*, irows, irows);
        sub_a := array2d_view(sub_a_raw, lda, row_major=false);
        sub_b := array2d_view(sub_b_raw, ldb, row_major=false);
        #if RT == float32 { sgghrd(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, *ierr); }
        else { dgghrd(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, *ierr); }
    }

    chtemp = ifx ilv || !wantsn then #char "S" else #char "E";
    #if RT == float32 { shgeqz(chtemp, jobvl, jobvr, n, ilo.*, ihi.*, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, *ierr); }
    else { dhgeqz(chtemp, jobvl, jobvr, n, ilo.*, ihi.*, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, *ierr); }
    if ierr != 0 {
        if ierr > 0 && ierr <= n { info.* = ierr; } else if ierr > n && ierr <= 2*n { info.* = ierr - n; } else { info.* = n + 1; }
        fset(work, 1, _ggevx_roundup_lwork(maxwrk, RT));
        return;
    }

    if ilv || !wantsn {
        if ilv {
            if ilvl { if ilvr then chtemp = #char "B" else chtemp = #char "L"; } else { chtemp = #char "R"; }
            #if RT == float32 { stgevc(chtemp, #char "B", ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, *in, work, *ierr); }
            else { dtgevc(chtemp, #char "B", ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, *in, work, *ierr); }
            if ierr != 0 { info.* = n + 2; }
        }

        if info.* == 0 && !wantsn {
            pair = false;
            for i:1..n {
                if pair { pair = false; continue; }
                mm = 1;
                if i < n && fget(a.*, i + 1, i) != zero_r {
                    pair = true;
                    mm = 2;
                }
                for j:1..n fset(bwork, j, false);
                fset(bwork, i, true);
                if mm == 2 fset(bwork, i + 1, true);
                iwrk = mm*n + 1;
                iwrk1 = iwrk + mm*n;

                if wantse || wantsb {
                    left_sel := _ggevx_raw_slice(work.data, mm*n);
                    right_sel := _ggevx_raw_slice(work.data + (iwrk - 1), mm*n);
                    left_mat := array2d_view(left_sel, n, row_major=false);
                    right_mat := array2d_view(right_sel, n, row_major=false);
                    work_sel := _ggevx_subvector(work.*, iwrk1, max(0, lwork - iwrk1 + 1));
                    #if RT == float32 { stgevc(#char "B", #char "S", bwork, n, a, lda, b, ldb, *left_mat, n, *right_mat, n, mm, *m, work_sel, *ierr); }
                    else { dtgevc(#char "B", #char "S", bwork, n, a, lda, b, ldb, *left_mat, n, *right_mat, n, mm, *m, work_sel, *ierr); }
                    if ierr != 0 { info.* = n + 2; break; }
                }

                left_sel := _ggevx_raw_slice(work.data, mm*n);
                right_sel := _ggevx_raw_slice(work.data + (iwrk - 1), mm*n);
                left_mat := array2d_view(left_sel, n, row_major=false);
                right_mat := array2d_view(right_sel, n, row_major=false);
                cond_work := _ggevx_subvector(work.*, iwrk1, max(0, lwork - iwrk1 + 1));
                rconde_i := _ggevx_subvector(rconde.*, i, max(0, n - i + 1));
                rcondv_i := _ggevx_subvector(rcondv.*, i, max(0, n - i + 1));
                #if RT == float32 { stgsna(sense, #char "S", bwork, n, a, lda, b, ldb, *left_mat, n, *right_mat, n, rconde_i, rcondv_i, mm, *m, cond_work, lwork - iwrk1 + 1, iwork, *ierr); }
                else { dtgsna(sense, #char "S", bwork, n, a, lda, b, ldb, *left_mat, n, *right_mat, n, rconde_i, rcondv_i, mm, *m, cond_work, lwork - iwrk1 + 1, iwork, *ierr); }
            }
        }
    }

    if ilvl {
        #if RT == float32 { sggbak(balanc, #char "L", n, ilo.*, ihi.*, lscale, rscale, n, vl, ldvl, *ierr); }
        else { dggbak(balanc, #char "L", n, ilo.*, ihi.*, lscale, rscale, n, vl, ldvl, *ierr); }
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
        #if RT == float32 { sggbak(balanc, #char "R", n, ilo.*, ihi.*, lscale, rscale, n, vr, ldvr, *ierr); }
        else { dggbak(balanc, #char "R", n, ilo.*, ihi.*, lscale, rscale, n, vr, ldvr, *ierr); }
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

    if ilascl { _ggevx_scale_real_segment(alphar, 1, n, anrm / anrmto, RT); _ggevx_scale_real_segment(alphai, 1, n, anrm / anrmto, RT); }
    if ilbscl {
        _ggevx_scale_real_segment(beta, 1, n, bnrm / bnrmto, RT);
    }
    fset(work, 1, _ggevx_roundup_lwork(maxwrk, RT));
}

_ggevx_complex :: (
    $$balanc: u8, $$jobvl: u8, $$jobvr: u8, $$sense: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    alpha: *$ALPHA/ArrayType1D, beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D, $$ldvl: int, vr: *$VR/ArrayType2D, $$ldvr: int,
    ilo: *int, ihi: *int,
    lscale: *$LSCALE/ArrayType1D, rscale: *$RSCALE/ArrayType1D,
    abnrm: *RT, bbnrm: *RT,
    rconde: *$RCONDE/ArrayType1D, rcondv: *$RCONDV/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, bwork: *$BWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && B.data_type == MT &&
           ALPHA.data_type == MT && BETA.data_type == MT &&
           VL.data_type == MT && VR.data_type == MT &&
           LSCALE.data_type == RT && RSCALE.data_type == RT &&
           RCONDE.data_type == RT && RCONDV.data_type == RT &&
           WORK.data_type == MT && RWORK.data_type == RT &&
           IWORK.data_type == int && BWORK.data_type == bool;
} {
    one_r :: #run one(RT);
    zero_r :: #run zero(RT);
    zero_t :: #run cpx(zero(RT), zero(RT));
    one_t :: #run cpx(one(RT), zero(RT));
    routine_name :: #run if RT == float32 { "CGGEVX" } else { "ZGGEVX" };
    geqrf_name :: #run if RT == float32 { "CGEQRF" } else { "ZGEQRF" };
    unmqr_name :: #run if RT == float32 { "CUNMQR" } else { "ZUNMQR" };
    ungqr_name :: #run if RT == float32 { "CUNGQR" } else { "ZUNGQR" };

    ilascl, ilbscl, ilv, ilvl, ilvr, lquery, noscl, wantsb, wantse, wantsn, wantsv : bool = ---;
    chtemp : u8 = ---;
    anrm, anrmto, bignum, bnrm, bnrmto, eps, smlnum, temp : RT = ---;
    ierr, ihi0, ijobvl, ijobvr, ileft, ilo0, in, iright, irwrk, irows, icols, itau, iwrk, iwrk1, i, j, jc, jr, m, maxwrk, minwrk : int = ---;
    ldumma: [1] bool;
    query: [1] MT;

    info.* = 0;
    lquery = lwork == -1;
    if lsame(jobvl, #char "N") { ijobvl = 1; ilvl = false; } else if lsame(jobvl, #char "V") { ijobvl = 2; ilvl = true; } else { ijobvl = -1; ilvl = false; }
    if lsame(jobvr, #char "N") { ijobvr = 1; ilvr = false; } else if lsame(jobvr, #char "V") { ijobvr = 2; ilvr = true; } else { ijobvr = -1; ilvr = false; }
    ilv = ilvl || ilvr;
    noscl = lsame(balanc, #char "N") || lsame(balanc, #char "P");
    wantsn = lsame(sense, #char "N");
    wantse = lsame(sense, #char "E");
    wantsv = lsame(sense, #char "V");
    wantsb = lsame(sense, #char "B");

    if !(noscl || lsame(balanc, #char "S") || lsame(balanc, #char "B")) { info.* = -1; }
    else if ijobvl <= 0 { info.* = -2; }
    else if ijobvr <= 0 { info.* = -3; }
    else if !(wantsn || wantse || wantsv || wantsb) { info.* = -4; }
    else if n < 0 { info.* = -5; }
    else if lda < max(1, n) { info.* = -7; }
    else if ldb < max(1, n) { info.* = -9; }
    else if ldvl < 1 || (ilvl && ldvl < n) { info.* = -13; }
    else if ldvr < 1 || (ilvr && ldvr < n) { info.* = -15; }

    if info.* == 0 {
        if n == 0 {
            minwrk = 1; maxwrk = 1;
        } else {
            minwrk = 2*n;
            if wantse { minwrk = 4*n; }
            else if wantsv || wantsb { minwrk = 2*n*(n + 1); }
            maxwrk = minwrk;
            maxwrk = max(maxwrk, n + n*ilaenv(1, geqrf_name, " "s, n, 1, n, 0));
            maxwrk = max(maxwrk, n + n*ilaenv(1, unmqr_name, " "s, n, 1, n, 0));
            if ilvl {
                maxwrk = max(maxwrk, n + n*ilaenv(1, ungqr_name, " "s, n, 1, n, 0));
            }
        }
        fset(work, 1, ifx n == 0 then one_t else _ggevx_complex_from_real(_ggevx_roundup_lwork(maxwrk, RT), RT, MT));
        if lwork < minwrk && !lquery {
            info.* = -25;
        }
    }

    if info.* != 0 { xerbla(routine_name, -info.*); return; }
    if lquery then return;
    if n == 0 { abnrm.* = zero_r; bbnrm.* = zero_r; return; }

    eps = _ggevx_lamch_precision(RT);
    smlnum = _ggevx_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = one_r / smlnum;

    anrm = _ggevx_lange(#char "M", n, n, a.*, lda, rwork, RT, MT);
    ilascl = false;
    if anrm > zero_r && anrm < smlnum { anrmto = smlnum; ilascl = true; } else if anrm > bignum { anrmto = bignum; ilascl = true; }
    if ilascl {
        _ggevx_lascl_matrix(#char "G", 0, 0, anrm, anrmto, n, n, a, lda, *ierr, RT, MT);
    }

    bnrm = _ggevx_lange(#char "M", n, n, b.*, ldb, rwork, RT, MT);
    ilbscl = false;
    if bnrm > zero_r && bnrm < smlnum { bnrmto = smlnum; ilbscl = true; } else if bnrm > bignum { bnrmto = bignum; ilbscl = true; }
    if ilbscl {
        _ggevx_lascl_matrix(#char "G", 0, 0, bnrm, bnrmto, n, n, b, ldb, *ierr, RT, MT);
    }

    #if RT == float32 { cggbal(balanc, n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, rwork, *ierr); }
    else { zggbal(balanc, n, a, lda, b, ldb, *ilo, *ihi, lscale, rscale, rwork, *ierr); }

    abnrm.* = _ggevx_lange(#char "1", n, n, a.*, lda, rwork, RT, MT);
    if ilascl {
        abnrm.* = abnrm.* * anrm / anrmto;
    }
    bbnrm.* = _ggevx_lange(#char "1", n, n, b.*, ldb, rwork, RT, MT);
    if ilbscl {
        bbnrm.* = bbnrm.* * bnrm / bnrmto;
    }

    irows = ihi.* + 1 - ilo.*;
    if ilv || !wantsn { icols = n + 1 - ilo.*; } else { icols = irows; }
    itau = 1; iwrk = itau + irows;
    tau := _ggevx_subvector(work.*, itau, max(0, irows));
    wrk := _ggevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    rwrk := _ggevx_subvector(rwork.*, 1, max(0, rwork.count));
    b_raw := _ggevx_matrix_storage(b.*, ldb, ilo.*, ilo.*, irows, icols);
    b_sub := array2d_view(b_raw, ldb, row_major=false);
    #if RT == float32 { cgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr); }
    else { zgeqrf(irows, icols, *b_sub, ldb, tau, wrk, lwork - iwrk + 1, *ierr); }
    a_raw := _ggevx_matrix_storage(a.*, lda, ilo.*, ilo.*, irows, icols);
    a_sub := array2d_view(a_raw, lda, row_major=false);
    #if RT == float32 { cunmqr(#char "L", #char "C", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr); }
    else { zunmqr(#char "L", #char "C", irows, icols, irows, *b_sub, ldb, tau, *a_sub, lda, wrk, lwork - iwrk + 1, *ierr); }

    if ilvl {
        _ggevx_laset(#char "F", n, n, zero_t, one_t, vl, ldvl, MT);
        if irows > 1 {
            b_lower_raw := _ggevx_matrix_storage(b.*, ldb, ilo.* + 1, ilo.*, irows - 1, irows - 1);
            vl_lower_raw := _ggevx_matrix_storage(vl.*, ldvl, ilo.* + 1, ilo.*, irows - 1, irows - 1);
            b_lower := array2d_view(b_lower_raw, ldb, row_major=false);
            vl_lower := array2d_view(vl_lower_raw, ldvl, row_major=false);
            _ggevx_lacpy(#char "L", irows - 1, irows - 1, *b_lower, ldb, vl_lower, ldvl, MT);
        }
        vl_raw := _ggevx_matrix_storage(vl.*, ldvl, ilo.*, ilo.*, irows, irows);
        vl_sub := array2d_view(vl_raw, ldvl, row_major=false);
        #if RT == float32 { cungqr(irows, irows, irows, *vl_sub, ldvl, tau, wrk, lwork - iwrk + 1, *ierr); }
        else { zungqr(irows, irows, irows, *vl_sub, ldvl, tau, wrk, lwork - iwrk + 1, *ierr); }
    }
    if ilvr {
        _ggevx_laset(#char "F", n, n, zero_t, one_t, vr, ldvr, MT);
    }

    if ilv || !wantsn {
        #if RT == float32 { cgghrd(jobvl, jobvr, n, ilo.*, ihi.*, a, lda, b, ldb, vl, ldvl, vr, ldvr, *ierr); }
        else { zgghrd(jobvl, jobvr, n, ilo.*, ihi.*, a, lda, b, ldb, vl, ldvl, vr, ldvr, *ierr); }
    } else {
        sub_a_raw := _ggevx_matrix_storage(a.*, lda, ilo.*, ilo.*, irows, irows);
        sub_b_raw := _ggevx_matrix_storage(b.*, ldb, ilo.*, ilo.*, irows, irows);
        sub_a := array2d_view(sub_a_raw, lda, row_major=false);
        sub_b := array2d_view(sub_b_raw, ldb, row_major=false);
        #if RT == float32 { cgghrd(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, *ierr); }
        else { zgghrd(#char "N", #char "N", irows, 1, irows, *sub_a, lda, *sub_b, ldb, vl, ldvl, vr, ldvr, *ierr); }
    }

    iwrk = itau;
    wrk = _ggevx_subvector(work.*, iwrk, max(0, lwork - iwrk + 1));
    chtemp = ifx ilv || !wantsn then #char "S" else #char "E";
    #if RT == float32 { chgeqz(chtemp, jobvl, jobvr, n, ilo.*, ihi.*, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, rwork, *ierr); }
    else { zhgeqz(chtemp, jobvl, jobvr, n, ilo.*, ihi.*, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, wrk, lwork - iwrk + 1, rwork, *ierr); }
    if ierr != 0 {
        if ierr > 0 && ierr <= n { info.* = ierr; } else if ierr > n && ierr <= 2*n { info.* = ierr - n; } else { info.* = n + 1; }
        fset(work, 1, _ggevx_complex_from_real(_ggevx_roundup_lwork(maxwrk, RT), RT, MT));
        return;
    }

    if ilv || !wantsn {
        if ilv {
            if ilvl { if ilvr then chtemp = #char "B" else chtemp = #char "L"; } else { chtemp = #char "R"; }
            #if RT == float32 { ctgevc(chtemp, #char "B", ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, *in, wrk, rwork, *ierr); }
            else { ztgevc(chtemp, #char "B", ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, *in, wrk, rwork, *ierr); }
            if ierr != 0 { info.* = n + 2; }
        }

        if info.* == 0 && !wantsn {
            for i:1..n {
                for j:1..n fset(bwork, j, false);
                fset(bwork, i, true);
                iwrk = n + 1;
                iwrk1 = iwrk + n;

                if wantse || wantsb {
                    left_sel := _ggevx_raw_slice(work.data, n);
                    right_sel := _ggevx_raw_slice(work.data + (iwrk - 1), n);
                    left_mat := array2d_view(left_sel, n, row_major=false);
                    right_mat := array2d_view(right_sel, n, row_major=false);
                    sel_work := _ggevx_subvector(work.*, iwrk1, max(0, lwork - iwrk1 + 1));
                    #if RT == float32 { ctgevc(#char "B", #char "S", bwork, n, a, lda, b, ldb, *left_mat, n, *right_mat, n, 1, *m, sel_work, rwork, *ierr); }
                    else { ztgevc(#char "B", #char "S", bwork, n, a, lda, b, ldb, *left_mat, n, *right_mat, n, 1, *m, sel_work, rwork, *ierr); }
                    if ierr != 0 { info.* = n + 2; break; }
                }

                left_sel := _ggevx_raw_slice(work.data, n);
                right_sel := _ggevx_raw_slice(work.data + (iwrk - 1), n);
                left_mat := array2d_view(left_sel, n, row_major=false);
                right_mat := array2d_view(right_sel, n, row_major=false);
                cond_work := _ggevx_subvector(work.*, iwrk1, max(0, lwork - iwrk1 + 1));
                rconde_i := _ggevx_subvector(rconde.*, i, max(0, n - i + 1));
                rcondv_i := _ggevx_subvector(rcondv.*, i, max(0, n - i + 1));
                #if RT == float32 { ctgsna(sense, #char "S", bwork, n, a, lda, b, ldb, *left_mat, n, *right_mat, n, rconde_i, rcondv_i, 1, *m, cond_work, lwork - iwrk1 + 1, iwork, *ierr); }
                else { ztgsna(sense, #char "S", bwork, n, a, lda, b, ldb, *left_mat, n, *right_mat, n, rconde_i, rcondv_i, 1, *m, cond_work, lwork - iwrk1 + 1, iwork, *ierr); }
            }
        }
    }

    if ilvl {
        #if RT == float32 { cggbak(balanc, #char "L", n, ilo.*, ihi.*, lscale, rscale, n, vl, ldvl, *ierr); }
        else { zggbak(balanc, #char "L", n, ilo.*, ihi.*, lscale, rscale, n, vl, ldvl, *ierr); }
        for jc:1..n {
            temp = zero_r;
            for jr:1..n temp = max(temp, _ggevx_abs1(fget(vl.*, jr, jc), RT, MT));
            if temp < smlnum then continue;
            scale_t := cpx(one_r / temp, zero_r);
            for jr:1..n fset(vl, jr, jc, fget(vl.*, jr, jc)*scale_t);
        }
    }

    if ilvr {
        #if RT == float32 { cggbak(balanc, #char "R", n, ilo.*, ihi.*, lscale, rscale, n, vr, ldvr, *ierr); }
        else { zggbak(balanc, #char "R", n, ilo.*, ihi.*, lscale, rscale, n, vr, ldvr, *ierr); }
        for jc:1..n {
            temp = zero_r;
            for jr:1..n temp = max(temp, _ggevx_abs1(fget(vr.*, jr, jc), RT, MT));
            if temp < smlnum then continue;
            scale_t := cpx(one_r / temp, zero_r);
            for jr:1..n fset(vr, jr, jc, fget(vr.*, jr, jc)*scale_t);
        }
    }

    if ilascl {
        _ggevx_scale_complex_segment(alpha, 1, n, anrm / anrmto, RT, MT);
    }
    if ilbscl {
        _ggevx_scale_complex_segment(beta, 1, n, bnrm / bnrmto, RT, MT);
    }
    fset(work, 1, _ggevx_complex_from_real(_ggevx_roundup_lwork(maxwrk, RT), RT, MT));
}
