// ? LAPACK *gbsvx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbsvx/cgbsvx/dgbsvx/zgbsvx
// *GBSVX uses the LU factorization to compute the solution to a
// system of linear equations for GB matrices, and provides error
// bounds plus a condition estimate.

sgbsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    afb: *$AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: *$IPIV/ArrayType1D,
    equed: *u8,
    r: *$R/ArrayType1D,
    c: *$C/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *float32,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, RT=float32, MT=float32);
}

dgbsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    afb: *$AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: *$IPIV/ArrayType1D,
    equed: *u8,
    r: *$R/ArrayType1D,
    c: *$C/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *float64,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, RT=float64, MT=float64);
}

cgbsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    afb: *$AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: *$IPIV/ArrayType1D,
    equed: *u8,
    r: *$R/ArrayType1D,
    c: *$C/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *float32,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, RT=float32, MT=C32);
}

zgbsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    afb: *$AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: *$IPIV/ArrayType1D,
    equed: *u8,
    r: *$R/ArrayType1D,
    c: *$C/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *float64,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gbsvx_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gbsvx_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_gbsvx_abs_for_growth :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    return abs(z);
}

_gbsvx_copy_band_to_afb :: (
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    afb: *$AFB/ArrayType2D
) #modify {
    return AB.data_type == AFB.data_type;
} {
    j1, j2 : int = ---;
    for j:1..n {
        j1 = max(j-ku, 1);
        j2 = min(j+kl, n);
        for i:j1..j2 {
            fset(afb, kl+ku+1-j+i, j, fget(ab, ku+1-j+i, j));
        }
    }
}

_gbsvx_rpvgrw_partial :: (
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$ncols: int,
    ab: $AB/ArrayType2D,
    afb: $AFB/ArrayType2D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           AFB.data_type == MT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    anorm, umax : RT = zero_r, zero_r;
    kband, row : int = ---;

    if ncols <= 0 {
        return one_r;
    }

    for j:1..ncols {
        for i:max(ku+2-j, 1)..min(n+ku+1-j, kl+ku+1) {
            anorm = max(anorm, _gbsvx_abs_for_growth(fget(ab, i, j), RT, MT));
        }
    }

    kband = min(ncols-1, kl+ku);
    for j:1..ncols {
        for i:max(1, j-kband)..j {
            row = kl + ku + 1 + i - j;
            umax = max(umax, _gbsvx_abs_for_growth(fget(afb, row, j), RT, MT));
        }
    }

    if umax == zero_r {
        return one_r;
    }
    return anorm / umax;
}

_gbsvx_langb :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           WORK.data_type == RT;
} {
    #if MT == float32 {
        return slangb(norm, n, kl, ku, ab, ldab, work);
    } else if MT == float64 {
        return dlangb(norm, n, kl, ku, ab, ldab, work);
    } else if MT == C32 {
        return clangb(norm, n, kl, ku, ab, ldab, work);
    } else {
        return zlangb(norm, n, kl, ku, ab, ldab, work);
    }
}

_gbsvx_lacpy :: (
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

_gbsvx_laqgb :: (
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    r: $R/ArrayType1D,
    c: $C/ArrayType1D,
    rowcnd: $RT,
    colcnd: $RT,
    amax: $RT,
    equed: *u8,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           R.data_type == RT &&
           C.data_type == RT;
} {
    #if MT == float32 {
        slaqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
    } else if MT == float64 {
        dlaqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
    } else if MT == C32 {
        claqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
    } else {
        zlaqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
    }
}

_gbsvx_gbtrf :: (
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AB.data_type == MT &&
           IPIV.data_type == int;
} {
    #if MT == float32 {
        sgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    } else if MT == float64 {
        dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    } else if MT == C32 {
        cgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    } else {
        zgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    }
}

_gbsvx_gbtrs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    #if MT == float32 {
        sgbtrs(trans, n, kl, ku, nrhs, afb, ldafb, ipiv, b, ldb, info);
    } else if MT == float64 {
        dgbtrs(trans, n, kl, ku, nrhs, afb, ldafb, ipiv, b, ldb, info);
    } else if MT == C32 {
        cgbtrs(trans, n, kl, ku, nrhs, afb, ldafb, ipiv, b, ldb, info);
    } else {
        zgbtrs(trans, n, kl, ku, nrhs, afb, ldafb, ipiv, b, ldb, info);
    }
}

_gbsvx_gbequ :: (
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    r: *$R/ArrayType1D,
    c: *$C/ArrayType1D,
    rowcnd: *$RT,
    colcnd: *$RT,
    amax: *$RT,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           R.data_type == RT &&
           C.data_type == RT;
} {
    #if MT == float32 {
        sgbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    } else if MT == float64 {
        dgbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    } else if MT == C32 {
        cgbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    } else {
        zgbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    }
}

_gbsvx_gbcon :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    anorm: $RT,
    rcond: *$RT,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        sgbcon(norm, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, aux, info);
    } else if MT == float64 {
        dgbcon(norm, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, aux, info);
    } else if MT == C32 {
        cgbcon(norm, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, aux, info);
    } else {
        zgbcon(norm, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, aux, info);
    }
}

_gbsvx_gbrfs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT &&
           X.data_type == MT &&
           FERR.data_type == RT &&
           BERR.data_type == RT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        sgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    } else if MT == float64 {
        dgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    } else if MT == C32 {
        cgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    } else {
        zgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    }
}

_gbsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    afb: *$AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: *$IPIV/ArrayType1D,
    equed: *u8,
    r: *$R/ArrayType1D,
    c: *$C/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *$RT,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           R.data_type == RT &&
           C.data_type == RT &&
           B.data_type == MT &&
           X.data_type == MT &&
           FERR.data_type == RT &&
           BERR.data_type == RT &&
           ((MT == RT && WORK.data_type == RT && AUX.data_type == int) ||
            (MT == complex_with_inner_type(RT) && WORK.data_type == MT && AUX.data_type == RT));
} {
    #if CHECKS {
        #run {
            assert(RT == FERR.data_type);
            assert(RT == BERR.data_type);
        }
    }

    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == float32 {
        "SGBSVX"
    } else if MT == float64 {
        "DGBSVX"
    } else if MT == C32 {
        "CGBSVX"
    } else {
        "ZGBSVX"
    };

    amax, anorm, bignum, colcnd, rcmax, rcmin, rowcnd, rpvgrw, smlnum : RT = ---;
    colequ, equil, nofact, notran, rowequ : bool = ---;
    infequ, i, j : int = ---;
    norm : u8 = ---;

    info.* = 0;
    nofact = lsame(fact, #char "N");
    equil  = lsame(fact, #char "E");
    notran = lsame(trans, #char "N");
    if nofact || equil {
        equed.* = #char "N";
        rowequ = false;
        colequ = false;
    } else {
        rowequ = lsame(equed.*, #char "R") || lsame(equed.*, #char "B");
        colequ = lsame(equed.*, #char "C") || lsame(equed.*, #char "B");
        smlnum = _gbsvx_lamch_safe_minimum(RT);
        bignum = one_r / smlnum;
    }

    if !nofact && !equil && !lsame(fact, #char "F") {
        info.* = -1;
    } else if !notran && !lsame(trans, #char "T") && !lsame(trans, #char "C") {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if kl < 0 {
        info.* = -4;
    } else if ku < 0 {
        info.* = -5;
    } else if nrhs < 0 {
        info.* = -6;
    } else if ldab < kl + ku + 1 {
        info.* = -8;
    } else if ldafb < 2*kl + ku + 1 {
        info.* = -10;
    } else if lsame(fact, #char "F") && !(rowequ || colequ || lsame(equed.*, #char "N")) {
        info.* = -12;
    } else {
        if rowequ {
            rcmin = bignum;
            rcmax = zero_r;
            for j:1..n {
                rcmin = min(rcmin, fget(r.*, j));
                rcmax = max(rcmax, fget(r.*, j));
            }
            if rcmin <= zero_r {
                info.* = -13;
            } else if n > 0 {
                rowcnd = max(rcmin, smlnum) / min(rcmax, bignum);
            } else {
                rowcnd = one_r;
            }
        }
        if colequ && info.* == 0 {
            rcmin = bignum;
            rcmax = zero_r;
            for j:1..n {
                rcmin = min(rcmin, fget(c.*, j));
                rcmax = max(rcmax, fget(c.*, j));
            }
            if rcmin <= zero_r {
                info.* = -14;
            } else if n > 0 {
                colcnd = max(rcmin, smlnum) / min(rcmax, bignum);
            } else {
                colcnd = one_r;
            }
        }
        if info.* == 0 {
            if ldb < max(1, n) {
                info.* = -16;
            } else if ldx < max(1, n) {
                info.* = -18;
            }
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if equil {
        _gbsvx_gbequ(n, n, kl, ku, ab.*, ldab, r, c, *rowcnd, *colcnd, *amax, *infequ, RT, MT);
        if infequ == 0 {
            _gbsvx_laqgb(n, n, kl, ku, ab, ldab, r.*, c.*, rowcnd, colcnd, amax, equed, RT, MT);
            rowequ = lsame(equed.*, #char "R") || lsame(equed.*, #char "B");
            colequ = lsame(equed.*, #char "C") || lsame(equed.*, #char "B");
        }
    }

    if notran {
        if rowequ {
            for j:1..nrhs {
                for i:1..n {
                    fset(b, i, j, fget(r.*, i) * fget(b.*, i, j));
                }
            }
        }
    } else if colequ {
        for j:1..nrhs {
            for i:1..n {
                fset(b, i, j, fget(c.*, i) * fget(b.*, i, j));
            }
        }
    }

    if nofact || equil {
        _gbsvx_copy_band_to_afb(n, kl, ku, ab.*, afb);
        _gbsvx_gbtrf(n, n, kl, ku, afb, ldafb, ipiv, info, MT);

        if info.* > 0 {
            rpvgrw = _gbsvx_rpvgrw_partial(n, kl, ku, info.*, ab.*, afb.*, RT, MT);
            #if MT == float32 || MT == float64 {
                fset(work, 1, rpvgrw);
            } else {
                fset(aux, 1, rpvgrw);
            }
            rcond.* = zero_r;
            return;
        }
    }

    if notran {
        norm = #char "1";
    } else {
        norm = #char "I";
    }

    #if MT == float32 || MT == float64 {
        anorm = _gbsvx_langb(norm, n, kl, ku, ab.*, ldab, work, RT, MT);
    } else {
        anorm = _gbsvx_langb(norm, n, kl, ku, ab.*, ldab, aux, RT, MT);
    }
    rpvgrw = _gbsvx_rpvgrw_partial(n, kl, ku, n, ab.*, afb.*, RT, MT);

    _gbsvx_gbcon(norm, n, kl, ku, afb.*, ldafb, ipiv.*, anorm, rcond, work, aux, info, RT, MT);
    _gbsvx_lacpy(#char "F", n, nrhs, b.*, ldb, x, ldx, MT);
    _gbsvx_gbtrs(trans, n, kl, ku, nrhs, afb.*, ldafb, ipiv.*, x, ldx, info, MT);
    _gbsvx_gbrfs(trans, n, kl, ku, nrhs, ab.*, ldab, afb.*, ldafb, ipiv.*, b.*, ldb, x, ldx, ferr, berr, work, aux, info, RT, MT);

    if notran {
        if colequ {
            for j:1..nrhs {
                for i:1..n {
                    fset(x, i, j, fget(c.*, i) * fget(x.*, i, j));
                }
            }
            for j:1..nrhs {
                fset(ferr, j, fget(ferr.*, j) / colcnd);
            }
        }
    } else if rowequ {
        for j:1..nrhs {
            for i:1..n {
                fset(x, i, j, fget(r.*, i) * fget(x.*, i, j));
            }
        }
        for j:1..nrhs {
            fset(ferr, j, fget(ferr.*, j) / rowcnd);
        }
    }

    if rcond.* < _gbsvx_lamch_epsilon(RT) {
        info.* = n + 1;
    }

    #if MT == float32 || MT == float64 {
        fset(work, 1, rpvgrw);
    } else {
        fset(aux, 1, rpvgrw);
    }
}
