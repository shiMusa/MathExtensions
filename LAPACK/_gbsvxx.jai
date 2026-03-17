// ? LAPACK *gbsvxx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbsvxx/cgbsvxx/dgbsvxx/zgbsvxx
// *GBSVXX uses the LU factorization to compute the solution to a
// system of linear equations for GB matrices and can return
// normwise and componentwise error bounds.

sgbsvxx :: (
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
    rpvgrw: *float32,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info, RT=float32, MT=float32);
}

dgbsvxx :: (
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
    rpvgrw: *float64,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info, RT=float64, MT=float64);
}

cgbsvxx :: (
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
    rpvgrw: *float32,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float32, MT=C32);
}

zgbsvxx :: (
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
    rpvgrw: *float64,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gbsvxx_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gbsvxx_copy_ab_to_afb :: (
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    afb: *$AFB/ArrayType2D
) #modify {
    return AB.data_type == AFB.data_type;
} {
    for j:1..n {
        for i:kl+1..2*kl+ku+1 {
            fset(afb, i, j, fget(ab, i-kl, j));
        }
    }
}

_gbsvxx_lacpy :: (
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

_gbsvxx_laqgb :: (
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

_gbsvxx_gbequb :: (
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
        sgbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    } else if MT == float64 {
        dgbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    } else if MT == C32 {
        cgbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    } else {
        zgbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    }
}

_gbsvxx_gbtrf :: (
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

_gbsvxx_gbtrs :: (
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

_gbsvxx_lascl2 :: (
    $$m: int,
    $$n: int,
    d: $D/ArrayType1D,
    x: *$X/ArrayType2D,
    $$ldx: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           D.data_type == RT &&
           X.data_type == MT;
} {
    #if MT == float32 {
        slascl2(m, n, d, x, ldx);
    } else if MT == float64 {
        dlascl2(m, n, d, x, ldx);
    } else if MT == C32 {
        clascl2(m, n, d, x, ldx);
    } else {
        zlascl2(m, n, d, x, ldx);
    }
}

_gbsvxx_gbrpvgrw :: (
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$ncols: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           AFB.data_type == MT;
} {
    #if MT == float32 {
        return sla_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb);
    } else if MT == float64 {
        return dla_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb);
    } else if MT == C32 {
        return cla_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb);
    } else {
        return zla_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb);
    }
}

_gbsvxx_gbrfsx :: (
    $$trans: u8,
    $$equed: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    r: $R/ArrayType1D,
    c: $C/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *$RT,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
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
           BERR.data_type == RT &&
           ERR_BNDS_NORM.data_type == RT &&
           ERR_BNDS_COMP.data_type == RT &&
           PARAMS.data_type == RT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        sgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info);
    } else if MT == float64 {
        dgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info);
    } else if MT == C32 {
        cgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info);
    } else {
        zgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info);
    }
}

_gbsvxx :: (
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
    rpvgrw: *$RT,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
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
           BERR.data_type == RT &&
           ERR_BNDS_NORM.data_type == RT &&
           ERR_BNDS_COMP.data_type == RT &&
           PARAMS.data_type == RT &&
           ((MT == RT && WORK.data_type == RT && AUX.data_type == int) ||
            (MT == complex_with_inner_type(RT) && WORK.data_type == MT && AUX.data_type == RT));
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == float32 {
        "SGBSVXX"
    } else if MT == float64 {
        "DGBSVXX"
    } else if MT == C32 {
        "CGBSVXX"
    } else {
        "ZGBSVXX"
    };

    amax, bignum, colcnd, rcmax, rcmin, rowcnd, smlnum : RT = ---;
    colequ, equil, nofact, notran, rowequ : bool = ---;
    infequ, j : int = ---;

    info.* = 0;
    nofact = lsame(fact, #char "N");
    equil  = lsame(fact, #char "E");
    notran = lsame(trans, #char "N");
    smlnum = _gbsvxx_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    if nofact || equil {
        equed.* = #char "N";
        rowequ = false;
        colequ = false;
    } else {
        rowequ = lsame(equed.*, #char "R") || lsame(equed.*, #char "B");
        colequ = lsame(equed.*, #char "C") || lsame(equed.*, #char "B");
    }

    rpvgrw.* = zero_r;

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
                info.* = -15;
            } else if ldx < max(1, n) {
                info.* = -16;
            }
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if equil {
        _gbsvxx_gbequb(n, n, kl, ku, ab.*, ldab, r, c, *rowcnd, *colcnd, *amax, *infequ, RT, MT);
        if infequ == 0 {
            _gbsvxx_laqgb(n, n, kl, ku, ab, ldab, r.*, c.*, rowcnd, colcnd, amax, equed, RT, MT);
            rowequ = lsame(equed.*, #char "R") || lsame(equed.*, #char "B");
            colequ = lsame(equed.*, #char "C") || lsame(equed.*, #char "B");
        }

        if !rowequ {
            for j:1..n {
                fset(r, j, one_r);
            }
        }
        if !colequ {
            for j:1..n {
                fset(c, j, one_r);
            }
        }
    }

    if notran {
        if rowequ {
            _gbsvxx_lascl2(n, nrhs, r.*, b, ldb, RT, MT);
        }
    } else if colequ {
        _gbsvxx_lascl2(n, nrhs, c.*, b, ldb, RT, MT);
    }

    if nofact || equil {
        _gbsvxx_copy_ab_to_afb(n, kl, ku, ab.*, afb);
        _gbsvxx_gbtrf(n, n, kl, ku, afb, ldafb, ipiv, info, MT);

        if info.* > 0 {
            rpvgrw.* = _gbsvxx_gbrpvgrw(n, kl, ku, info.*, ab.*, ldab, afb.*, ldafb, RT, MT);
            return;
        }
    }

    rpvgrw.* = _gbsvxx_gbrpvgrw(n, kl, ku, n, ab.*, ldab, afb.*, ldafb, RT, MT);

    _gbsvxx_lacpy(#char "F", n, nrhs, b.*, ldb, x, ldx, MT);
    _gbsvxx_gbtrs(trans, n, kl, ku, nrhs, afb.*, ldafb, ipiv.*, x, ldx, info, MT);
    _gbsvxx_gbrfsx(trans, equed.*, n, kl, ku, nrhs, ab.*, ldab, afb.*, ldafb, ipiv.*, r.*, c.*, b.*, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info, RT, MT);

    if colequ && notran {
        _gbsvxx_lascl2(n, nrhs, c.*, x, ldx, RT, MT);
    } else if rowequ && !notran {
        _gbsvxx_lascl2(n, nrhs, r.*, x, ldx, RT, MT);
    }
}
