// ? LAPACK *gesvxx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgesvxx/cgesvxx/dgesvxx/zgesvxx
// *GESVXX uses the LU factorization to compute the solution to a
// system of linear equations for GE matrices and can return
// normwise and componentwise error bounds.

sgesvxx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    af: *$AF/ArrayType2D,
    $$ldaf: int,
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
    _gesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info, RT=float32, MT=float32);
}

dgesvxx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    af: *$AF/ArrayType2D,
    $$ldaf: int,
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
    _gesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info, RT=float64, MT=float64);
}

cgesvxx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    af: *$AF/ArrayType2D,
    $$ldaf: int,
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
    _gesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float32, MT=C32);
}

zgesvxx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    af: *$AF/ArrayType2D,
    $$ldaf: int,
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
    _gesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gesvxx_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gesvxx_lacpy :: (
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

_gesvxx_laqge :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
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
           A.data_type == MT &&
           R.data_type == RT &&
           C.data_type == RT;
} {
    #if MT == float32 {
        slaqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
    } else if MT == float64 {
        dlaqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
    } else if MT == C32 {
        claqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
    } else {
        zlaqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
    }
}

_gesvxx_geequb :: (
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
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
           A.data_type == MT &&
           R.data_type == RT &&
           C.data_type == RT;
} {
    #if MT == float32 {
        sgeequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    } else if MT == float64 {
        dgeequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    } else if MT == C32 {
        cgeequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    } else {
        zgeequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    }
}

_gesvxx_getrf :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int;
} {
    #if MT == float32 {
        sgetrf(m, n, a, lda, ipiv, info);
    } else if MT == float64 {
        dgetrf(m, n, a, lda, ipiv, info);
    } else if MT == C32 {
        cgetrf(m, n, a, lda, ipiv, info);
    } else {
        zgetrf(m, n, a, lda, ipiv, info);
    }
}

_gesvxx_getrs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    x: *$X/ArrayType2D,
    $$ldx: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sgetrs(trans, n, nrhs, af, ldaf, ipiv, x, ldx, info);
    } else if MT == float64 {
        dgetrs(trans, n, nrhs, af, ldaf, ipiv, x, ldx, info);
    } else if MT == C32 {
        cgetrs(trans, n, nrhs, af, ldaf, ipiv, x, ldx, info);
    } else {
        zgetrs(trans, n, nrhs, af, ldaf, ipiv, x, ldx, info);
    }
}

_gesvxx_lascl2 :: (
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

_gesvxx_gerfsx :: (
    $$trans: u8,
    $$equed: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
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
           A.data_type == MT &&
           AF.data_type == MT &&
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
        sgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info);
    } else if MT == float64 {
        dgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info);
    } else if MT == C32 {
        cgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info);
    } else {
        zgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info);
    }
}

_gesvxx_gerpvgrw :: (
    $$n: int,
    $$ncols: int,
    a: $A/ArrayType2D,
    af: $AF/ArrayType2D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           AF.data_type == MT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    amax, umax : RT = zero_r, zero_r;

    if ncols <= 0 {
        return one_r;
    }

    for j:1..ncols {
        for i:1..n {
            amax = max(amax, abs(fget(a, i, j)));
        }
    }

    for j:1..ncols {
        for i:1..j {
            umax = max(umax, abs(fget(af, i, j)));
        }
    }

    if umax == zero_r {
        return one_r;
    }
    return amax / umax;
}

_gesvxx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    af: *$AF/ArrayType2D,
    $$ldaf: int,
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
           A.data_type == MT &&
           AF.data_type == MT &&
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
        "SGESVXX"
    } else if MT == float64 {
        "DGESVXX"
    } else if MT == C32 {
        "CGESVXX"
    } else {
        "ZGESVXX"
    };

    amax, bignum, colcnd, rcmax, rcmin, rowcnd, smlnum : RT = ---;
    colequ, equil, nofact, notran, rowequ : bool = ---;
    infequ, j : int = ---;

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
        smlnum = _gesvxx_lamch_safe_minimum(RT);
        bignum = one_r / smlnum;
    }

    if !nofact && !equil && !lsame(fact, #char "F") {
        info.* = -1;
    } else if !notran && !lsame(trans, #char "T") && !lsame(trans, #char "C") {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if nrhs < 0 {
        info.* = -4;
    } else if lda < max(1, n) {
        info.* = -6;
    } else if ldaf < max(1, n) {
        info.* = -8;
    } else if lsame(fact, #char "F") && !(rowequ || colequ || lsame(equed.*, #char "N")) {
        info.* = -10;
    } else {
        if rowequ {
            rcmin = bignum;
            rcmax = zero_r;
            for j:1..n {
                rcmin = min(rcmin, fget(r.*, j));
                rcmax = max(rcmax, fget(r.*, j));
            }
            if rcmin <= zero_r {
                info.* = -11;
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
                info.* = -12;
            } else if n > 0 {
                colcnd = max(rcmin, smlnum) / min(rcmax, bignum);
            } else {
                colcnd = one_r;
            }
        }
        if info.* == 0 {
            if ldb < max(1, n) {
                info.* = -14;
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
        _gesvxx_geequb(n, n, a.*, lda, r, c, *rowcnd, *colcnd, *amax, *infequ, RT, MT);
        if infequ == 0 {
            _gesvxx_laqge(n, n, a, lda, r.*, c.*, rowcnd, colcnd, amax, equed, RT, MT);
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
            _gesvxx_lascl2(n, nrhs, r.*, b, ldb, RT, MT);
        }
    } else if colequ {
        _gesvxx_lascl2(n, nrhs, c.*, b, ldb, RT, MT);
    }

    if nofact || equil {
        _gesvxx_lacpy(#char "F", n, n, a.*, lda, af, ldaf, MT);
        _gesvxx_getrf(n, n, af, ldaf, ipiv, info, MT);
        if info.* > 0 {
            rpvgrw.* = _gesvxx_gerpvgrw(n, info.*, a.*, af.*, RT, MT);
            return;
        }
    }

    rpvgrw.* = _gesvxx_gerpvgrw(n, n, a.*, af.*, RT, MT);
    _gesvxx_lacpy(#char "F", n, nrhs, b.*, ldb, x, ldx, MT);
    _gesvxx_getrs(trans, n, nrhs, af.*, ldaf, ipiv.*, x, ldx, info, MT);
    _gesvxx_gerfsx(trans, equed.*, n, nrhs, a.*, lda, af.*, ldaf, ipiv.*, r.*, c.*, b.*, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, aux, info, RT, MT);

    if colequ && notran {
        _gesvxx_lascl2(n, nrhs, c.*, x, ldx, RT, MT);
    } else if rowequ && !notran {
        _gesvxx_lascl2(n, nrhs, r.*, x, ldx, RT, MT);
    }
}
