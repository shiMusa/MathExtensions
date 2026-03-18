// ? LAPACK *gesvx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgesvx/cgesvx/dgesvx/zgesvx
// *GESVX uses the LU factorization to compute the solution to a
// system of linear equations for GE matrices, and provides error
// bounds plus a condition estimate.

sgesvx :: (
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
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, RT=float32, MT=float32);
}

dgesvx :: (
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
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, RT=float64, MT=float64);
}

cgesvx :: (
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
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, RT=float32, MT=C32);
}

zgesvx :: (
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
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gesvx_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gesvx_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_gesvx_lange :: (
    $$norm: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    aux: *$AUX/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           AUX.data_type == RT;
} {
    #if MT == float32 {
        return slange(norm, m, n, a, lda, aux);
    } else if MT == float64 {
        return dlange(norm, m, n, a, lda, aux);
    } else if MT == C32 {
        return clange(norm, m, n, a, lda, aux);
    } else {
        return zlange(norm, m, n, a, lda, aux);
    }
}

_gesvx_lacpy :: (
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

_gesvx_geequ :: (
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
        sgeequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    } else if MT == float64 {
        dgeequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    } else if MT == C32 {
        cgeequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    } else {
        zgeequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    }
}

_gesvx_laqge :: (
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

_gesvx_getrf :: (
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

_gesvx_getrs :: (
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

_gesvx_gecon :: (
    $$norm: u8,
    $$n: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    anorm: $RT,
    rcond: *RT,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AF.data_type == MT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        sgecon(norm, n, af, ldaf, anorm, rcond, work, aux, info);
    } else if MT == float64 {
        dgecon(norm, n, af, ldaf, anorm, rcond, work, aux, info);
    } else if MT == C32 {
        cgecon(norm, n, af, ldaf, anorm, rcond, work, aux, info);
    } else {
        zgecon(norm, n, af, ldaf, anorm, rcond, work, aux, info);
    }
}

_gesvx_gerfs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
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
           A.data_type == MT &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT &&
           X.data_type == MT &&
           FERR.data_type == RT &&
           BERR.data_type == RT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        sgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    } else if MT == float64 {
        dgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    } else if MT == C32 {
        cgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    } else {
        zgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    }
}

_gesvx_rpvgrw_partial :: (
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

_gesvx :: (
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
           A.data_type == MT &&
           AF.data_type == MT &&
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
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == float32 {
        "SGESVX"
    } else if MT == float64 {
        "DGESVX"
    } else if MT == C32 {
        "CGESVX"
    } else {
        "ZGESVX"
    };

    amax, anorm, bignum, colcnd, rcmax, rcmin, rowcnd, smlnum : RT = ---;
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
        smlnum = _gesvx_lamch_safe_minimum(RT);
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
        _gesvx_geequ(n, n, a.*, lda, r, c, *rowcnd, *colcnd, *amax, *infequ, RT, MT);
        if infequ == 0 {
            _gesvx_laqge(n, n, a, lda, r.*, c.*, rowcnd, colcnd, amax, equed, RT, MT);
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
        _gesvx_lacpy(#char "F", n, n, a.*, lda, af, ldaf, MT);
        _gesvx_getrf(n, n, af, ldaf, ipiv, info, MT);

        if info.* > 0 {
            rcond.* = zero_r;
            #if MT == float32 || MT == float64 {
                fset(work, 1, _gesvx_rpvgrw_partial(n, info.*, a.*, af.*, RT, MT));
            } else {
                fset(aux, 1, _gesvx_rpvgrw_partial(n, info.*, a.*, af.*, RT, MT));
            }
            return;
        }
    }

    if notran {
        norm = #char "1";
    } else {
        norm = #char "I";
    }

    #if MT == float32 || MT == float64 {
        anorm = _gesvx_lange(norm, n, n, a.*, lda, work, RT, MT);
    } else {
        anorm = _gesvx_lange(norm, n, n, a.*, lda, aux, RT, MT);
    }

    _gesvx_gecon(norm, n, af.*, ldaf, anorm, rcond, work, aux, info, RT, MT);
    _gesvx_lacpy(#char "F", n, nrhs, b.*, ldb, x, ldx, MT);
    _gesvx_getrs(trans, n, nrhs, af.*, ldaf, ipiv.*, x, ldx, info, MT);
    _gesvx_gerfs(trans, n, nrhs, a.*, lda, af.*, ldaf, ipiv.*, b.*, ldb, x, ldx, ferr, berr, work, aux, info, RT, MT);

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

    if rcond.* < _gesvx_lamch_epsilon(RT) {
        info.* = n + 1;
    }

    #if MT == float32 || MT == float64 {
        fset(work, 1, _gesvx_rpvgrw_partial(n, n, a.*, af.*, RT, MT));
    } else {
        fset(aux, 1, _gesvx_rpvgrw_partial(n, n, a.*, af.*, RT, MT));
    }
}
