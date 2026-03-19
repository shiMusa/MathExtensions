// ? LAPACK *gtsvx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgtsvx/cgtsvx/dgtsvx/zgtsvx
// *GTSVX uses the LU factorization to solve a tridiagonal system and
// provides error bounds plus a condition estimate.

sgtsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: *$DLF/ArrayType1D,
    df: *$DF/ArrayType1D,
    duf: *$DUF/ArrayType1D,
    du2: *$DU2/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    b: $B/ArrayType2D,
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
    _gtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, RT=float32, MT=float32);
}

dgtsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: *$DLF/ArrayType1D,
    df: *$DF/ArrayType1D,
    duf: *$DUF/ArrayType1D,
    du2: *$DU2/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    b: $B/ArrayType2D,
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
    _gtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, RT=float64, MT=float64);
}

cgtsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: *$DLF/ArrayType1D,
    df: *$DF/ArrayType1D,
    duf: *$DUF/ArrayType1D,
    du2: *$DU2/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    b: $B/ArrayType2D,
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
    _gtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, RT=float32, MT=C32);
}

zgtsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: *$DLF/ArrayType1D,
    df: *$DF/ArrayType1D,
    duf: *$DUF/ArrayType1D,
    du2: *$DU2/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    b: $B/ArrayType2D,
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
    _gtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gtsvx_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_gtsvx_langt :: (
    $$norm: u8,
    $$n: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT;
} {
    #if MT == float32 {
        return slangt(norm, n, dl, d, du);
    } else if MT == float64 {
        return dlangt(norm, n, dl, d, du);
    } else if MT == C32 {
        return clangt(norm, n, dl, d, du);
    } else {
        return zlangt(norm, n, dl, d, du);
    }
}

_gtsvx_copy :: (
    $$n: int,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        scopy(n, x, 1, y, 1);
    } else if MT == float64 {
        dcopy(n, x, 1, y, 1);
    } else if MT == C32 {
        ccopy(n, x, 1, y, 1);
    } else {
        zcopy(n, x, 1, y, 1);
    }
}

_gtsvx_lacpy :: (
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

_gtsvx_gttrf :: (
    $$n: int,
    dl: *$DL/ArrayType1D,
    d: *$D/ArrayType1D,
    du: *$DU/ArrayType1D,
    du2: *$DU2/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int;
} {
    #if MT == float32 {
        sgttrf(n, dl, d, du, du2, ipiv, info);
    } else if MT == float64 {
        dgttrf(n, dl, d, du, du2, ipiv, info);
    } else if MT == C32 {
        cgttrf(n, dl, d, du, du2, ipiv, info);
    } else {
        zgttrf(n, dl, d, du, du2, ipiv, info);
    }
}

_gtsvx_gttrs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    #if MT == float32 {
        sgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info);
    } else if MT == float64 {
        dgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info);
    } else if MT == C32 {
        cgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info);
    } else {
        zgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info);
    }
}

_gtsvx_gtcon :: (
    $$norm: u8,
    $$n: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    du2: $DU2/ArrayType1D,
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
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        sgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, aux, info);
    } else if MT == float64 {
        dgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, aux, info);
    } else if MT == C32 {
        cgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info);
    } else {
        zgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info);
    }
}

_gtsvx_gtrfs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: $DLF/ArrayType1D,
    df: $DF/ArrayType1D,
    duf: $DUF/ArrayType1D,
    du2: $DU2/ArrayType1D,
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
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DLF.data_type == MT &&
           DF.data_type == MT &&
           DUF.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT &&
           X.data_type == MT &&
           FERR.data_type == RT &&
           BERR.data_type == RT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        sgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    } else if MT == float64 {
        dgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    } else if MT == C32 {
        cgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    } else {
        zgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, aux, info);
    }
}

_gtsvx :: (
    $$fact: u8,
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: *$DLF/ArrayType1D,
    df: *$DF/ArrayType1D,
    duf: *$DUF/ArrayType1D,
    du2: *$DU2/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    b: $B/ArrayType2D,
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
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DLF.data_type == MT &&
           DF.data_type == MT &&
           DUF.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT &&
           X.data_type == MT &&
           FERR.data_type == RT &&
           BERR.data_type == RT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == float32 {
        "SGTSVX"
    } else if MT == float64 {
        "DGTSVX"
    } else if MT == C32 {
        "CGTSVX"
    } else {
        "ZGTSVX"
    };

    nofact, notran : bool = ---;
    norm : u8 = ---;
    anorm : RT = ---;

    info.* = 0;
    nofact = lsame(fact, #char "N");
    notran = lsame(trans, #char "N");
    if !nofact && !lsame(fact, #char "F") {
        info.* = -1;
    } else if !notran && !lsame(trans, #char "T") && !lsame(trans, #char "C") {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if nrhs < 0 {
        info.* = -4;
    } else if ldb < max(1, n) {
        info.* = -14;
    } else if ldx < max(1, n) {
        info.* = -16;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if nofact {
        _gtsvx_copy(n, d, df, MT);
        if n > 1 {
            _gtsvx_copy(n - 1, dl, dlf, MT);
            _gtsvx_copy(n - 1, du, duf, MT);
        }
        _gtsvx_gttrf(n, dlf, df, duf, du2, ipiv, info, RT, MT);
        if info.* > 0 {
            rcond.* = zero_r;
            return;
        }
    }

    if notran {
        norm = #char "1";
    } else {
        norm = #char "I";
    }
    anorm = _gtsvx_langt(norm, n, dl, d, du, RT, MT);

    _gtsvx_gtcon(norm, n, dlf.*, df.*, duf.*, du2.*, ipiv.*, anorm, rcond, work, aux, info, RT, MT);

    _gtsvx_lacpy(#char "F", n, nrhs, b, ldb, x, ldx, MT);
    _gtsvx_gttrs(trans, n, nrhs, dlf.*, df.*, duf.*, du2.*, ipiv.*, x, ldx, info, RT, MT);

    _gtsvx_gtrfs(trans, n, nrhs, dl, d, du, dlf.*, df.*, duf.*, du2.*, ipiv.*, b, ldb, x, ldx, ferr, berr, work, aux, info, RT, MT);

    if rcond.* < _gtsvx_lamch_epsilon(RT) {
        info.* = n + 1;
    }
}
