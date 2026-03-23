// ? LAPACK *hesvx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chesvx/zhesvx
// *HESVX uses the diagonal pivoting factorization to compute the
// solution to a system of linear equations for Hermitian matrices.

chesvx :: (
    $$fact: u8, $$uplo: u8, $$n: int, $$nrhs: int,
    a: $A/ArrayType2D, $$lda: int, af: *$AF/ArrayType2D, $$ldaf: int,
    ipiv: *$IPIV/ArrayType1D, b: $B/ArrayType2D, $$ldb: int,
    x: *$X/ArrayType2D, $$ldx: int, rcond: *float32,
    ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _hesvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info, RT=float32, MT=C32);
}

zhesvx :: (
    $$fact: u8, $$uplo: u8, $$n: int, $$nrhs: int,
    a: $A/ArrayType2D, $$lda: int, af: *$AF/ArrayType2D, $$ldaf: int,
    ipiv: *$IPIV/ArrayType1D, b: $B/ArrayType2D, $$ldb: int,
    x: *$X/ArrayType2D, $$ldx: int, rcond: *float64,
    ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _hesvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info, RT=float64, MT=C64);
}

#scope_file
_hesvx_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hesvx_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hesvx_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "e"); } else { return dlamch(#char "e"); }
}

_hesvx_lanhe :: (
    $$norm: u8, $$uplo: u8, $$n: int, a: $A/ArrayType2D, $$lda: int,
    rwork: *$RWORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 {
        return clanhe(norm, uplo, n, a, lda, rwork);
    } else {
        return zlanhe(norm, uplo, n, a, lda, rwork);
    }
}

_hesvx_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_hesvx_hetrf :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    #if MT == C32 {
        chetrf(uplo, n, a, lda, ipiv, work, lwork, info);
    } else {
        zhetrf(uplo, n, a, lda, ipiv, work, lwork, info);
    }
}

_hesvx_hetrs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, af: $AF/ArrayType2D, $$ldaf: int,
    ipiv: $IPIV/ArrayType1D, x: *$X/ArrayType2D, $$ldx: int, info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           AF.data_type == MT && IPIV.data_type == int && X.data_type == MT;
} {
    #if MT == C32 {
        chetrs(uplo, n, nrhs, af, ldaf, ipiv, x, ldx, info);
    } else {
        zhetrs(uplo, n, nrhs, af, ldaf, ipiv, x, ldx, info);
    }
}

_hesvx_hecon :: (
    $$uplo: u8, $$n: int, af: $AF/ArrayType2D, $$ldaf: int, ipiv: $IPIV/ArrayType1D,
    anorm: $RT, rcond: *RT, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AF.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    #if MT == C32 {
        checon(uplo, n, af, ldaf, ipiv, anorm, rcond, work, info);
    } else {
        zhecon(uplo, n, af, ldaf, ipiv, anorm, rcond, work, info);
    }
}

_hesvx_herfs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: $A/ArrayType2D, $$lda: int,
    af: $AF/ArrayType2D, $$ldaf: int, ipiv: $IPIV/ArrayType1D,
    b: $B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && AF.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && X.data_type == MT &&
           FERR.data_type == RT && BERR.data_type == RT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 {
        cherfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    } else {
        zherfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    }
}

_hesvx :: (
    $$fact: u8, $$uplo: u8, $$n: int, $$nrhs: int,
    a: $A/ArrayType2D, $$lda: int, af: *$AF/ArrayType2D, $$ldaf: int,
    ipiv: *$IPIV/ArrayType1D, b: $B/ArrayType2D, $$ldb: int,
    x: *$X/ArrayType2D, $$ldx: int, rcond: *$RT,
    ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && AF.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && X.data_type == MT &&
           FERR.data_type == RT && BERR.data_type == RT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == C32 { "CHESVX" } else { "ZHESVX" };

    anorm : RT = ---;
    lquery, nofact, upper : bool = ---;
    lwkmin, lwkopt, nb : int = ---;

    info.* = 0;
    nofact = lsame(fact, #char "N");
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;
    lwkmin = max(1, 2*n);
    if !nofact && !lsame(fact, #char "F") {
        info.* = -1;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if nrhs < 0 {
        info.* = -4;
    } else if lda < max(1, n) {
        info.* = -6;
    } else if ldaf < max(1, n) {
        info.* = -8;
    } else if ldb < max(1, n) {
        info.* = -11;
    } else if ldx < max(1, n) {
        info.* = -13;
    } else if lwork < lwkmin && !lquery {
        info.* = -18;
    }

    if info.* == 0 {
        lwkopt = lwkmin;
        if nofact {
            nb = ilaenv(1, ifx MT == C32 then "CHETRF" else "ZHETRF", cast(*u8) *uplo, n, -1, -1, -1);
            lwkopt = max(lwkopt, n*nb);
        }
        fset(work, 1, _hesvx_complex_from_real(_hesvx_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if nofact {
        _hesvx_lacpy(uplo, n, n, a, lda, af, ldaf, MT);
        _hesvx_hetrf(uplo, n, af, ldaf, ipiv, work, lwork, info, MT);
        if info.* > 0 {
            rcond.* = zero_r;
            return;
        }
    }

    anorm = _hesvx_lanhe(#char "I", uplo, n, a, lda, rwork, RT, MT);
    _hesvx_hecon(uplo, n, af.*, ldaf, ipiv.*, anorm, rcond, work, info, RT, MT);
    _hesvx_lacpy(#char "F", n, nrhs, b, ldb, x, ldx, MT);
    _hesvx_hetrs(uplo, n, nrhs, af.*, ldaf, ipiv.*, x, ldx, info, MT);
    _hesvx_herfs(uplo, n, nrhs, a, lda, af.*, ldaf, ipiv.*, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT, MT);

    if rcond.* < _hesvx_lamch_epsilon(RT) {
        info.* = n + 1;
    }

    fset(work, 1, _hesvx_complex_from_real(_hesvx_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
}
