// ? LAPACK *hpsvx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chpsvx/zhpsvx
// *HPSVX uses the diagonal pivoting factorization to compute the
// solution to a system of linear equations for Hermitian packed matrices.

chpsvx :: (
    $$fact: u8, $$uplo: u8, $$n: int, $$nrhs: int,
    ap: $AP/ArrayType1D, afp: *$AFP/ArrayType1D, ipiv: *$IPIV/ArrayType1D,
    b: $B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    rcond: *float32, ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D, info: *int
) {
    _hpsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, RT=float32, MT=C32);
}

zhpsvx :: (
    $$fact: u8, $$uplo: u8, $$n: int, $$nrhs: int,
    ap: $AP/ArrayType1D, afp: *$AFP/ArrayType1D, ipiv: *$IPIV/ArrayType1D,
    b: $B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    rcond: *float64, ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D, info: *int
) {
    _hpsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_hpsvx_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(#char "e"); } else { return dlamch(#char "e"); }
}

_hpsvx_lanhp :: (
    $$norm: u8, $$uplo: u8, $$n: int, ap: $AP/ArrayType1D, rwork: *$RWORK/ArrayType1D, $RT: Type, $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { return clanhp(norm, uplo, n, ap, rwork); } else { return zlanhp(norm, uplo, n, ap, rwork); }
}

_hpsvx_copy_packed :: (
    $$n: int, x: $X/ArrayType1D, y: *$Y/ArrayType1D, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    count := (n*(n + 1)) / 2;
    #if MT == C32 { ccopy(count, x, 1, y, 1); } else { zcopy(count, x, 1, y, 1); }
}

_hpsvx_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 { clacpy(uplo, m, n, a, lda, b, ldb); } else { zlacpy(uplo, m, n, a, lda, b, ldb); }
}

_hpsvx_hptrf :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, ipiv: *$IPIV/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && IPIV.data_type == int;
} {
    #if MT == C32 { chptrf(uplo, n, ap, ipiv, info); } else { zhptrf(uplo, n, ap, ipiv, info); }
}

_hpsvx_hptrs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, ap: $AP/ArrayType1D, ipiv: $IPIV/ArrayType1D,
    x: *$X/ArrayType2D, $$ldx: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && IPIV.data_type == int && X.data_type == MT;
} {
    #if MT == C32 { chptrs(uplo, n, nrhs, ap, ipiv, x, ldx, info); } else { zhptrs(uplo, n, nrhs, ap, ipiv, x, ldx, info); }
}

_hpsvx_hpcon :: (
    $$uplo: u8, $$n: int, ap: $AP/ArrayType1D, ipiv: $IPIV/ArrayType1D, anorm: $RT, rcond: *RT,
    work: *$WORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    #if MT == C32 { chpcon(uplo, n, ap, ipiv, anorm, rcond, work, info); } else { zhpcon(uplo, n, ap, ipiv, anorm, rcond, work, info); }
}

_hpsvx_hprfs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, ap: $AP/ArrayType1D, afp: $AFP/ArrayType1D, ipiv: $IPIV/ArrayType1D,
    b: $B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int, ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D, work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && AFP.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && X.data_type == MT && FERR.data_type == RT && BERR.data_type == RT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { chprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info); }
    else { zhprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info); }
}

_hpsvx :: (
    $$fact: u8, $$uplo: u8, $$n: int, $$nrhs: int,
    ap: $AP/ArrayType1D, afp: *$AFP/ArrayType1D, ipiv: *$IPIV/ArrayType1D,
    b: $B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    rcond: *$RT, ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && AFP.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && X.data_type == MT &&
           FERR.data_type == RT && BERR.data_type == RT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    anorm : RT = ---;
    nofact : bool = ---;
    routine_name :: #run if MT == C32 { "CHPSVX" } else { "ZHPSVX" };

    info.* = 0;
    nofact = lsame(fact, #char "N");
    if !nofact && !lsame(fact, #char "F") {
        info.* = -1;
    } else if !lsame(uplo, #char "U") && !lsame(uplo, #char "L") {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if nrhs < 0 {
        info.* = -4;
    } else if ldb < max(1, n) {
        info.* = -9;
    } else if ldx < max(1, n) {
        info.* = -11;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if nofact {
        _hpsvx_copy_packed(n, ap, afp, MT);
        _hpsvx_hptrf(uplo, n, afp, ipiv, info, MT);
        if info.* > 0 {
            rcond.* = zero_r;
            return;
        }
    }

    anorm = _hpsvx_lanhp(#char "I", uplo, n, ap, rwork, RT, MT);
    _hpsvx_hpcon(uplo, n, afp.*, ipiv.*, anorm, rcond, work, info, RT, MT);
    _hpsvx_lacpy(#char "F", n, nrhs, b, ldb, x, ldx, MT);
    _hpsvx_hptrs(uplo, n, nrhs, afp.*, ipiv.*, x, ldx, info, MT);
    _hpsvx_hprfs(uplo, n, nrhs, ap, afp.*, ipiv.*, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT, MT);

    if rcond.* < _hpsvx_lamch_epsilon(RT) {
        info.* = n + 1;
    }
}
