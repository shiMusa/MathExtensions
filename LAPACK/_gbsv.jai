// ? LAPACK *gbsv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbsv/cgbsv/dgbsv/zgbsv
// *GBSV computes the solution to a system of linear equations
// A * X = B for GB matrices.

sgbsv :: (
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info, T=float32);
}

dgbsv :: (
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info, T=float64);
}

cgbsv :: (
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info, T=C32);
}

zgbsv :: (
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info, T=C64);
}

#scope_file
_gbsv_gbtrf :: (
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $T: Type
) #modify {
    return (T == float32 || T == float64 || T == C32 || T == C64) &&
           AB.data_type == T &&
           IPIV.data_type == int;
} {
    #if T == float32 {
        sgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    } else if T == float64 {
        dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    } else if T == C32 {
        cgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    } else {
        zgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    }
}

_gbsv_gbtrs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $T: Type
) #modify {
    return (T == float32 || T == float64 || T == C32 || T == C64) &&
           AB.data_type == T &&
           IPIV.data_type == int &&
           B.data_type == T;
} {
    #if T == float32 {
        sgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    } else if T == float64 {
        dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    } else if T == C32 {
        cgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    } else {
        zgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    }
}

_gbsv :: (
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: *$AB/ArrayType2D,
    $$ldab: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $T: Type
) #modify {
    return (T == float32 || T == float64 || T == C32 || T == C64) &&
           AB.data_type == T &&
           IPIV.data_type == int &&
           B.data_type == T;
} {
    routine_name :: #run if T == float32 {
        "SGBSV "
    } else if T == float64 {
        "DGBSV "
    } else if T == C32 {
        "CGBSV "
    } else {
        "ZGBSV "
    };

    info.* = 0;
    if n < 0 {
        info.* = -1;
    } else if kl < 0 {
        info.* = -2;
    } else if ku < 0 {
        info.* = -3;
    } else if nrhs < 0 {
        info.* = -4;
    } else if ldab < 2*kl + ku + 1 {
        info.* = -6;
    } else if ldb < max(n, 1) {
        info.* = -9;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    _gbsv_gbtrf(n, n, kl, ku, ab, ldab, ipiv, info, T);
    if info.* == 0 {
        _gbsv_gbtrs(#char "N", n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info, T);
    }
}
