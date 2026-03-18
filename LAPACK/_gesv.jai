// ? LAPACK *gesv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgesv/cgesv/dgesv/zgesv
// *GESV computes the solution to a system of linear equations
// A * X = B for GE matrices.

sgesv :: (
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gesv(n, nrhs, a, lda, ipiv, b, ldb, info, T=float32);
}

dgesv :: (
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gesv(n, nrhs, a, lda, ipiv, b, ldb, info, T=float64);
}

cgesv :: (
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gesv(n, nrhs, a, lda, ipiv, b, ldb, info, T=C32);
}

zgesv :: (
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int
) {
    _gesv(n, nrhs, a, lda, ipiv, b, ldb, info, T=C64);
}

#scope_file
_gesv_getrf :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $T: Type
) #modify {
    return (T == float32 || T == float64 || T == C32 || T == C64) &&
           A.data_type == T &&
           IPIV.data_type == int;
} {
    #if T == float32 {
        sgetrf(m, n, a, lda, ipiv, info);
    } else if T == float64 {
        dgetrf(m, n, a, lda, ipiv, info);
    } else if T == C32 {
        cgetrf(m, n, a, lda, ipiv, info);
    } else {
        zgetrf(m, n, a, lda, ipiv, info);
    }
}

_gesv_getrs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $T: Type
) #modify {
    return (T == float32 || T == float64 || T == C32 || T == C64) &&
           A.data_type == T &&
           IPIV.data_type == int &&
           B.data_type == T;
} {
    #if T == float32 {
        sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
    } else if T == float64 {
        dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
    } else if T == C32 {
        cgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
    } else {
        zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
    }
}

_gesv :: (
    $$n: int,
    $$nrhs: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $T: Type
) #modify {
    return (T == float32 || T == float64 || T == C32 || T == C64) &&
           A.data_type == T &&
           IPIV.data_type == int &&
           B.data_type == T;
} {
    routine_name :: #run if T == float32 {
        "SGESV "
    } else if T == float64 {
        "DGESV "
    } else if T == C32 {
        "CGESV "
    } else {
        "ZGESV "
    };

    info.* = 0;
    if n < 0 {
        info.* = -1;
    } else if nrhs < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if ldb < max(1, n) {
        info.* = -7;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    _gesv_getrf(n, n, a, lda, ipiv, info, T);
    if info.* == 0 {
        _gesv_getrs(#char "N", n, nrhs, a, lda, ipiv, b, ldb, info, T);
    }
}
