// ? LAPACK *getrs()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgetrs/cgetrs/dgetrs/zgetrs
// *GETRS solves a system of linear equations using the LU
// factorization computed by *GETRF.

sgetrs :: #bake_arguments _getrs(MT=float32);
dgetrs :: #bake_arguments _getrs(MT=float64);
cgetrs :: #bake_arguments _getrs(MT=C32);
zgetrs :: #bake_arguments _getrs(MT=C64);

#scope_file
_getrs_laswp :: (
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $$k1: int,
    $$k2: int,
    ipiv: $IPIV/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int;
} {
    #if MT == float32 {
        slaswp(n, a, lda, k1, k2, ipiv, incx);
    } else if MT == float64 {
        dlaswp(n, a, lda, k1, k2, ipiv, incx);
    } else if MT == C32 {
        claswp(n, a, lda, k1, k2, ipiv, incx);
    } else {
        zlaswp(n, a, lda, k1, k2, ipiv, incx);
    }
}

_getrs_trsm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
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
        strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == float64 {
        dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_getrs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    one_t :: #run one(MT);
    routine_name :: #run if MT == float32 {
        "SGETRS"
    } else if MT == float64 {
        "DGETRS"
    } else if MT == C32 {
        "CGETRS"
    } else {
        "ZGETRS"
    };

    notran := lsame(trans, #char "N");

    info.* = 0;
    if !notran && !lsame(trans, #char "T") && !lsame(trans, #char "C") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ldb < max(1, n) {
        info.* = -8;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 || nrhs == 0 then return;

    if notran {
        _getrs_laswp(nrhs, b, ldb, 1, n, ipiv, 1, MT);
        _getrs_trsm(#char "L", #char "L", #char "N", #char "U", n, nrhs, one_t, a, lda, b, ldb, MT);
        _getrs_trsm(#char "L", #char "U", #char "N", #char "N", n, nrhs, one_t, a, lda, b, ldb, MT);
    } else {
        _getrs_trsm(#char "L", #char "U", trans, #char "N", n, nrhs, one_t, a, lda, b, ldb, MT);
        _getrs_trsm(#char "L", #char "L", trans, #char "U", n, nrhs, one_t, a, lda, b, ldb, MT);
        _getrs_laswp(nrhs, b, ldb, 1, n, ipiv, -1, MT);
    }
}
