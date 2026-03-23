// ? LAPACK *hpsv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chpsv/zhpsv
// *HPSV computes the solution to a system of linear equations
// A * X = B for Hermitian packed matrices.

chpsv :: #bake_arguments _hpsv(MT=C32);
zhpsv :: #bake_arguments _hpsv(MT=C64);

#scope_file
_hpsv_hptrf :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, ipiv: *$IPIV/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && IPIV.data_type == int;
} {
    #if MT == C32 { chptrf(uplo, n, ap, ipiv, info); } else { zhptrf(uplo, n, ap, ipiv, info); }
}

_hpsv_hptrs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, ap: *$AP/ArrayType1D, ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    #if MT == C32 { chptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info); } else { zhptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info); }
}

_hpsv :: (
    $$uplo: u8, $$n: int, $$nrhs: int, ap: *$AP/ArrayType1D, ipiv: *$IPIV/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    routine_name :: #run if MT == C32 { "CHPSV " } else { "ZHPSV " };

    info.* = 0;
    if !lsame(uplo, #char "U") && !lsame(uplo, #char "L") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if ldb < max(1, n) {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    _hpsv_hptrf(uplo, n, ap, ipiv, info, MT);
    if info.* == 0 {
        _hpsv_hptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info, MT);
    }
}
