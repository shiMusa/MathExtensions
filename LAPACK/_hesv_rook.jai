// ? LAPACK *hesv_rook()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chesv_rook/zhesv_rook
// *HESV_ROOK computes the solution to a system of linear equations
// A * X = B for Hermitian matrices using rook pivoting.

chesv_rook :: #bake_arguments _hesv_rook(RT=float32, MT=C32);
zhesv_rook :: #bake_arguments _hesv_rook(RT=float64, MT=C64);

#scope_file
_hesv_rook_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hesv_rook_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hesv_rook_hetrf :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    #if MT == C32 {
        chetrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
    } else {
        zhetrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
    }
}

_hesv_rook_hetrs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int,
    info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    #if MT == C32 {
        chetrs_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, info);
    } else {
        zhetrs_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, info);
    }
}

_hesv_rook :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && WORK.data_type == MT;
} {
    routine_name :: #run if MT == C32 { "CHESV_ROOK " } else { "ZHESV_ROOK " };

    lquery, upper : bool = ---;
    lwkopt, nb : int = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ldb < max(1, n) {
        info.* = -8;
    } else if lwork < 1 && !lquery {
        info.* = -10;
    }

    if info.* == 0 {
        if n == 0 {
            lwkopt = 1;
        } else {
            nb = ilaenv(1, ifx MT == C32 then "CHETRF_ROOK" else "ZHETRF_ROOK", cast(*u8) *uplo, n, -1, -1, -1);
            lwkopt = n*nb;
        }
        fset(work, 1, _hesv_rook_complex_from_real(_hesv_rook_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    _hesv_rook_hetrf(uplo, n, a, lda, ipiv, work, lwork, info, MT);
    if info.* == 0 {
        _hesv_rook_hetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info, MT);
    }

    fset(work, 1, _hesv_rook_complex_from_real(_hesv_rook_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
}
