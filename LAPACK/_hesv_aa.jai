// ? LAPACK *hesv_aa()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chesv_aa/zhesv_aa
// *HESV_AA computes the solution to a system of linear equations
// A * X = B for Hermitian matrices using Aasen's algorithm.

chesv_aa :: #bake_arguments _hesv_aa(RT=float32, MT=C32);
zhesv_aa :: #bake_arguments _hesv_aa(RT=float64, MT=C64);

#scope_file
_hesv_aa_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hesv_aa_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hesv_aa_hetrf :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    #if MT == C32 {
        chetrf_aa(uplo, n, a, lda, ipiv, work, lwork, info);
    } else {
        zhetrf_aa(uplo, n, a, lda, ipiv, work, lwork, info);
    }
}

_hesv_aa_hetrs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        chetrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
    } else {
        zhetrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
    }
}

_hesv_aa :: (
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
    routine_name :: #run if MT == C32 { "CHESV_AA " } else { "ZHESV_AA " };

    lquery, upper : bool = ---;
    lwkmin, lwkopt, lwkopt_hetrf, lwkopt_hetrs : int = ---;
    hetrf_info, hetrs_info : int = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;
    lwkmin = max(1, 2*n, 3*n - 2);
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
    } else if lwork < lwkmin && !lquery {
        info.* = -10;
    }

    if info.* == 0 {
        _hesv_aa_hetrf(uplo, n, a, lda, ipiv, work, -1, *hetrf_info, MT);
        lwkopt_hetrf = cast(int) real(fget(work.*, 1));
        _hesv_aa_hetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, -1, *hetrs_info, MT);
        lwkopt_hetrs = cast(int) real(fget(work.*, 1));
        lwkopt = max(lwkmin, lwkopt_hetrf, lwkopt_hetrs);
        fset(work, 1, _hesv_aa_complex_from_real(_hesv_aa_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    _hesv_aa_hetrf(uplo, n, a, lda, ipiv, work, lwork, info, MT);
    if info.* == 0 {
        _hesv_aa_hetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, MT);
    }

    fset(work, 1, _hesv_aa_complex_from_real(_hesv_aa_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
}
