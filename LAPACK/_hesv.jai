// ? LAPACK *hesv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chesv/zhesv
// *HESV computes the solution to a system of linear equations
// A * X = B for Hermitian matrices.

chesv :: #bake_arguments _hesv(RT=float32, MT=C32);
zhesv :: #bake_arguments _hesv(RT=float64, MT=C64);

#scope_file
_hesv_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hesv_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hesv_hetrf :: (
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

_hesv_hetrs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int,
    info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    #if MT == C32 {
        chetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info);
    } else {
        zhetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info);
    }
}

_hesv_hetrs2 :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int,
    work: *$WORK/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        chetrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info);
    } else {
        zhetrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info);
    }
}

_hesv :: (
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
    routine_name :: #run if MT == C32 { "CHESV " } else { "ZHESV " };

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
            nb = ilaenv(1, ifx MT == C32 then "CHETRF" else "ZHETRF", cast(*u8) *uplo, n, -1, -1, -1);
            lwkopt = n*nb;
        }
        fset(work, 1, _hesv_complex_from_real(_hesv_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    _hesv_hetrf(uplo, n, a, lda, ipiv, work, lwork, info, MT);
    if info.* == 0 {
        if lwork < n {
            _hesv_hetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info, MT);
        } else {
            _hesv_hetrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info, MT);
        }
    }

    fset(work, 1, _hesv_complex_from_real(_hesv_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
}
