// ? LAPACK *hetri_3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetri_3/zhetri_3
// *HETRI_3 computes the inverse of a complex Hermitian indefinite matrix
// using the factorization from *HETRF_RK or *HETRF_BK and dispatches to *HETRI_3X.

chetri_3 :: #bake_arguments _hetri_3(RT=float32, MT=C32);
zhetri_3 :: #bake_arguments _hetri_3(RT=float64, MT=C64);

#scope_file
_hetri_3_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetri_3_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetri_3 :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    e: $E/ArrayType1D, ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D,
    $$lwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && E.data_type == MT &&
           IPIV.data_type == int && WORK.data_type == MT;
} {
    routine_name :: #run if MT == C32 { "CHETRI_3" } else { "ZHETRI_3" };

    lquery, upper : bool = ---;
    lwkopt, nb : int = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;

    if n == 0 {
        lwkopt = 1;
    } else {
        nb = max(1, ilaenv(1, ifx MT == C32 then "CHETRI_3" else "ZHETRI_3", cast(*u8) *uplo, n, -1, -1, -1));
        lwkopt = (n + nb + 1)*(nb + 3);
    }
    fset(work, 1, _hetri_3_complex_from_real(_hetri_3_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));

    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if lwork < lwkopt && !lquery {
        info.* = -8;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        return;
    }

    #if MT == C32 {
        chetri_3x(uplo, n, a, lda, e, ipiv, work, nb, info);
    } else {
        zhetri_3x(uplo, n, a, lda, e, ipiv, work, nb, info);
    }

    fset(work, 1, _hetri_3_complex_from_real(_hetri_3_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
}
