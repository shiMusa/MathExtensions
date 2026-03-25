// ? LAPACK *lanhf()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clanhf/zlanhf
// *LANHF returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// Hermitian matrix in RFP format.

clanhf :: #bake_arguments _lanhf(RT=float32, MT=C32);
zlanhf :: #bake_arguments _lanhf(RT=float64, MT=C64);

#scope_file
_lanhf_tfttr :: (
    $$transr: u8,
    $$uplo: u8,
    $$n: int,
    arf: $ARF/ArrayType1D,
    a: *$A/ArrayType2D,
    $$lda: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           ARF.data_type == MT &&
           A.data_type == MT;
} {
    #if MT == C32 {
        ctfttr(transr, uplo, n, arf, a, lda, info);
    } else {
        ztfttr(transr, uplo, n, arf, a, lda, info);
    }
}

_lanhf_lanhe :: (
    $$norm: u8,
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           WORK.data_type == RT;
} {
    #if MT == C32 {
        return clanhe(norm, uplo, n, a, lda, work);
    } else {
        return zlanhe(norm, uplo, n, a, lda, work);
    }
}

_lanhf :: (
    $$norm: u8,
    $$transr: u8,
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           WORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == C32 { "CLANHF" } else { "ZLANHF" };
    info: int = 0;

    if n <= 0 {
        return zero_r;
    }

    tmp := dhmat(MT, n, n, false);
    defer free(tmp);

    _lanhf_tfttr(transr, uplo, n, a, *tmp, n, *info, MT);
    if info != 0 {
        xerbla(routine_name, -info);
        return zero_r;
    }

    return _lanhf_lanhe(norm, uplo, n, tmp, n, work, RT, MT);
}
