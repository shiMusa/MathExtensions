// ? LAPACK *hetri_3x()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetri_3x/zhetri_3x
// *HETRI_3X computes the inverse of a complex Hermitian indefinite matrix
// using the blocked workspace kernel behind *HETRI_3.

chetri_3x :: #bake_arguments _hetri_3x(RT=float32, MT=C32);
zhetri_3x :: #bake_arguments _hetri_3x(RT=float64, MT=C64);

#scope_file
_hetri_3x :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    e: $E/ArrayType1D, ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D,
    $$nb: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && E.data_type == MT &&
           IPIV.data_type == int && WORK.data_type == MT;
} {
    routine_name :: #run if MT == C32 { "CHETRI_3X" } else { "ZHETRI_3X" };
    upper := lsame(uplo, #char "U");

    info.* = 0;
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if nb < 1 {
        info.* = -8;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }

    // Staged port note:
    // The blocked *HETRI_3X kernel still needs its own literal follow-up pass.
    // This file establishes the family interface and validation shape.
}
