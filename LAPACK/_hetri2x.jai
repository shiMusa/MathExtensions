// ? LAPACK *hetri2x()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetri2x/zhetri2x
// *HETRI2X computes the inverse of a complex Hermitian indefinite matrix
// using the blocked workspace kernel behind *HETRI2.

chetri2x :: #bake_arguments _hetri2x(RT=float32, MT=C32);
zhetri2x :: #bake_arguments _hetri2x(RT=float64, MT=C64);

#scope_file
_hetri2x :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D, $$nb: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    routine_name :: #run if MT == C32 { "CHETRI2X" } else { "ZHETRI2X" };
    upper := lsame(uplo, #char "U");

    info.* = 0;
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if nb < 1 {
        info.* = -7;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 {
        return;
    }

    // Staged port note:
    // The blocked *HETRI2X kernel still needs a literal follow-up pass.
    // This file establishes the family interface and validation shape.
}
