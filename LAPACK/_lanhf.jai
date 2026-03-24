// ? LAPACK *lanhf()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clanhf/zlanhf
// *LANHF returns the value of the 1-norm, Frobenius norm,
// infinity-norm, or the largest absolute value of any element of a
// Hermitian matrix in RFP format.

clanhf :: #bake_arguments _lanhf(RT=float32, MT=C32);
zlanhf :: #bake_arguments _lanhf(RT=float64, MT=C64);

#scope_file
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

    if n <= 0 {
        return zero_r;
    }

    // Staged port note:
    // *LANHF is much larger than the surrounding norm helpers because
    // it expands all RFP-format cases across TRANSR/UPLO and odd/even N.
    // This file lands the public c/z family interface that current
    // callers expect; the full RFP case split is still pending.
    _ = norm; _ = transr; _ = uplo; _ = a; _ = work;
    return zero_r;
}
