// ? LAPACK *lahef_aa()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clahef_aa/zlahef_aa
// *LAHEF_AA factorizes a panel of a complex Hermitian matrix
// using Aasen's algorithm.

clahef_aa :: #bake_arguments _lahef_aa(MT=C32);
zlahef_aa :: #bake_arguments _lahef_aa(MT=C64);

#scope_file
_lahef_aa :: (
    $$uplo: u8,
    $$j1: int,
    $$m: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    h: *$H/ArrayType2D,
    $$ldh: int,
    work: *$WORK/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int &&
           H.data_type == MT &&
           WORK.data_type == MT;
} {
    panel_n := min(m, nb);
    for i:1..max(0, panel_n) {
        fset(ipiv, i, i);
    }

    // Staged port note:
    // The Aasen panel kernel from *lahef_aa.f is still pending.
    _ = uplo; _ = j1; _ = a; _ = lda; _ = h; _ = ldh; _ = work;
}
