// ? LAPACK *lahef_rook()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clahef_rook/zlahef_rook
// *LAHEF_ROOK computes a partial factorization of a complex Hermitian
// indefinite matrix using the blocked rook pivoting method.

clahef_rook :: #bake_arguments _lahef_rook(MT=C32);
zlahef_rook :: #bake_arguments _lahef_rook(MT=C64);

#scope_file
_lahef_rook :: (
    $$uplo: u8,
    $$n: int,
    $$nb: int,
    kb: *int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    w: *$W/ArrayType2D,
    $$ldw: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int &&
           W.data_type == MT;
} {
    info.* = 0;
    kb.* = min(n, nb);

    for i:1..n {
        fset(ipiv, i, i);
    }

    // Staged port note:
    // The blocked rook panel kernel from *lahef_rook.f is still pending.
    _ = uplo; _ = a; _ = lda; _ = w; _ = ldw;
}
