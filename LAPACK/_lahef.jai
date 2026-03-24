// ? LAPACK *lahef()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clahef/zlahef
// *LAHEF computes a partial factorization of a complex Hermitian
// indefinite matrix using the blocked Bunch-Kaufman method.

clahef :: #bake_arguments _lahef(MT=C32);
zlahef :: #bake_arguments _lahef(MT=C64);

#scope_file
_lahef :: (
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
    // The blocked Bunch-Kaufman panel kernel from *lahef.f is still pending.
    _ = uplo; _ = a; _ = lda; _ = w; _ = ldw;
}
