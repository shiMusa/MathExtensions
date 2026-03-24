// ? LAPACK *lahef_rk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clahef_rk/zlahef_rk
// *LAHEF_RK computes a partial factorization of a complex Hermitian
// indefinite matrix using the bounded Bunch-Kaufman rook method.

clahef_rk :: #bake_arguments _lahef_rk(MT=C32);
zlahef_rk :: #bake_arguments _lahef_rk(MT=C64);

#scope_file
_lahef_rk :: (
    $$uplo: u8,
    $$n: int,
    $$nb: int,
    kb: *int,
    a: *$A/ArrayType2D,
    $$lda: int,
    e: *$E/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    w: *$W/ArrayType2D,
    $$ldw: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           E.data_type == MT &&
           IPIV.data_type == int &&
           W.data_type == MT;
} {
    zero_t :: #run zero(MT);

    info.* = 0;
    kb.* = min(n, nb);

    for i:1..n {
        fset(ipiv, i, i);
        fset(e, i, zero_t);
    }

    // Staged port note:
    // The blocked rook panel kernel from *lahef_rk.f is still pending.
    _ = uplo; _ = a; _ = lda; _ = w; _ = ldw;
}
