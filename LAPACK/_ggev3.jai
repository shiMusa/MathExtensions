// ? LAPACK *ggev3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sggev3/cggev3/dggev3/zggev3
// *GGEV3 computes the generalized eigenvalues and, optionally, the
// left and/or right generalized eigenvectors of a matrix pair using
// the blocked generalized Hessenberg/QZ path.

sggev3 :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D,
    alphai: *$ALPHAI/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int
) {
    _ggev_real(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info, blocked=true, RT=float32);
}

dggev3 :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alphar: *$ALPHAR/ArrayType1D,
    alphai: *$ALPHAI/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    info: *int
) {
    _ggev_real(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info, blocked=true, RT=float64);
}

cggev3 :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alpha: *$ALPHA/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _ggev_complex(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info, blocked=true, RT=float32, MT=C32);
}

zggev3 :: (
    $$jobvl: u8,
    $$jobvr: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    alpha: *$ALPHA/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vl: *$VL/ArrayType2D,
    $$ldvl: int,
    vr: *$VR/ArrayType2D,
    $$ldvr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _ggev_complex(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info, blocked=true, RT=float64, MT=C64);
}
