// ? LAPACK *gges3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgges3/cgges3/dgges3/zgges3
// *GGES3 computes the generalized eigenvalues, generalized Schur form,
// and optionally the left and/or right Schur vectors of a matrix pair,
// using the blocked generalized Hessenberg/QZ path.

sgges3 :: (
    $$jobvsl: u8,
    $$jobvsr: u8,
    $$sort: u8,
    selctg: (float32, float32, float32) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    sdim: *int,
    alphar: *$ALPHAR/ArrayType1D,
    alphai: *$ALPHAI/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D,
    $$ldvsl: int,
    vsr: *$VSR/ArrayType2D,
    $$ldvsr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _gges_real(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info, blocked=true, RT=float32);
}

dgges3 :: (
    $$jobvsl: u8,
    $$jobvsr: u8,
    $$sort: u8,
    selctg: (float64, float64, float64) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    sdim: *int,
    alphar: *$ALPHAR/ArrayType1D,
    alphai: *$ALPHAI/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D,
    $$ldvsl: int,
    vsr: *$VSR/ArrayType2D,
    $$ldvsr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _gges_real(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info, blocked=true, RT=float64);
}

cgges3 :: (
    $$jobvsl: u8,
    $$jobvsr: u8,
    $$sort: u8,
    selctg: (C32, C32) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    sdim: *int,
    alpha: *$ALPHA/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D,
    $$ldvsl: int,
    vsr: *$VSR/ArrayType2D,
    $$ldvsr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _gges_complex(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info, blocked=true, RT=float32, MT=C32);
}

zgges3 :: (
    $$jobvsl: u8,
    $$jobvsr: u8,
    $$sort: u8,
    selctg: (C64, C64) -> bool,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    sdim: *int,
    alpha: *$ALPHA/ArrayType1D,
    beta: *$BETA/ArrayType1D,
    vsl: *$VSL/ArrayType2D,
    $$ldvsl: int,
    vsr: *$VSR/ArrayType2D,
    $$ldvsr: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    bwork: *$BWORK/ArrayType1D,
    info: *int
) {
    _gges_complex(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info, blocked=true, RT=float64, MT=C64);
}
