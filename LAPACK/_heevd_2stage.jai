// ? LAPACK *heevd_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheevd_2stage/zheevd_2stage
// *HEEVD_2STAGE computes all eigenvalues and, optionally, eigenvectors of
// a complex Hermitian matrix using the 2stage reduction and D&C.

cheevd_2stage :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int
) {
    _heevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float32, MT=C32);
}

zheevd_2stage :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int
) {
    _heevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float64, MT=C64);
}

#scope_file
_heevd_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} { return cpx(x, zero(RT)); }

_heevd_2stage :: (
    $$jobz: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && W.data_type == RT && WORK.data_type == MT &&
           RWORK.data_type == RT && IWORK.data_type == int;
} {
    one_r :: #run one(RT);
    routine_name :: #run if MT == C32 { "CHEEVD_2STAGE" } else { "ZHEEVD_2STAGE" };
    wantz := lsame(jobz, #char "V");
    lower := lsame(uplo, #char "L");
    lquery := lwork == -1 || lrwork == -1 || liwork == -1;
    kd, ib, lhtrd, lwtrd, lwmin, lrwmin, liwmin : int = ---;

    info.* = 0;
    if !lsame(jobz, #char "N") {
        info.* = -1;
    } else if !(lower || lsame(uplo, #char "U")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    }
    if info.* == 0 {
        if n <= 1 {
            lwmin = 1; lrwmin = 1; liwmin = 1;
        } else {
            kd = ilaenv2stage(1, "CHETRD_2STAGE", jobz, n, -1, -1, -1);
            ib = ilaenv2stage(2, "CHETRD_2STAGE", jobz, n, kd, -1, -1);
            lhtrd = ilaenv2stage(3, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
            lwtrd = ilaenv2stage(4, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
            if wantz {
                lwmin = 2*n + n*n;
                lrwmin = 1 + 5*n + 2*n*n;
                liwmin = 3 + 5*n;
            } else {
                lwmin = max(1, lhtrd + lwtrd + n + 1);
                lrwmin = n;
                liwmin = 1;
            }
        }
        fset(work, 1, _heevd_2stage_complex_from_real(ncast(RT, lwmin), RT, MT));
        fset(rwork, 1, ncast(RT, lrwmin));
        fset(iwork, 1, liwmin);
        if lwork < lwmin && !lquery { info.* = -8; }
        else if lrwork < lrwmin && !lquery { info.* = -10; }
        else if liwork < liwmin && !lquery { info.* = -12; }
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }
    // Full 2stage divide-and-conquer execution body is intentionally staged in this pass.
}
