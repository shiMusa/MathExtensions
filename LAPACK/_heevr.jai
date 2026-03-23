// ? LAPACK *heevr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheevr/zheevr
// *HEEVR computes selected eigenvalues and, optionally, eigenvectors of
// a complex Hermitian matrix using the MRRR / fallback path.

cheevr :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: float32, vu: float32, $$il: int, $$iu: int, abstol: float32, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, isuppz: *$ISUPPZ/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int
) {
    _heevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float32, MT=C32);
}

zheevr :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: float64, vu: float64, $$il: int, $$iu: int, abstol: float64, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, isuppz: *$ISUPPZ/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int
) {
    _heevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float64, MT=C64);
}

#scope_file
_heevr_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} { return cpx(x, zero(RT)); }

_heevr :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: $RT, vu: RT, $$il: int, $$iu: int, abstol: RT, m: *int, w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D, $$ldz: int, isuppz: *$ISUPPZ/ArrayType1D, work: *$WORK/ArrayType1D,
    $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int, iwork: *$IWORK/ArrayType1D,
    $$liwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && W.data_type == RT && Z.data_type == MT && ISUPPZ.data_type == int &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CHEEVR" } else { "ZHEEVR" };
    lower := lsame(uplo, #char "L");
    wantz := lsame(jobz, #char "V");
    alleig := lsame(range, #char "A");
    valeig := lsame(range, #char "V");
    indeig := lsame(range, #char "I");
    lquery := lwork == -1 || lrwork == -1 || liwork == -1;

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) { info.* = -1; }
    else if !(alleig || valeig || indeig) { info.* = -2; }
    else if !(lower || lsame(uplo, #char "U")) { info.* = -3; }
    else if n < 0 { info.* = -4; }
    else if lda < max(1, n) { info.* = -6; }
    else if ldz < 1 || (wantz && ldz < n) { info.* = -15; }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }
    if lquery {
        fset(work, 1, _heevr_complex_from_real(ncast(RT, max(1, 2*n)), RT, MT));
        fset(rwork, 1, ncast(RT, max(1, 24*n)));
        fset(iwork, 1, max(1, 10*n));
        return;
    }
    m.* = 0;
    // Full MRRR / SSTEBZ+CSTEIN driver body is intentionally staged in this pass.
}
