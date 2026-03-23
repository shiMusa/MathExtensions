// ? LAPACK *heevx_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheevx_2stage/zheevx_2stage
// *HEEVX_2STAGE computes selected eigenvalues and, optionally,
// eigenvectors of a complex Hermitian matrix using the 2stage reduction.

cheevx_2stage :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: float32, vu: float32, $$il: int, $$iu: int, abstol: float32, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int
) {
    _heevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info, RT=float32, MT=C32);
}

zheevx_2stage :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: float64, vu: float64, $$il: int, $$iu: int, abstol: float64, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int
) {
    _heevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info, RT=float64, MT=C64);
}

#scope_file
_heevx_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} { return cpx(x, zero(RT)); }

_heevx_2stage :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: $RT, vu: RT, $$il: int, $$iu: int, abstol: RT, m: *int, w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int && IFAIL.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CHEEVX_2STAGE" } else { "ZHEEVX_2STAGE" };
    wantz := lsame(jobz, #char "V");
    lower := lsame(uplo, #char "L");
    alleig := lsame(range, #char "A");
    valeig := lsame(range, #char "V");
    indeig := lsame(range, #char "I");
    lquery := lwork == -1;
    kd, ib, lhtrd, lwtrd, lwmin : int = ---;

    info.* = 0;
    if !lsame(jobz, #char "N") { info.* = -1; }
    else if !(alleig || valeig || indeig) { info.* = -2; }
    else if !(lower || lsame(uplo, #char "U")) { info.* = -3; }
    else if n < 0 { info.* = -4; }
    else if lda < max(1, n) { info.* = -6; }
    else {
        if valeig {
            if n > 0 && vu <= vl { info.* = -8; }
        } else if indeig {
            if il < 1 || il > max(1, n) { info.* = -9; }
            else if iu < min(n, il) || iu > n { info.* = -10; }
        }
    }
    if info.* == 0 && (ldz < 1 || (wantz && ldz < n)) { info.* = -15; }
    if info.* == 0 {
        if n <= 1 {
            lwmin = 1;
        } else {
            kd = ilaenv2stage(1, "CHETRD_2STAGE", jobz, n, -1, -1, -1);
            ib = ilaenv2stage(2, "CHETRD_2STAGE", jobz, n, kd, -1, -1);
            lhtrd = ilaenv2stage(3, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
            lwtrd = ilaenv2stage(4, "CHETRD_2STAGE", jobz, n, kd, ib, -1);
            lwmin = ifx wantz then 0 else max(1, max(8*n, lhtrd + lwtrd + (kd + 1)*n + n));
        }
        fset(work, 1, _heevx_2stage_complex_from_real(ncast(RT, lwmin), RT, MT));
        if lwork < lwmin && !lquery { info.* = -17; }
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }
    m.* = 0;
    // Full selective 2stage body is intentionally staged in this pass.
}
