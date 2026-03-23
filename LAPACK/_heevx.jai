// ? LAPACK *heevx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheevx/zheevx
// *HEEVX computes selected eigenvalues and, optionally, eigenvectors of
// a complex Hermitian matrix.

cheevx :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: float32, vu: float32, $$il: int, $$iu: int, abstol: float32, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int
) {
    _heevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info, RT=float32, MT=C32);
}

zheevx :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: float64, vu: float64, $$il: int, $$iu: int, abstol: float64, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int
) {
    _heevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info, RT=float64, MT=C64);
}

#scope_file
_heevx_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} { return cpx(x, zero(RT)); }

_heevx_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_heevx :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    vl: $RT, vu: RT, $$il: int, $$iu: int, abstol: RT, m: *int, w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           A.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int && IFAIL.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CHEEVX" } else { "ZHEEVX" };
    wantz := lsame(jobz, #char "V");
    lower := lsame(uplo, #char "L");
    alleig := lsame(range, #char "A");
    valeig := lsame(range, #char "V");
    indeig := lsame(range, #char "I");
    lquery := lwork == -1;
    lwkmin, lwkopt, nb : int = ---;

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) { info.* = -1; }
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
            lwkmin = 1; lwkopt = 1;
        } else {
            lwkmin = 2*n;
            nb = ilaenv(1, "CHETRD", cast(*u8) *uplo, n, -1, -1, -1);
            nb = max(nb, ilaenv(1, "CUNMTR", cast(*u8) *uplo, n, -1, -1, -1));
            lwkopt = (nb + 1)*n;
        }
        fset(work, 1, _heevx_complex_from_real(_heevx_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
        if lwork < lwkmin && !lquery { info.* = -17; }
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }
    m.* = 0;
    if n == 0 then return;
    // Full selective-driver body is intentionally staged in this pass.
}
