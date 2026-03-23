// ? LAPACK *hegvx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chegvx/zhegvx
// *HEGVX computes selected eigenvalues, and optionally, eigenvectors
// of a complex generalized Hermitian-definite eigenproblem.

chegvx :: (
    $$itype: int, $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    vl: float32, vu: float32, $$il: int, $$iu: int, abstol: float32, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D,
    $$lwork: int, rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D,
    ifail: *$IFAIL/ArrayType1D, info: *int
) {
    _hegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info, RT=float32, MT=C32);
}

zhegvx :: (
    $$itype: int, $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    vl: float64, vu: float64, $$il: int, $$iu: int, abstol: float64, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D,
    $$lwork: int, rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D,
    ifail: *$IFAIL/ArrayType1D, info: *int
) {
    _hegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info, RT=float64, MT=C64);
}

#scope_file
_hegvx_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hegvx_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return x;
    }
}

_hegvx_potrf :: (
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 {
        cpotrf(uplo, n, a, lda, info);
    } else {
        zpotrf(uplo, n, a, lda, info);
    }
}

_hegvx_hegst :: (
    $$itype: int,
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == C32 {
        chegst(itype, uplo, n, a, lda, b, ldb, info);
    } else {
        zhegst(itype, uplo, n, a, lda, b, ldb, info);
    }
}

_hegvx_heevx :: (
    $$jobz: u8,
    $$range: u8,
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    vl: $RT,
    vu: RT,
    $$il: int,
    $$iu: int,
    abstol: RT,
    m: *int,
    w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    ifail: *$IFAIL/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           W.data_type == RT &&
           Z.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int &&
           IFAIL.data_type == int;
} {
    #if MT == C32 {
        cheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
    } else {
        zheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
    }
}

_hegvx_trmm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == C32 {
        ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_hegvx_trsm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_hegvx :: (
    $$itype: int, $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    vl: $RT, vu: RT, $$il: int, $$iu: int, abstol: RT, m: *int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D,
    $$lwork: int, rwork: *$RWORK/ArrayType1D, iwork: *$IWORK/ArrayType1D,
    ifail: *$IFAIL/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT &&
           W.data_type == RT &&
           Z.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int &&
           IFAIL.data_type == int;
} {
    cone_t :: #run _hegvx_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHEGVX" } else { "ZHEGVX" };

    alleig, indeig, lquery, upper, valeig, wantz : bool = ---;
    lwkopt, nb : int = ---;
    trans : u8 = ---;

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");
    alleig = lsame(range, #char "A");
    valeig = lsame(range, #char "V");
    indeig = lsame(range, #char "I");
    lquery = lwork == -1;

    info.* = 0;
    if itype < 1 || itype > 3 {
        info.* = -1;
    } else if !(wantz || lsame(jobz, #char "N")) {
        info.* = -2;
    } else if !(alleig || valeig || indeig) {
        info.* = -3;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -4;
    } else if n < 0 {
        info.* = -5;
    } else if lda < max(1, n) {
        info.* = -7;
    } else if ldb < max(1, n) {
        info.* = -9;
    } else {
        if valeig {
            if n > 0 && vu <= vl {
                info.* = -11;
            }
        } else if indeig {
            if il < 1 || il > max(1, n) {
                info.* = -12;
            } else if iu < min(n, il) || iu > n {
                info.* = -13;
            }
        }
    }
    if info.* == 0 {
        if ldz < 1 || (wantz && ldz < n) {
            info.* = -18;
        }
    }

    if info.* == 0 {
        nb = ilaenv(1, ifx MT == C32 then "CHETRD" else "ZHETRD", cast(*u8) *uplo, n, -1, -1, -1);
        lwkopt = max(1, (nb + 1)*n);
        fset(work, 1, _hegvx_complex_from_real(_hegvx_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
        if lwork < max(1, 2*n) && !lquery {
            info.* = -20;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    m.* = 0;
    if n == 0 then return;

    _hegvx_potrf(uplo, n, b, ldb, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    _hegvx_hegst(itype, uplo, n, a, lda, b, ldb, info, MT);
    _hegvx_heevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info, RT, MT);

    if wantz {
        if info.* > 0 {
            m.* = info.* - 1;
        }

        if itype == 1 || itype == 2 {
            if upper {
                trans = #char "N";
            } else {
                trans = #char "C";
            }
            _hegvx_trsm(#char "L", uplo, trans, #char "N", n, m.*, cone_t, b.*, ldb, z, ldz, MT);
        } else if itype == 3 {
            if upper {
                trans = #char "C";
            } else {
                trans = #char "N";
            }
            _hegvx_trmm(#char "L", uplo, trans, #char "N", n, m.*, cone_t, b.*, ldb, z, ldz, MT);
        }
    }

    fset(work, 1, _hegvx_complex_from_real(_hegvx_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
}
