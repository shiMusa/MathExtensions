// ? LAPACK *hegv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chegv/zhegv
// *HEGV computes all the eigenvalues and, optionally, the eigenvectors
// of a complex generalized Hermitian-definite eigenproblem.

chegv :: #bake_arguments _hegv(RT=float32, MT=C32);
zhegv :: #bake_arguments _hegv(RT=float64, MT=C64);

#scope_file
_hegv_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hegv_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return x;
    }
}

_hegv_potrf :: (
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

_hegv_hegst :: (
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

_hegv_heev :: (
    $$jobz: u8,
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           W.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    #if MT == C32 {
        cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
    } else {
        zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
    }
}

_hegv_trmm :: (
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

_hegv_trsm :: (
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

_hegv :: (
    $$itype: int,
    $$jobz: u8,
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT &&
           W.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    one_t :: #run _hegv_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHEGV " } else { "ZHEGV " };

    lquery, upper, wantz : bool = ---;
    lwkopt, nb, neig : int = ---;
    trans : u8 = ---;

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;

    info.* = 0;
    if itype < 1 || itype > 3 {
        info.* = -1;
    } else if !(wantz || lsame(jobz, #char "N")) {
        info.* = -2;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if lda < max(1, n) {
        info.* = -6;
    } else if ldb < max(1, n) {
        info.* = -8;
    }

    if info.* == 0 {
        nb = ilaenv(1, ifx MT == C32 then "CHETRD" else "ZHETRD", cast(*u8) *uplo, n, -1, -1, -1);
        lwkopt = max(1, (nb + 1)*n);
        fset(work, 1, _hegv_complex_from_real(_hegv_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
        if lwork < max(1, 2*n - 1) && !lquery {
            info.* = -11;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 then return;

    _hegv_potrf(uplo, n, b, ldb, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    _hegv_hegst(itype, uplo, n, a, lda, b, ldb, info, MT);
    _hegv_heev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info, RT, MT);

    if wantz {
        neig = n;
        if info.* > 0 {
            neig = info.* - 1;
        }

        if itype == 1 || itype == 2 {
            if upper {
                trans = #char "N";
            } else {
                trans = #char "C";
            }
            _hegv_trsm(#char "L", uplo, trans, #char "N", n, neig, one_t, b.*, ldb, a, lda, MT);
        } else if itype == 3 {
            if upper {
                trans = #char "C";
            } else {
                trans = #char "N";
            }
            _hegv_trmm(#char "L", uplo, trans, #char "N", n, neig, one_t, b.*, ldb, a, lda, MT);
        }
    }

    fset(work, 1, _hegv_complex_from_real(_hegv_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
}
