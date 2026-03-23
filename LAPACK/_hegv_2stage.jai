// ? LAPACK *hegv_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chegv_2stage/zhegv_2stage
// *HEGV_2STAGE computes all eigenvalues and, optionally, the eigenvectors
// of a complex generalized Hermitian-definite eigenproblem.

chegv_2stage :: #bake_arguments _hegv_2stage(RT=float32, MT=C32);
zhegv_2stage :: #bake_arguments _hegv_2stage(RT=float64, MT=C64);

#scope_file
_hegv_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hegv_2stage_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return x;
    }
}

_hegv_2stage_potrf :: (
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

_hegv_2stage_hegst :: (
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

_hegv_2stage_heev_2stage :: (
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
        cheev_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
    } else {
        zheev_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
    }
}

_hegv_2stage_trmm :: (
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

_hegv_2stage_trsm :: (
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

_hegv_2stage :: (
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
    one_t :: #run _hegv_2stage_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHEGV_2STAGE " } else { "ZHEGV_2STAGE " };

    ib, kd, lhtrd, lwmin, lwtrd, neig : int = ---;
    lquery, upper, wantz : bool = ---;
    trans : u8 = ---;

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;

    info.* = 0;
    if itype < 1 || itype > 3 {
        info.* = -1;
    } else if !lsame(jobz, #char "N") {
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
        kd = ilaenv2stage(1, ifx MT == C32 then "CHETRD_2STAGE" else "ZHETRD_2STAGE", jobz, n, -1, -1, -1);
        ib = ilaenv2stage(2, ifx MT == C32 then "CHETRD_2STAGE" else "ZHETRD_2STAGE", jobz, n, kd, -1, -1);
        lhtrd = ilaenv2stage(3, ifx MT == C32 then "CHETRD_2STAGE" else "ZHETRD_2STAGE", jobz, n, kd, ib, -1);
        lwtrd = ilaenv2stage(4, ifx MT == C32 then "CHETRD_2STAGE" else "ZHETRD_2STAGE", jobz, n, kd, ib, -1);
        lwmin = n + lhtrd + lwtrd;
        fset(work, 1, _hegv_2stage_complex_from_real(_hegv_2stage_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
        if lwork < lwmin && !lquery {
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

    _hegv_2stage_potrf(uplo, n, b, ldb, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    _hegv_2stage_hegst(itype, uplo, n, a, lda, b, ldb, info, MT);
    _hegv_2stage_heev_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, info, RT, MT);

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
            _hegv_2stage_trsm(#char "L", uplo, trans, #char "N", n, neig, one_t, b.*, ldb, a, lda, MT);
        } else if itype == 3 {
            if upper {
                trans = #char "C";
            } else {
                trans = #char "N";
            }
            _hegv_2stage_trmm(#char "L", uplo, trans, #char "N", n, neig, one_t, b.*, ldb, a, lda, MT);
        }
    }

    fset(work, 1, _hegv_2stage_complex_from_real(_hegv_2stage_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
}
