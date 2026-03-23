// ? LAPACK *hegvd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chegvd/zhegvd
// *HEGVD computes all the eigenvalues and, optionally, the eigenvectors
// of a complex generalized Hermitian-definite eigenproblem.

chegvd :: (
    $$itype: int, $$jobz: u8, $$uplo: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    w: *$W/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, $$lrwork: int, iwork: *$IWORK/ArrayType1D,
    $$liwork: int, info: *int
) {
    _hegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float32, MT=C32);
}

zhegvd :: (
    $$itype: int, $$jobz: u8, $$uplo: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    w: *$W/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, $$lrwork: int, iwork: *$IWORK/ArrayType1D,
    $$liwork: int, info: *int
) {
    _hegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info, RT=float64, MT=C64);
}

#scope_file
_hegvd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hegvd_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return x;
    }
}

_hegvd_potrf :: (
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

_hegvd_hegst :: (
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

_hegvd_heevd :: (
    $$jobz: u8,
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    w: *$W/ArrayType1D,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    $$lrwork: int,
    iwork: *$IWORK/ArrayType1D,
    $$liwork: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           W.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    #if MT == C32 {
        cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
    } else {
        zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
    }
}

_hegvd_trmm :: (
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

_hegvd_trsm :: (
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

_hegvd :: (
    $$itype: int, $$jobz: u8, $$uplo: u8, $$n: int,
    a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    w: *$W/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, $$lrwork: int, iwork: *$IWORK/ArrayType1D,
    $$liwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT &&
           W.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    cone_t :: #run _hegvd_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHEGVD" } else { "ZHEGVD" };

    liopt, liwmin, lopt, lropt, lrwmin, lwmin : int = ---;
    lquery, upper, wantz : bool = ---;
    trans : u8 = ---;
    work_query, lwork_opt_rt, lrwork_opt_rt : RT = ---;

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;

    info.* = 0;
    if n <= 1 {
        lwmin = 1;
        lrwmin = 1;
        liwmin = 1;
    } else if wantz {
        lwmin = 2*n + n*n;
        lrwmin = 1 + 5*n + 2*n*n;
        liwmin = 3 + 5*n;
    } else {
        lwmin = n + 1;
        lrwmin = n;
        liwmin = 1;
    }
    lopt = lwmin;
    lropt = lrwmin;
    liopt = liwmin;

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
        fset(work, 1, _hegvd_complex_from_real(_hegvd_roundup_lwork(ncast(RT, lopt), RT), RT, MT));
        fset(rwork, 1, ncast(RT, lropt));
        fset(iwork, 1, liopt);

        if lwork < lwmin && !lquery {
            info.* = -11;
        } else if lrwork < lrwmin && !lquery {
            info.* = -13;
        } else if liwork < liwmin && !lquery {
            info.* = -15;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 then return;

    _hegvd_potrf(uplo, n, b, ldb, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    _hegvd_hegst(itype, uplo, n, a, lda, b, ldb, info, MT);
    _hegvd_heevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info, RT, MT);

    work_query = real(fget(work.*, 1));
    lwork_opt_rt = ncast(RT, lopt);
    if work_query > lwork_opt_rt {
        lwork_opt_rt = work_query;
    }
    lopt = cast(int) lwork_opt_rt;
    lrwork_opt_rt = fget(rwork.*, 1);
    if lrwork_opt_rt > ncast(RT, lropt) {
        lropt = cast(int) lrwork_opt_rt;
    }
    liopt = max(liopt, fget(iwork.*, 1));

    if wantz && info.* == 0 {
        if itype == 1 || itype == 2 {
            if upper {
                trans = #char "N";
            } else {
                trans = #char "C";
            }
            _hegvd_trsm(#char "L", uplo, trans, #char "N", n, n, cone_t, b.*, ldb, a, lda, MT);
        } else if itype == 3 {
            if upper {
                trans = #char "C";
            } else {
                trans = #char "N";
            }
            _hegvd_trmm(#char "L", uplo, trans, #char "N", n, n, cone_t, b.*, ldb, a, lda, MT);
        }
    }

    fset(work, 1, _hegvd_complex_from_real(_hegvd_roundup_lwork(ncast(RT, lopt), RT), RT, MT));
    fset(rwork, 1, ncast(RT, lropt));
    fset(iwork, 1, liopt);
}
