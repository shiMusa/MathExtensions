// ? LAPACK *hbgvd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chbgvd/zhbgvd
// *HBGVD computes all eigenvalues and, optionally, eigenvectors of a
// complex generalized Hermitian-definite banded eigenproblem using
// divide and conquer.

chbgvd :: #bake_arguments _hbgvd(RT=float32, MT=C32);
zhbgvd :: #bake_arguments _hbgvd(RT=float64, MT=C64);

#scope_file
_hbgvd_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hbgvd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hbgvd_query_value :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        return cpx(sroundup_lwork(x), zero(RT));
    } else {
        return cpx(x, zero(RT));
    }
}

_hbgvd_pbstf :: (
    $$uplo: u8, $$n: int, $$kb: int, bb: *$BB/ArrayType2D, $$ldbb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && BB.data_type == MT;
} {
    #if MT == C32 { cpbstf(uplo, n, kb, bb, ldbb, info); } else { zpbstf(uplo, n, kb, bb, ldbb, info); }
}

_hbgvd_hbgst :: (
    $$vect: u8, $$uplo: u8, $$n: int, $$ka: int, $$kb: int,
    ab: *$AB/ArrayType2D, $$ldab: int, bb: *$BB/ArrayType2D, $$ldbb: int,
    x: *$X/ArrayType2D, $$ldx: int, work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AB.data_type == MT && BB.data_type == MT && X.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { chbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, rwork, info); }
    else { zhbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, rwork, info); }
}

_hbgvd_hbtrd :: (
    $$vect: u8, $$uplo: u8, $$n: int, $$kd: int, ab: *$AB/ArrayType2D, $$ldab: int,
    d: *$D/ArrayType1D, e: *$E/ArrayType1D, q: *$Q/ArrayType2D, $$ldq: int,
    work: *$WORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AB.data_type == MT && D.data_type == RT && E.data_type == RT &&
           Q.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 { chbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info); }
    else { zhbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info); }
}

_hbgvd_stedc :: (
    $$compz: u8, $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int;
} {
    #if MT == C32 { cstedc(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info); }
    else { zstedc(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info); }
}

_hbgvd_sterf :: (
    $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && D.data_type == RT && E.data_type == RT;
} {
    #if RT == float32 { ssterf(n, d, e, info); } else { dsterf(n, d, e, info); }
}

_hbgvd_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int, alpha: $MT,
    a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int, beta: MT,
    c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 { cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }
    else { zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }
}

_hbgvd_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 { clacpy(uplo, m, n, a, lda, b, ldb); } else { zlacpy(uplo, m, n, a, lda, b, ldb); }
}

_hbgvd :: (
    $$jobz: u8, $$uplo: u8, $$n: int, $$ka: int, $$kb: int,
    ab: *$AB/ArrayType2D, $$ldab: int, bb: *$BB/ArrayType2D, $$ldbb: int,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, $$lrwork: int,
    iwork: *$IWORK/ArrayType1D, $$liwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AB.data_type == MT && BB.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int;
} {
    zero_t :: #run _hbgvd_complex_from_real(zero(RT), RT, MT);
    cone_t :: #run _hbgvd_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHBGVD" } else { "ZHBGVD" };

    iinfo, inde, indwk2, indwrk, liwmin, llrwk, llwk2, lrwmin, lwmin : int = ---;
    lquery, upper, wantz : bool = ---;
    vect : u8 = ---;

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;

    info.* = 0;
    if n <= 1 {
        lwmin = 1 + n;
        lrwmin = 1 + n;
        liwmin = 1;
    } else if wantz {
        lwmin = 2*n*n;
        lrwmin = 1 + 5*n + 2*n*n;
        liwmin = 3 + 5*n;
    } else {
        lwmin = n;
        lrwmin = n;
        liwmin = 1;
    }

    if !(wantz || lsame(jobz, #char "N")) {
        info.* = -1;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if ka < 0 {
        info.* = -4;
    } else if kb < 0 || kb > ka {
        info.* = -5;
    } else if ldab < ka + 1 {
        info.* = -7;
    } else if ldbb < kb + 1 {
        info.* = -9;
    } else if ldz < 1 || (wantz && ldz < n) {
        info.* = -12;
    }

    if info.* == 0 {
        fset(work, 1, _hbgvd_query_value(ncast(RT, lwmin), RT, MT));
        fset(rwork, 1, ncast(RT, lrwmin));
        fset(iwork, 1, liwmin);

        if lwork < lwmin && !lquery {
            info.* = -14;
        } else if lrwork < lrwmin && !lquery {
            info.* = -16;
        } else if liwork < liwmin && !lquery {
            info.* = -18;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        return;
    }

    _hbgvd_pbstf(uplo, n, kb, bb, ldbb, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    inde = 1;
    indwrk = inde + n;
    indwk2 = 1 + n*n;
    llwk2 = lwork - indwk2 + 2;
    llrwk = lrwork - indwrk + 2;
    _hbgvd_hbgst(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, z, ldz, work, rwork, *iinfo, RT, MT);

    if wantz { vect = #char "U"; } else { vect = #char "N"; }
    e_view := _hbgvd_subvector(rwork.*, inde, n);
    _hbgvd_hbtrd(vect, uplo, n, ka, ab, ldab, w, *e_view, z, ldz, work, *iinfo, RT, MT);

    if !wantz {
        _hbgvd_sterf(n, w, *e_view, info, RT);
    } else {
        work_z := array2d_view(work.*, n, row_major=false);
        work2_view := _hbgvd_subvector(work.*, indwk2, llwk2);
        rwork2_view := _hbgvd_subvector(rwork.*, indwrk, llrwk);
        work2_mat_view := _hbgvd_subvector(work.*, indwk2, n*n);
        work2_mat := array2d_view(work2_mat_view, n, row_major=false);

        _hbgvd_stedc(#char "I", n, w, *e_view, *work_z, n, *work2_view, llwk2, *rwork2_view, llrwk, iwork, liwork, info, RT, MT);
        _hbgvd_gemm(#char "N", #char "N", n, n, n, cone_t, z.*, ldz, *work_z, n, zero_t, *work2_mat, n, MT);
        _hbgvd_lacpy(#char "A", n, n, *work2_mat, n, z, ldz, MT);
    }

    fset(work, 1, _hbgvd_query_value(ncast(RT, lwmin), RT, MT));
    fset(rwork, 1, ncast(RT, lrwmin));
    fset(iwork, 1, liwmin);
}
