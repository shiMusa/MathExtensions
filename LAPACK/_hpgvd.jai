// ? LAPACK *hpgvd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chpgvd/zhpgvd
// *HPGVD computes all eigenvalues and, optionally, eigenvectors of a
// complex generalized Hermitian-definite eigenproblem in packed storage.

chpgvd :: #bake_arguments _hpgvd(RT=float32, MT=C32);
zhpgvd :: #bake_arguments _hpgvd(RT=float64, MT=C64);

#scope_file
_hpgvd_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld;
    view.count = count;
    return view;
}

_hpgvd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hpgvd_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hpgvd_pptrf :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT;
} {
    #if MT == C32 { cpptrf(uplo, n, ap, info); } else { zpptrf(uplo, n, ap, info); }
}

_hpgvd_hpevd :: (
    $$jobz: u8, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, $$lrwork: int, iwork: *$IWORK/ArrayType1D, $$liwork: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int;
} {
    #if MT == C32 { chpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info); }
    else { zhpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info); }
}

_hpgvd_tpmv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int, ap: $AP/ArrayType1D, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 { ctpmv(uplo, trans, diag, n, ap, x, incx); } else { ztpmv(uplo, trans, diag, n, ap, x, incx); }
}

_hpgvd_tpsv :: (
    $$uplo: u8, $$trans: u8, $$diag: u8, $$n: int, ap: $AP/ArrayType1D, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT;
} {
    #if MT == C32 { ctpsv(uplo, trans, diag, n, ap, x, incx); } else { ztpsv(uplo, trans, diag, n, ap, x, incx); }
}

_hpgvd :: (
    $$itype: int, $$jobz: u8, $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, bp: *$BP/ArrayType1D,
    w: *$W/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, $$lwork: int,
    rwork: *$RWORK/ArrayType1D, $$lrwork: int, iwork: *$IWORK/ArrayType1D, $$liwork: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && BP.data_type == MT && W.data_type == RT && Z.data_type == MT &&
           WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int;
} {
    lquery, upper, wantz : bool = ---;
    j, liwmin, lrwmin, lwmin, neig : int = ---;
    trans : u8 = ---;
    routine_name :: #run if MT == C32 { "CHPGVD" } else { "ZHPGVD" };

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;

    info.* = 0;
    if itype < 1 || itype > 3 {
        info.* = -1;
    } else if !(wantz || lsame(jobz, #char "N")) {
        info.* = -2;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if ldz < 1 || (wantz && ldz < n) {
        info.* = -9;
    }

    if info.* == 0 {
        if n <= 1 {
            lwmin = 1; lrwmin = 1; liwmin = 1;
        } else if wantz {
            lwmin = 2*n;
            lrwmin = 1 + 5*n + 2*n*n;
            liwmin = 3 + 5*n;
        } else {
            lwmin = n;
            lrwmin = n;
            liwmin = 1;
        }

        fset(work, 1, _hpgvd_complex_from_real(_hpgvd_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
        fset(rwork, 1, ncast(RT, lrwmin));
        fset(iwork, 1, liwmin);
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

    _hpgvd_pptrf(uplo, n, bp, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    _hpgst(itype, uplo, n, ap, bp.*, info, RT, MT);
    _hpgvd_hpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info, RT, MT);
    lwmin = max(lwmin, ncast(int, real(fget(work.*, 1))));
    lrwmin = max(lrwmin, ncast(int, fget(rwork.*, 1)));
    liwmin = max(liwmin, fget(iwork.*, 1));

    if wantz {
        neig = n;
        if info.* > 0 { neig = info.* - 1; }

        if itype == 1 || itype == 2 {
            if upper { trans = #char "N"; } else { trans = #char "C"; }
            for j:1..neig {
                zcol := _hpgvd_matrix_column(z.*, ldz, j, n);
                _hpgvd_tpsv(uplo, trans, #char "N", n, bp.*, *zcol, 1, MT);
            }
        } else if itype == 3 {
            if upper { trans = #char "C"; } else { trans = #char "N"; }
            for j:1..neig {
                zcol := _hpgvd_matrix_column(z.*, ldz, j, n);
                _hpgvd_tpmv(uplo, trans, #char "N", n, bp.*, *zcol, 1, MT);
            }
        }
    }

    fset(work, 1, _hpgvd_complex_from_real(_hpgvd_roundup_lwork(ncast(RT, lwmin), RT), RT, MT));
    fset(rwork, 1, ncast(RT, lrwmin));
    fset(iwork, 1, liwmin);
}
