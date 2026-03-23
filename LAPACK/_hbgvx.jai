// ? LAPACK *hbgvx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chbgvx/zhbgvx
// *HBGVX computes selected eigenvalues and, optionally, eigenvectors of
// a complex generalized Hermitian-definite banded eigenproblem.

chbgvx :: #bake_arguments _hbgvx(RT=float32, MT=C32);
zhbgvx :: #bake_arguments _hbgvx(RT=float64, MT=C64);

#scope_file
_hbgvx_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hbgvx_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld;
    view.count = count;
    return view;
}

_hbgvx_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hbgvx_pbstf :: (
    $$uplo: u8, $$n: int, $$kb: int, bb: *$BB/ArrayType2D, $$ldbb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && BB.data_type == MT;
} {
    #if MT == C32 { cpbstf(uplo, n, kb, bb, ldbb, info); } else { zpbstf(uplo, n, kb, bb, ldbb, info); }
}

_hbgvx_hbgst :: (
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

_hbgvx_hbtrd :: (
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

_hbgvx_steqr :: (
    $$compz: u8, $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$RWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && Z.data_type == MT && RWORK.data_type == RT;
} {
    #if MT == C32 { csteqr(compz, n, d, e, z, ldz, work, info); } else { zsteqr(compz, n, d, e, z, ldz, work, info); }
}

_hbgvx_sterf :: (
    $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && D.data_type == RT && E.data_type == RT;
} {
    #if RT == float32 { ssterf(n, d, e, info); } else { dsterf(n, d, e, info); }
}

_hbgvx_copy_real :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT && Y.data_type == RT;
} {
    #if RT == float32 { scopy(n, x, incx, y, incy); } else { dcopy(n, x, incx, y, incy); }
}

_hbgvx_copy_complex :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { ccopy(n, x, incx, y, incy); } else { zcopy(n, x, incx, y, incy); }
}

_hbgvx_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 { clacpy(uplo, m, n, a, lda, b, ldb); } else { zlacpy(uplo, m, n, a, lda, b, ldb); }
}

_hbgvx_gemv :: (
    $$trans: u8, $$m: int, $$n: int, alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    x: *$X/ArrayType1D, $$incx: int, beta: MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }
    else { zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }
}

_hbgvx_stebz :: (
    $$range: u8, $$order: u8, $$n: int, vl: $RT, vu: RT, $$il: int, $$iu: int, abstol: RT,
    d: *$D/ArrayType1D, e: *$E/ArrayType1D, m: *int, nsplit: *int, w: *$W/ArrayType1D,
    iblock: *$IBLOCK/ArrayType1D, isplit: *$ISPLIT/ArrayType1D, work: *$WORK/ArrayType1D, iwork: *$IWORK/ArrayType1D,
    info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && D.data_type == RT && E.data_type == RT &&
           W.data_type == RT && WORK.data_type == RT && IBLOCK.data_type == int && ISPLIT.data_type == int && IWORK.data_type == int;
} {
    #if RT == float32 { sstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info); }
    else { dstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info); }
}

_hbgvx_stein :: (
    $$n: int, d: *$D/ArrayType1D, e: *$E/ArrayType1D, $$m: int, w: *$W/ArrayType1D,
    iblock: *$IBLOCK/ArrayType1D, isplit: *$ISPLIT/ArrayType1D, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           D.data_type == RT && E.data_type == RT && W.data_type == RT && WORK.data_type == RT &&
           IBLOCK.data_type == int && ISPLIT.data_type == int && Z.data_type == MT && IWORK.data_type == int && IFAIL.data_type == int;
} {
    #if MT == C32 { cstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info); }
    else { zstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info); }
}

_hbgvx_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hbgvx :: (
    $$jobz: u8, $$range: u8, $$uplo: u8, $$n: int, $$ka: int, $$kb: int,
    ab: *$AB/ArrayType2D, $$ldab: int, bb: *$BB/ArrayType2D, $$ldbb: int, q: *$Q/ArrayType2D, $$ldq: int,
    vl: $RT, vu: RT, $$il: int, $$iu: int, abstol: RT, m: *int, w: *$W/ArrayType1D,
    z: *$Z/ArrayType2D, $$ldz: int, work: *$WORK/ArrayType1D, rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D, ifail: *$IFAIL/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) &&
           AB.data_type == MT && BB.data_type == MT && Q.data_type == MT && W.data_type == RT &&
           Z.data_type == MT && WORK.data_type == MT && RWORK.data_type == RT && IWORK.data_type == int && IFAIL.data_type == int;
} {
    zero_r  :: #run zero(RT);
    czero_t :: #run _hbgvx_complex_from_real(zero_r, RT, MT);
    cone_t  :: #run _hbgvx_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHBGVX" } else { "ZHBGVX" };

    alleig, indeig, test, upper, wantz, valeig : bool = ---;
    i, iinfo, indd, inde, indee, indisp, indiwk, indrwk, indwrk, itmp1, j, jj, nsplit : int = ---;
    order, vect : u8 = ---;
    tmp1 : RT = ---;

    wantz = lsame(jobz, #char "V");
    upper = lsame(uplo, #char "U");
    alleig = lsame(range, #char "A");
    valeig = lsame(range, #char "V");
    indeig = lsame(range, #char "I");

    info.* = 0;
    if !(wantz || lsame(jobz, #char "N")) {
        info.* = -1;
    } else if !(alleig || valeig || indeig) {
        info.* = -2;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -3;
    } else if n < 0 {
        info.* = -4;
    } else if ka < 0 {
        info.* = -5;
    } else if kb < 0 || kb > ka {
        info.* = -6;
    } else if ldab < ka + 1 {
        info.* = -8;
    } else if ldbb < kb + 1 {
        info.* = -10;
    } else if ldq < 1 || (wantz && ldq < n) {
        info.* = -12;
    } else {
        if valeig {
            if n > 0 && vu <= vl {
                info.* = -14;
            }
        } else if indeig {
            if il < 1 || il > max(1, n) {
                info.* = -15;
            } else if iu < min(n, il) || iu > n {
                info.* = -16;
            }
        }
    }
    if info.* == 0 {
        if ldz < 1 || (wantz && ldz < n) {
            info.* = -21;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    m.* = 0;
    if n == 0 {
        return;
    }

    _hbgvx_pbstf(uplo, n, kb, bb, ldbb, info, MT);
    if info.* != 0 {
        info.* = n + info.*;
        return;
    }

    _hbgvx_hbgst(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, work, rwork, *iinfo, RT, MT);

    indd = 1;
    inde = indd + n;
    indrwk = inde + n;
    indwrk = 1;
    if wantz { vect = #char "U"; } else { vect = #char "N"; }

    d_view := _hbgvx_subvector(rwork.*, indd, n);
    e_view := _hbgvx_subvector(rwork.*, inde, n);
    work_view := _hbgvx_subvector(work.*, indwrk, n);
    _hbgvx_hbtrd(vect, uplo, n, ka, ab, ldab, *d_view, *e_view, q, ldq, *work_view, *iinfo, RT, MT);

    test = false;
    if indeig && il == 1 && iu == n {
        test = true;
    }

    if (alleig || test) && (abstol <= zero_r) {
        _hbgvx_copy_real(n, *d_view, 1, w, 1, RT);
        indee = indrwk + 2*n;
        ee_view := _hbgvx_subvector(rwork.*, indee, max(0, n - 1));
        _hbgvx_copy_real(n - 1, *e_view, 1, *ee_view, 1, RT);
        if !wantz {
            _hbgvx_sterf(n, w, *ee_view, info, RT);
        } else {
            _hbgvx_lacpy(#char "A", n, n, q.*, ldq, z, ldz, MT);
            rwork_tail := _hbgvx_subvector(rwork.*, indrwk, 7*n - indrwk + 1);
            _hbgvx_steqr(jobz, n, w, *ee_view, z, ldz, *rwork_tail, info, RT, MT);
            if info.* == 0 {
                for i: 1..n {
                    fset(ifail, i, 0);
                }
            }
        }
        if info.* == 0 {
            m.* = n;
        } else {
            info.* = 0;
        }
    }

    if m.* == 0 {
        if wantz { order = #char "B"; } else { order = #char "E"; }
        indisp = 1 + n;
        indiwk = indisp + n;
        iblock_view := _hbgvx_subvector(iwork.*, 1, n);
        isplit_view := _hbgvx_subvector(iwork.*, indisp, n);
        iwork_tail := _hbgvx_subvector(iwork.*, indiwk, 5*n - indiwk + 1);
        rwork_tail := _hbgvx_subvector(rwork.*, indrwk, 7*n - indrwk + 1);

        _hbgvx_stebz(range, order, n, vl, vu, il, iu, abstol, *d_view, *e_view, m, *nsplit, w, *iblock_view, *isplit_view, *rwork_tail, *iwork_tail, info, RT);

        if wantz {
            _hbgvx_stein(n, *d_view, *e_view, m.*, w, *iblock_view, *isplit_view, z, ldz, *rwork_tail, *iwork_tail, ifail, info, RT, MT);
            for j: 1..m.* {
                work_col := _hbgvx_subvector(work.*, 1, n);
                z_col := _hbgvx_matrix_column(z.*, ldz, j, n);
                _hbgvx_copy_complex(n, *z_col, 1, *work_col, 1, MT);
                _hbgvx_gemv(#char "N", n, n, cone_t, q.*, ldq, *work_col, 1, czero_t, *z_col, 1, MT);
            }
        }
    }

    if wantz {
        for j: 1..m.* - 1 {
            i = 0;
            tmp1 = fget(w, j);
            for jj: j + 1..m.* {
                if fget(w, jj) < tmp1 {
                    i = jj;
                    tmp1 = fget(w, jj);
                }
            }
            if i != 0 {
                itmp1 = fget(iwork, i);
                fset(w, i, fget(w, j));
                fset(iwork, i, fget(iwork, j));
                fset(w, j, tmp1);
                fset(iwork, j, itmp1);

                zi_col := _hbgvx_matrix_column(z.*, ldz, i, n);
                zj_col := _hbgvx_matrix_column(z.*, ldz, j, n);
                _hbgvx_swap(n, *zi_col, 1, *zj_col, 1, MT);

                if info.* != 0 {
                    itmp1 = fget(ifail, i);
                    fset(ifail, i, fget(ifail, j));
                    fset(ifail, j, itmp1);
                }
            }
        }
    }
}
