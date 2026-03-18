// ? LAPACK *geqp3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeqp3/cgeqp3/dgeqp3/zgeqp3
// *GEQP3 computes a QR factorization with column pivoting.

sgeqp3 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int
) {
    _geqp3_real(m, n, a, lda, jpvt, tau, work, lwork, info, RT=float32);
}

dgeqp3 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int
) {
    _geqp3_real(m, n, a, lda, jpvt, tau, work, lwork, info, RT=float64);
}

cgeqp3 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, info: *int
) {
    _geqp3_complex(m, n, a, lda, jpvt, tau, work, lwork, rwork, info, RT=float32, MT=C32);
}

zgeqp3 :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, info: *int
) {
    _geqp3_complex(m, n, a, lda, jpvt, tau, work, lwork, rwork, info, RT=float64, MT=C64);
}

#scope_file
_geqp3_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_geqp3_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_geqp3_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _geqp3_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_geqp3_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _geqp3_raw_slice(a.data + (col - 1)*lda + (row - 1), max(0, count));
}

_geqp3_roundup_lwork :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_geqp3_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && (MT == C32 || MT == C64);
} {
    return cpx(x, zero(RT));
}

_geqp3_set_work1_real :: inline (work: *$WORK/ArrayType1D, value: $RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && WORK.data_type == RT;
} {
    fset(work, 1, value);
}

_geqp3_set_work1_complex :: inline (work: *$WORK/ArrayType1D, value: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && (MT == C32 || MT == C64) && WORK.data_type == MT;
} {
    fset(work, 1, _geqp3_complex_from_real(value, RT, MT));
}

_geqp3_work1_to_int :: inline (work1: $T, $RT: Type, $T: Type) -> int #modify {
    return (RT == float32 || RT == float64) && (T == RT || T == complex_with_inner_type(RT));
} {
    #if T == float32 || T == float64 {
        return ncast(int, work1);
    } else {
        return ncast(int, real(work1));
    }
}

_geqp3_column_swap :: (
    $$m: int, a: *$A/ArrayType2D, $$lda: int, $$j1: int, $$j2: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) && A.data_type == MT;
} {
    if j1 == j2 || m <= 0 return;

    col1 := _geqp3_matrix_column(a.*, lda, 1, j1, m);
    col2 := _geqp3_matrix_column(a.*, lda, 1, j2, m);
    #if MT == float32 {
        sswap(m, col1, 1, col2, 1);
    } else if MT == float64 {
        dswap(m, col1, 1, col2, 1);
    } else if MT == C32 {
        cswap(m, col1, 1, col2, 1);
    } else {
        zswap(m, col1, 1, col2, 1);
    }
}

_geqp3_nrm2 :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return snrm2(n, x, incx);
    } else {
        return dnrm2(n, x, incx);
    }
}

_geqp3_cnrm2 :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        return scnrm2(n, x, incx);
    } else {
        return dznrm2(n, x, incx);
    }
}

_geqp3_geqrf :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tau: *$TAU/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == float64 {
        dgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else if MT == C32 {
        cgeqrf(m, n, a, lda, tau, work, lwork, info);
    } else {
        zgeqrf(m, n, a, lda, tau, work, lwork, info);
    }
}

_geqp3_mqr :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$k: int,
    a: $A/ArrayType2D, $$lda: int, tau: $TAU/ArrayType1D,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        sormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == float64 {
        dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else if MT == C32 {
        cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_geqp3_laqps_real :: (
    $$m: int, $$n: int, $$offset: int, $$nb: int, kb: *int,
    a: *$A/ArrayType2D, $$lda: int, jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, auxv: *$AUXV/ArrayType1D,
    f: *$F/ArrayType2D, $$ldf: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPVT.data_type == int &&
           TAU.data_type == RT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           AUXV.data_type == RT &&
           F.data_type == RT;
} {
    #if RT == float32 {
        slaqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
    } else {
        dlaqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
    }
}

_geqp3_laqp2_real :: (
    $$m: int, $$n: int, $$offset: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D, vn1: *$VN1/ArrayType1D,
    vn2: *$VN2/ArrayType1D, work: *$WORK/ArrayType1D, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPVT.data_type == int &&
           TAU.data_type == RT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           WORK.data_type == RT;
} {
    #if RT == float32 {
        slaqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
    } else {
        dlaqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
    }
}

_geqp3_laqps_complex :: (
    $$m: int, $$n: int, $$offset: int, $$nb: int, kb: *int,
    a: *$A/ArrayType2D, $$lda: int, jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, auxv: *$AUXV/ArrayType1D,
    f: *$F/ArrayType2D, $$ldf: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           JPVT.data_type == int &&
           TAU.data_type == MT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           AUXV.data_type == MT &&
           F.data_type == MT;
} {
    #if MT == C32 {
        claqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
    } else {
        zlaqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
    }
}

_geqp3_laqp2_complex :: (
    $$m: int, $$n: int, $$offset: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D, vn1: *$VN1/ArrayType1D,
    vn2: *$VN2/ArrayType1D, work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           JPVT.data_type == int &&
           TAU.data_type == MT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           WORK.data_type == MT;
} {
    #if MT == C32 {
        claqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
    } else {
        zlaqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
    }
}

_geqp3_real :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPVT.data_type == int &&
           TAU.data_type == RT &&
           WORK.data_type == RT;
} {
    routine_name :: #run if RT == float32 then "SGEQP3" else "DGEQP3";
    qrf_name :: #run if RT == float32 then "SGEQRF" else "DGEQRF";

    lquery : bool = ---;
    fjb, iws, j, jb, lwkopt, minmn, minws, na, nb, nbmin, nfxd, nx, sm, sminmn, sn, topbmn : int = ---;

    info.* = 0;
    lquery = lwork == -1;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    }

    if info.* == 0 {
        minmn = min(m, n);
        if minmn == 0 {
            iws = 1;
            lwkopt = 1;
        } else {
            iws = 3*n + 1;
            nb = ilaenv(1, qrf_name, " "s, m, n, -1, -1);
            lwkopt = 2*n + (n + 1)*nb;
        }
        _geqp3_set_work1_real(work, _geqp3_roundup_lwork(ncast(RT, lwkopt), RT), RT);

        if lwork < iws && !lquery {
            info.* = -8;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    nfxd = 1;
    for j:1..n {
        if fget(jpvt.*, j) != 0 {
            if j != nfxd {
                _geqp3_column_swap(m, a, lda, j, nfxd, RT);
                jp := fget(jpvt.*, nfxd);
                fset(jpvt, j, jp);
                fset(jpvt, nfxd, j);
            } else {
                fset(jpvt, j, j);
            }
            nfxd += 1;
        } else {
            fset(jpvt, j, j);
        }
    }
    nfxd -= 1;

    if nfxd > 0 {
        na = min(m, nfxd);
        _geqp3_geqrf(m, na, a, lda, tau, work, lwork, info, RT, RT);
        iws = max(iws, _geqp3_work1_to_int(fget(work.*, 1), RT, RT));
        if na < n {
            c_raw := _geqp3_matrix_storage(a.*, lda, 1, na + 1, m, n - na);
            c_view := array2d_view(c_raw, lda, row_major=false);
            _geqp3_mqr(#char "L", #char "T", m, n - na, na, a.*, lda, tau.*, *c_view, lda, work, lwork, info, RT, RT);
            iws = max(iws, _geqp3_work1_to_int(fget(work.*, 1), RT, RT));
        }
    }

    if nfxd < minmn {
        sm = m - nfxd;
        sn = n - nfxd;
        sminmn = minmn - nfxd;

        nb = ilaenv(1, qrf_name, " "s, sm, sn, -1, -1);
        nbmin = 2;
        nx = 0;
        if nb > 1 && nb < sminmn {
            nx = max(0, ilaenv(3, qrf_name, " "s, sm, sn, -1, -1));
            if nx < sminmn {
                minws = 2*sn + (sn + 1)*nb;
                iws = max(iws, minws);
                if lwork < minws {
                    nb = (lwork - 2*sn) / (sn + 1);
                    nbmin = max(2, ilaenv(2, qrf_name, " "s, sm, sn, -1, -1));
                }
            }
        }

        for j:nfxd + 1..n {
            col := _geqp3_matrix_column(a.*, lda, nfxd + 1, j, sm);
            norm := _geqp3_nrm2(sm, col, 1, RT);
            fset(work, j, norm);
            fset(work, n + j, norm);
        }

        if nb >= nbmin && nb < sminmn && nx < sminmn {
            j = nfxd + 1;
            topbmn = minmn - nx;
            while j <= topbmn {
                jb = min(nb, topbmn - j + 1);

                a_sub_raw := _geqp3_matrix_storage(a.*, lda, 1, j, m, n - j + 1);
                a_sub := array2d_view(a_sub_raw, lda, row_major=false);
                jpvt_sub := _geqp3_subvector(jpvt.*, j, n - j + 1);
                tau_sub := _geqp3_subvector(tau.*, j, n - j + 1);
                vn1 := _geqp3_subvector(work.*, j, n - j + 1);
                vn2 := _geqp3_subvector(work.*, n + j, n - j + 1);
                auxv := _geqp3_subvector(work.*, 2*n + 1, jb);
                f_raw := _geqp3_raw_slice(work.data + (2*n + jb), max(0, (n - j + 1)*jb));
                f_view := array2d_view(f_raw, n - j + 1, row_major=false);

                _geqp3_laqps_real(m, n - j + 1, j - 1, jb, *fjb, *a_sub, lda, *jpvt_sub, *tau_sub, *vn1, *vn2, *auxv, *f_view, n - j + 1, RT);
                j += fjb;
            }
        } else {
            j = nfxd + 1;
        }

        if j <= minmn {
            a_sub_raw := _geqp3_matrix_storage(a.*, lda, 1, j, m, n - j + 1);
            a_sub := array2d_view(a_sub_raw, lda, row_major=false);
            jpvt_sub := _geqp3_subvector(jpvt.*, j, n - j + 1);
            tau_sub := _geqp3_subvector(tau.*, j, n - j + 1);
            vn1 := _geqp3_subvector(work.*, j, n - j + 1);
            vn2 := _geqp3_subvector(work.*, n + j, n - j + 1);
            auxv := _geqp3_subvector(work.*, 2*n + 1, max(0, n - j + 1));
            _geqp3_laqp2_real(m, n - j + 1, j - 1, *a_sub, lda, *jpvt_sub, *tau_sub, *vn1, *vn2, *auxv, RT);
        }
    }

    _geqp3_set_work1_real(work, _geqp3_roundup_lwork(ncast(RT, iws), RT), RT);
}

_geqp3_complex :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    jpvt: *$JPVT/ArrayType1D, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, rwork: *$RWORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           JPVT.data_type == int &&
           TAU.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    routine_name :: #run if MT == C32 then "CGEQP3" else "ZGEQP3";
    qrf_name :: #run if MT == C32 then "CGEQRF" else "ZGEQRF";

    lquery : bool = ---;
    fjb, iws, j, jb, lwkopt, minmn, minws, na, nb, nbmin, nfxd, nx, sm, sminmn, sn, topbmn : int = ---;

    info.* = 0;
    lquery = lwork == -1;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    }

    if info.* == 0 {
        minmn = min(m, n);
        if minmn == 0 {
            iws = 1;
            lwkopt = 1;
        } else {
            iws = n + 1;
            nb = ilaenv(1, qrf_name, " "s, m, n, -1, -1);
            lwkopt = (n + 1)*nb;
        }
        _geqp3_set_work1_complex(work, ncast(RT, lwkopt), RT, MT);

        if lwork < iws && !lquery {
            info.* = -8;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    nfxd = 1;
    for j:1..n {
        if fget(jpvt.*, j) != 0 {
            if j != nfxd {
                _geqp3_column_swap(m, a, lda, j, nfxd, MT);
                jp := fget(jpvt.*, nfxd);
                fset(jpvt, j, jp);
                fset(jpvt, nfxd, j);
            } else {
                fset(jpvt, j, j);
            }
            nfxd += 1;
        } else {
            fset(jpvt, j, j);
        }
    }
    nfxd -= 1;

    if nfxd > 0 {
        na = min(m, nfxd);
        _geqp3_geqrf(m, na, a, lda, tau, work, lwork, info, RT, MT);
        iws = max(iws, _geqp3_work1_to_int(fget(work.*, 1), RT, MT));
        if na < n {
            c_raw := _geqp3_matrix_storage(a.*, lda, 1, na + 1, m, n - na);
            c_view := array2d_view(c_raw, lda, row_major=false);
            _geqp3_mqr(#char "L", #char "C", m, n - na, na, a.*, lda, tau.*, *c_view, lda, work, lwork, info, RT, MT);
            iws = max(iws, _geqp3_work1_to_int(fget(work.*, 1), RT, MT));
        }
    }

    if nfxd < minmn {
        sm = m - nfxd;
        sn = n - nfxd;
        sminmn = minmn - nfxd;

        nb = ilaenv(1, qrf_name, " "s, sm, sn, -1, -1);
        nbmin = 2;
        nx = 0;
        if nb > 1 && nb < sminmn {
            nx = max(0, ilaenv(3, qrf_name, " "s, sm, sn, -1, -1));
            if nx < sminmn {
                minws = (sn + 1)*nb;
                iws = max(iws, minws);
                if lwork < minws {
                    nb = lwork / (sn + 1);
                    nbmin = max(2, ilaenv(2, qrf_name, " "s, sm, sn, -1, -1));
                }
            }
        }

        for j:nfxd + 1..n {
            col := _geqp3_matrix_column(a.*, lda, nfxd + 1, j, sm);
            norm := _geqp3_cnrm2(sm, col, 1, RT, MT);
            fset(rwork, j, norm);
            fset(rwork, n + j, norm);
        }

        if nb >= nbmin && nb < sminmn && nx < sminmn {
            j = nfxd + 1;
            topbmn = minmn - nx;
            while j <= topbmn {
                jb = min(nb, topbmn - j + 1);

                a_sub_raw := _geqp3_matrix_storage(a.*, lda, 1, j, m, n - j + 1);
                a_sub := array2d_view(a_sub_raw, lda, row_major=false);
                jpvt_sub := _geqp3_subvector(jpvt.*, j, n - j + 1);
                tau_sub := _geqp3_subvector(tau.*, j, n - j + 1);
                vn1 := _geqp3_subvector(rwork.*, j, n - j + 1);
                vn2 := _geqp3_subvector(rwork.*, n + j, n - j + 1);
                auxv := _geqp3_subvector(work.*, 1, jb);
                f_raw := _geqp3_raw_slice(work.data + jb, max(0, (n - j + 1)*jb));
                f_view := array2d_view(f_raw, n - j + 1, row_major=false);

                _geqp3_laqps_complex(m, n - j + 1, j - 1, jb, *fjb, *a_sub, lda, *jpvt_sub, *tau_sub, *vn1, *vn2, *auxv, *f_view, n - j + 1, RT, MT);
                j += fjb;
            }
        } else {
            j = nfxd + 1;
        }

        if j <= minmn {
            a_sub_raw := _geqp3_matrix_storage(a.*, lda, 1, j, m, n - j + 1);
            a_sub := array2d_view(a_sub_raw, lda, row_major=false);
            jpvt_sub := _geqp3_subvector(jpvt.*, j, n - j + 1);
            tau_sub := _geqp3_subvector(tau.*, j, n - j + 1);
            vn1 := _geqp3_subvector(rwork.*, j, n - j + 1);
            vn2 := _geqp3_subvector(rwork.*, n + j, n - j + 1);
            auxv := _geqp3_subvector(work.*, 1, max(0, n - j + 1));
            _geqp3_laqp2_complex(m, n - j + 1, j - 1, *a_sub, lda, *jpvt_sub, *tau_sub, *vn1, *vn2, *auxv, RT, MT);
        }
    }

    _geqp3_set_work1_complex(work, ncast(RT, lwkopt), RT, MT);
}
