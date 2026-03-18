// ? LAPACK *gesdd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgesdd/cgesdd/dgesdd/zgesdd
// *GESDD computes the singular value decomposition using
// a divide-and-conquer algorithm.

sgesdd :: #bake_arguments _gesdd_real(RT=float32);
dgesdd :: #bake_arguments _gesdd_real(RT=float64);
cgesdd :: #bake_arguments _gesdd_complex(RT=float32, MT=C32);
zgesdd :: #bake_arguments _gesdd_complex(RT=float64, MT=C64);

#scope_file
_gesdd_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_gesdd_roundup_lwork :: (x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(x);
    } else {
        return droundup_lwork(x);
    }
}

_gesdd_set_work1 :: inline (work: *$WORK/ArrayType1D, x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           ((MT == RT) || MT == complex_with_inner_type(RT)) &&
           WORK.data_type == MT;
} {
    #if MT == RT {
        fset(work, 1, x);
    } else {
        fset(work, 1, _gesdd_complex_from_real(x, RT, MT));
    }
}

_gesdd_query_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == RT {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_gesdd_real_query_gebrd :: ($$m: int, $$n: int, $RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    dum := array1d(1, zero(RT));
    a := array2d_view(dum, 1, row_major=false);
    info := 0;
    #if RT == float32 {
        sgebrd(m, n, *a, 1, *dum, *dum, *dum, *dum, *dum, -1, *info);
    } else {
        dgebrd(m, n, *a, 1, *dum, *dum, *dum, *dum, *dum, -1, *info);
    }
    return _gesdd_query_to_int(fget(dum, 1), RT, RT);
}

_gesdd_real_query_geqrf :: ($$m: int, $$n: int, $RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    dum := array1d(1, zero(RT));
    a := array2d_view(dum, 1, row_major=false);
    info := 0;
    #if RT == float32 {
        sgeqrf(m, n, *a, 1, *dum, *dum, -1, *info);
    } else {
        dgeqrf(m, n, *a, 1, *dum, *dum, -1, *info);
    }
    return _gesdd_query_to_int(fget(dum, 1), RT, RT);
}

_gesdd_real_query_gelqf :: ($$m: int, $$n: int, $RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    dum := array1d(1, zero(RT));
    a := array2d_view(dum, 1, row_major=false);
    info := 0;
    #if RT == float32 {
        sgelqf(m, n, *a, 1, *dum, *dum, -1, *info);
    } else {
        dgelqf(m, n, *a, 1, *dum, *dum, -1, *info);
    }
    return _gesdd_query_to_int(fget(dum, 1), RT, RT);
}

_gesdd_real_query_orgbr :: ($$vect: u8, $$m: int, $$n: int, $$k: int, $RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    dum := array1d(1, zero(RT));
    a := array2d_view(dum, 1, row_major=false);
    info := 0;
    #if RT == float32 {
        sorgbr(vect, m, n, k, *a, 1, *dum, *dum, -1, *info);
    } else {
        dorgbr(vect, m, n, k, *a, 1, *dum, *dum, -1, *info);
    }
    return _gesdd_query_to_int(fget(dum, 1), RT, RT);
}

_gesdd_real_query_orgqr :: ($$m: int, $$n: int, $$k: int, $RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    dum := array1d(1, zero(RT));
    a := array2d_view(dum, 1, row_major=false);
    info := 0;
    #if RT == float32 {
        sorgqr(m, n, k, *a, 1, *dum, *dum, -1, *info);
    } else {
        dorgqr(m, n, k, *a, 1, *dum, *dum, -1, *info);
    }
    return _gesdd_query_to_int(fget(dum, 1), RT, RT);
}

_gesdd_real_query_orglq :: ($$m: int, $$n: int, $$k: int, $RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    dum := array1d(1, zero(RT));
    a := array2d_view(dum, 1, row_major=false);
    info := 0;
    #if RT == float32 {
        sorglq(m, n, k, *a, 1, *dum, *dum, -1, *info);
    } else {
        dorglq(m, n, k, *a, 1, *dum, *dum, -1, *info);
    }
    return _gesdd_query_to_int(fget(dum, 1), RT, RT);
}

_gesdd_real_query_ormbr :: (
    $$vect: u8, $$side: u8, $$trans: u8,
    $$m: int, $$n: int, $$k: int,
    $RT: Type
) -> int #modify {
    return RT == float32 || RT == float64;
} {
    dum := array1d(1, zero(RT));
    a := array2d_view(dum, 1, row_major=false);
    c := array2d_view(dum, 1, row_major=false);
    info := 0;
    #if RT == float32 {
        sormbr(vect, side, trans, m, n, k, *a, 1, *dum, *c, 1, *dum, -1, *info);
    } else {
        dormbr(vect, side, trans, m, n, k, *a, 1, *dum, *c, 1, *dum, -1, *info);
    }
    return _gesdd_query_to_int(fget(dum, 1), RT, RT);
}

_gesdd_complex_query_gebrd :: ($$m: int, $$n: int, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    cdum := array1d(1, zero(MT));
    rdum := array1d(1, zero(RT));
    a := array2d_view(cdum, 1, row_major=false);
    info := 0;
    #if MT == C32 {
        cgebrd(m, n, *a, 1, *rdum, *rdum, *cdum, *cdum, *cdum, -1, *info);
    } else {
        zgebrd(m, n, *a, 1, *rdum, *rdum, *cdum, *cdum, *cdum, -1, *info);
    }
    return _gesdd_query_to_int(fget(cdum, 1), RT, MT);
}

_gesdd_complex_query_geqrf :: ($$m: int, $$n: int, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    cdum := array1d(1, zero(MT));
    a := array2d_view(cdum, 1, row_major=false);
    info := 0;
    #if MT == C32 {
        cgeqrf(m, n, *a, 1, *cdum, *cdum, -1, *info);
    } else {
        zgeqrf(m, n, *a, 1, *cdum, *cdum, -1, *info);
    }
    return _gesdd_query_to_int(fget(cdum, 1), RT, MT);
}

_gesdd_complex_query_gelqf :: ($$m: int, $$n: int, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    cdum := array1d(1, zero(MT));
    a := array2d_view(cdum, 1, row_major=false);
    info := 0;
    #if MT == C32 {
        cgelqf(m, n, *a, 1, *cdum, *cdum, -1, *info);
    } else {
        zgelqf(m, n, *a, 1, *cdum, *cdum, -1, *info);
    }
    return _gesdd_query_to_int(fget(cdum, 1), RT, MT);
}

_gesdd_complex_query_ungbr :: ($$vect: u8, $$m: int, $$n: int, $$k: int, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    cdum := array1d(1, zero(MT));
    a := array2d_view(cdum, 1, row_major=false);
    info := 0;
    #if MT == C32 {
        cungbr(vect, m, n, k, *a, 1, *cdum, *cdum, -1, *info);
    } else {
        zungbr(vect, m, n, k, *a, 1, *cdum, *cdum, -1, *info);
    }
    return _gesdd_query_to_int(fget(cdum, 1), RT, MT);
}

_gesdd_complex_query_ungqr :: ($$m: int, $$n: int, $$k: int, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    cdum := array1d(1, zero(MT));
    a := array2d_view(cdum, 1, row_major=false);
    info := 0;
    #if MT == C32 {
        cungqr(m, n, k, *a, 1, *cdum, *cdum, -1, *info);
    } else {
        zungqr(m, n, k, *a, 1, *cdum, *cdum, -1, *info);
    }
    return _gesdd_query_to_int(fget(cdum, 1), RT, MT);
}

_gesdd_complex_query_unglq :: ($$m: int, $$n: int, $$k: int, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    cdum := array1d(1, zero(MT));
    a := array2d_view(cdum, 1, row_major=false);
    info := 0;
    #if MT == C32 {
        cunglq(m, n, k, *a, 1, *cdum, *cdum, -1, *info);
    } else {
        zunglq(m, n, k, *a, 1, *cdum, *cdum, -1, *info);
    }
    return _gesdd_query_to_int(fget(cdum, 1), RT, MT);
}

_gesdd_complex_query_unmbr :: (
    $$vect: u8, $$side: u8, $$trans: u8,
    $$m: int, $$n: int, $$k: int,
    $RT: Type, $MT: Type
) -> int #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    cdum := array1d(1, zero(MT));
    a := array2d_view(cdum, 1, row_major=false);
    c := array2d_view(cdum, 1, row_major=false);
    info := 0;
    #if MT == C32 {
        cunmbr(vect, side, trans, m, n, k, *a, 1, *cdum, *c, 1, *cdum, -1, *info);
    } else {
        zunmbr(vect, side, trans, m, n, k, *a, 1, *cdum, *c, 1, *cdum, -1, *info);
    }
    return _gesdd_query_to_int(fget(cdum, 1), RT, MT);
}

_gesdd_real :: (
    $$jobz: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    s: *$S/ArrayType1D,
    u: *$U/ArrayType2D,
    $$ldu: int,
    vt: *$VT/ArrayType2D,
    $$ldvt: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           S.data_type == RT &&
           U.data_type == RT &&
           VT.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if RT == float32 { "SGESDD" } else { "DGESDD" };

    // First-stage staging note:
    // This file preserves the full LAPACK workspace-query logic and the
    // path split for *GESDD. The remaining source-shaped execution body
    // still needs to be filled in path by path:
    // - M >= N: Paths 1..5
    // - N >  M: Paths 1t..5t
    // including the QR/LQ reductions, bidiagonalization, DBDSDC,
    // ORGBR/ORGQR/ORGLQ, ORMBR backtransforms, and the fast/slow
    // overwrite cases.

    lquery, wntqa, wntqas, wntqn, wntqo, wntqs : bool = ---;
    bdspac, maxwrk, minmn, minwrk, mnthr, wrkbl : int = ---;
    lwork_gebrd_mn, lwork_gebrd_mm, lwork_gebrd_nn, lwork_gelqf_mn, lwork_geqrf_mn : int = ---;
    lwork_orgbr_q_nn, lwork_orglq_mn, lwork_orglq_nn, lwork_orgqr_mm, lwork_orgqr_mn : int = ---;
    lwork_ormbr_prt_mm, lwork_ormbr_prt_mn, lwork_ormbr_prt_nn, lwork_ormbr_qln_mm, lwork_ormbr_qln_mn, lwork_ormbr_qln_nn : int = ---;

    info.* = 0;
    minmn = min(m, n);
    wntqa = char_eq(jobz, #char "a");
    wntqs = char_eq(jobz, #char "s");
    wntqas = wntqa || wntqs;
    wntqo = char_eq(jobz, #char "o");
    wntqn = char_eq(jobz, #char "n");
    lquery = lwork == -1;

    if !(wntqa || wntqs || wntqo || wntqn) {
        info.* = -1;
    } else if m < 0 {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, m) {
        info.* = -5;
    } else if ldu < 1 || (wntqas && ldu < m) || (wntqo && m < n && ldu < m) {
        info.* = -8;
    } else if ldvt < 1 || (wntqa && ldvt < n) || (wntqs && ldvt < minmn) || (wntqo && m >= n && ldvt < n) {
        info.* = -10;
    }

    if info.* == 0 {
        minwrk = 1;
        maxwrk = 1;
        bdspac = 0;
        mnthr = cast(int) (ncast(RT, minmn) * ncast(RT, 11) / ncast(RT, 6));

        if m >= n && minmn > 0 {
            bdspac = ifx wntqn then 7*n else 3*n*n + 4*n;

            lwork_gebrd_mn = _gesdd_real_query_gebrd(m, n, RT);
            lwork_gebrd_nn = _gesdd_real_query_gebrd(n, n, RT);
            lwork_geqrf_mn = _gesdd_real_query_geqrf(m, n, RT);
            lwork_orgbr_q_nn = _gesdd_real_query_orgbr(#char "q", n, n, n, RT);
            lwork_orgqr_mm = _gesdd_real_query_orgqr(m, m, n, RT);
            lwork_orgqr_mn = _gesdd_real_query_orgqr(m, n, n, RT);
            lwork_ormbr_prt_nn = _gesdd_real_query_ormbr(#char "p", #char "r", #char "t", n, n, n, RT);
            lwork_ormbr_qln_nn = _gesdd_real_query_ormbr(#char "q", #char "l", #char "n", n, n, n, RT);
            lwork_ormbr_qln_mn = _gesdd_real_query_ormbr(#char "q", #char "l", #char "n", m, n, n, RT);
            lwork_ormbr_qln_mm = _gesdd_real_query_ormbr(#char "q", #char "l", #char "n", m, m, n, RT);

            if m >= mnthr {
                if wntqn {
                    wrkbl = n + lwork_geqrf_mn;
                    wrkbl = max(wrkbl, 3*n + lwork_gebrd_nn);
                    maxwrk = max(wrkbl, bdspac + n);
                    minwrk = bdspac + n;
                } else if wntqo {
                    wrkbl = n + lwork_geqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_orgqr_mn);
                    wrkbl = max(wrkbl, 3*n + lwork_gebrd_nn);
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_qln_nn);
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_prt_nn);
                    wrkbl = max(wrkbl, 3*n + bdspac);
                    maxwrk = wrkbl + 2*n*n;
                    minwrk = bdspac + 2*n*n + 3*n;
                } else if wntqs {
                    wrkbl = n + lwork_geqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_orgqr_mn);
                    wrkbl = max(wrkbl, 3*n + lwork_gebrd_nn);
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_qln_nn);
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_prt_nn);
                    wrkbl = max(wrkbl, 3*n + bdspac);
                    maxwrk = wrkbl + n*n;
                    minwrk = bdspac + n*n + 3*n;
                } else {
                    wrkbl = n + lwork_geqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_orgqr_mm);
                    wrkbl = max(wrkbl, 3*n + lwork_gebrd_nn);
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_qln_nn);
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_prt_nn);
                    wrkbl = max(wrkbl, 3*n + bdspac);
                    maxwrk = wrkbl + n*n;
                    minwrk = n*n + max(3*n + bdspac, n + m);
                }
            } else {
                wrkbl = 3*n + lwork_gebrd_mn;
                if wntqn {
                    maxwrk = max(wrkbl, 3*n + bdspac);
                    minwrk = 3*n + max(m, bdspac);
                } else if wntqo {
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_prt_nn);
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_qln_mn);
                    wrkbl = max(wrkbl, 3*n + bdspac);
                    maxwrk = wrkbl + m*n;
                    minwrk = 3*n + max(m, n*n + bdspac);
                } else if wntqs {
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_qln_mn);
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_prt_nn);
                    maxwrk = max(wrkbl, 3*n + bdspac);
                    minwrk = 3*n + max(m, bdspac);
                } else {
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_qln_mm);
                    wrkbl = max(wrkbl, 3*n + lwork_ormbr_prt_nn);
                    maxwrk = max(wrkbl, 3*n + bdspac);
                    minwrk = 3*n + max(m, bdspac);
                }
            }
        } else if minmn > 0 {
            bdspac = ifx wntqn then 7*m else 3*m*m + 4*m;

            lwork_gebrd_mn = _gesdd_real_query_gebrd(m, n, RT);
            lwork_gebrd_mm = _gesdd_real_query_gebrd(m, m, RT);
            lwork_gelqf_mn = _gesdd_real_query_gelqf(m, n, RT);
            lwork_orglq_nn = _gesdd_real_query_orglq(n, n, m, RT);
            lwork_orglq_mn = _gesdd_real_query_orglq(m, n, m, RT);
            lwork_ormbr_prt_mm = _gesdd_real_query_ormbr(#char "p", #char "r", #char "t", m, m, m, RT);
            lwork_ormbr_prt_mn = _gesdd_real_query_ormbr(#char "p", #char "r", #char "t", m, n, m, RT);
            lwork_ormbr_prt_nn = _gesdd_real_query_ormbr(#char "p", #char "r", #char "t", n, n, m, RT);
            lwork_ormbr_qln_mm = _gesdd_real_query_ormbr(#char "q", #char "l", #char "n", m, m, m, RT);

            if n >= mnthr {
                if wntqn {
                    wrkbl = m + lwork_gelqf_mn;
                    wrkbl = max(wrkbl, 3*m + lwork_gebrd_mm);
                    maxwrk = max(wrkbl, bdspac + m);
                    minwrk = bdspac + m;
                } else if wntqo {
                    wrkbl = m + lwork_gelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_orglq_mn);
                    wrkbl = max(wrkbl, 3*m + lwork_gebrd_mm);
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_qln_mm);
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_prt_mm);
                    wrkbl = max(wrkbl, 3*m + bdspac);
                    maxwrk = wrkbl + 2*m*m;
                    minwrk = bdspac + 2*m*m + 3*m;
                } else if wntqs {
                    wrkbl = m + lwork_gelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_orglq_mn);
                    wrkbl = max(wrkbl, 3*m + lwork_gebrd_mm);
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_qln_mm);
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_prt_mm);
                    wrkbl = max(wrkbl, 3*m + bdspac);
                    maxwrk = wrkbl + m*m;
                    minwrk = bdspac + m*m + 3*m;
                } else {
                    wrkbl = m + lwork_gelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_orglq_nn);
                    wrkbl = max(wrkbl, 3*m + lwork_gebrd_mm);
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_qln_mm);
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_prt_mm);
                    wrkbl = max(wrkbl, 3*m + bdspac);
                    maxwrk = wrkbl + m*m;
                    minwrk = m*m + max(3*m + bdspac, m + n);
                }
            } else {
                wrkbl = 3*m + lwork_gebrd_mn;
                if wntqn {
                    maxwrk = max(wrkbl, 3*m + bdspac);
                    minwrk = 3*m + max(n, bdspac);
                } else if wntqo {
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_qln_mm);
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_prt_mn);
                    wrkbl = max(wrkbl, 3*m + bdspac);
                    maxwrk = wrkbl + m*n;
                    minwrk = 3*m + max(n, m*m + bdspac);
                } else if wntqs {
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_qln_mm);
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_prt_mn);
                    maxwrk = max(wrkbl, 3*m + bdspac);
                    minwrk = 3*m + max(n, bdspac);
                } else {
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_qln_mm);
                    wrkbl = max(wrkbl, 3*m + lwork_ormbr_prt_nn);
                    maxwrk = max(wrkbl, 3*m + bdspac);
                    minwrk = 3*m + max(n, bdspac);
                }
            }
        }

        maxwrk = max(maxwrk, minwrk);
        _gesdd_set_work1(work, _gesdd_roundup_lwork(ncast(RT, maxwrk), RT), RT, RT);
        if lwork < minwrk && !lquery {
            info.* = -12;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if m == 0 || n == 0 {
        return;
    }
}

_gesdd_complex :: (
    $$jobz: u8,
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    s: *$S/ArrayType1D,
    u: *$U/ArrayType2D,
    $$ldu: int,
    vt: *$VT/ArrayType2D,
    $$ldvt: int,
    work: *$WORK/ArrayType1D,
    $$lwork: int,
    rwork: *$RWORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           S.data_type == RT &&
           U.data_type == MT &&
           VT.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT &&
           IWORK.data_type == int;
} {
    routine_name :: #run if MT == C32 { "CGESDD" } else { "ZGESDD" };

    // First-stage staging note:
    // This keeps the full complex workspace-query path split from
    // LAPACK, including the MNTHR1/MNTHR2 thresholds and the separate
    // complex-vs-real workspace formulas. The execution body still
    // needs the source-shaped QR/LQ, bidiagonal, DBDSDC, UNGBR/UNGQR/
    // UNGLQ, UNMBR, and overwrite fast/slow path translations.

    lquery, wntqa, wntqas, wntqn, wntqo, wntqs : bool = ---;
    maxwrk, minmn, minwrk, mnthr1, mnthr2, wrkbl : int = ---;
    lwork_zgebrd_mn, lwork_zgebrd_mm, lwork_zgebrd_nn, lwork_zgelqf_mn, lwork_zgeqrf_mn : int = ---;
    lwork_zungbr_p_mn, lwork_zungbr_p_nn, lwork_zungbr_q_mn, lwork_zungbr_q_mm : int = ---;
    lwork_zunglq_mn, lwork_zunglq_nn, lwork_zungqr_mm, lwork_zungqr_mn : int = ---;
    lwork_zunmbr_prc_mm, lwork_zunmbr_prc_mn, lwork_zunmbr_prc_nn : int = ---;
    lwork_zunmbr_qln_mm, lwork_zunmbr_qln_mn, lwork_zunmbr_qln_nn : int = ---;

    info.* = 0;
    minmn = min(m, n);
    mnthr1 = cast(int) (ncast(RT, minmn) * ncast(RT, 17) / ncast(RT, 9));
    mnthr2 = cast(int) (ncast(RT, minmn) * ncast(RT, 5) / ncast(RT, 3));
    wntqa = char_eq(jobz, #char "a");
    wntqs = char_eq(jobz, #char "s");
    wntqas = wntqa || wntqs;
    wntqo = char_eq(jobz, #char "o");
    wntqn = char_eq(jobz, #char "n");
    lquery = lwork == -1;
    minwrk = 1;
    maxwrk = 1;

    if !(wntqa || wntqs || wntqo || wntqn) {
        info.* = -1;
    } else if m < 0 {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, m) {
        info.* = -5;
    } else if ldu < 1 || (wntqas && ldu < m) || (wntqo && m < n && ldu < m) {
        info.* = -8;
    } else if ldvt < 1 || (wntqa && ldvt < n) || (wntqs && ldvt < minmn) || (wntqo && m >= n && ldvt < n) {
        info.* = -10;
    }

    if info.* == 0 {
        if m >= n && minmn > 0 {
            lwork_zgebrd_mn = _gesdd_complex_query_gebrd(m, n, RT, MT);
            lwork_zgebrd_nn = _gesdd_complex_query_gebrd(n, n, RT, MT);
            lwork_zgeqrf_mn = _gesdd_complex_query_geqrf(m, n, RT, MT);
            lwork_zungbr_p_nn = _gesdd_complex_query_ungbr(#char "p", n, n, n, RT, MT);
            lwork_zungbr_q_mm = _gesdd_complex_query_ungbr(#char "q", m, m, n, RT, MT);
            lwork_zungbr_q_mn = _gesdd_complex_query_ungbr(#char "q", m, n, n, RT, MT);
            lwork_zungqr_mm = _gesdd_complex_query_ungqr(m, m, n, RT, MT);
            lwork_zungqr_mn = _gesdd_complex_query_ungqr(m, n, n, RT, MT);
            lwork_zunmbr_prc_nn = _gesdd_complex_query_unmbr(#char "p", #char "r", #char "c", n, n, n, RT, MT);
            lwork_zunmbr_qln_mm = _gesdd_complex_query_unmbr(#char "q", #char "l", #char "n", m, m, n, RT, MT);
            lwork_zunmbr_qln_mn = _gesdd_complex_query_unmbr(#char "q", #char "l", #char "n", m, n, n, RT, MT);
            lwork_zunmbr_qln_nn = _gesdd_complex_query_unmbr(#char "q", #char "l", #char "n", n, n, n, RT, MT);

            if m >= mnthr1 {
                if wntqn {
                    maxwrk = n + lwork_zgeqrf_mn;
                    maxwrk = max(maxwrk, 2*n + lwork_zgebrd_nn);
                    minwrk = 3*n;
                } else if wntqo {
                    wrkbl = n + lwork_zgeqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_zungqr_mn);
                    wrkbl = max(wrkbl, 2*n + lwork_zgebrd_nn);
                    wrkbl = max(wrkbl, 2*n + lwork_zunmbr_qln_nn);
                    wrkbl = max(wrkbl, 2*n + lwork_zunmbr_prc_nn);
                    maxwrk = m*n + n*n + wrkbl;
                    minwrk = 2*n*n + 3*n;
                } else if wntqs {
                    wrkbl = n + lwork_zgeqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_zungqr_mn);
                    wrkbl = max(wrkbl, 2*n + lwork_zgebrd_nn);
                    wrkbl = max(wrkbl, 2*n + lwork_zunmbr_qln_nn);
                    wrkbl = max(wrkbl, 2*n + lwork_zunmbr_prc_nn);
                    maxwrk = n*n + wrkbl;
                    minwrk = n*n + 3*n;
                } else {
                    wrkbl = n + lwork_zgeqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_zungqr_mm);
                    wrkbl = max(wrkbl, 2*n + lwork_zgebrd_nn);
                    wrkbl = max(wrkbl, 2*n + lwork_zunmbr_qln_nn);
                    wrkbl = max(wrkbl, 2*n + lwork_zunmbr_prc_nn);
                    maxwrk = n*n + wrkbl;
                    minwrk = n*n + max(3*n, n + m);
                }
            } else if m >= mnthr2 {
                maxwrk = 2*n + lwork_zgebrd_mn;
                minwrk = 2*n + m;
                if wntqo {
                    maxwrk = max(maxwrk, 2*n + lwork_zungbr_p_nn);
                    maxwrk = max(maxwrk, 2*n + lwork_zungbr_q_mn);
                    maxwrk += m*n;
                    minwrk += n*n;
                } else if wntqs {
                    maxwrk = max(maxwrk, 2*n + lwork_zungbr_p_nn);
                    maxwrk = max(maxwrk, 2*n + lwork_zungbr_q_mn);
                } else if wntqa {
                    maxwrk = max(maxwrk, 2*n + lwork_zungbr_p_nn);
                    maxwrk = max(maxwrk, 2*n + lwork_zungbr_q_mm);
                }
            } else {
                maxwrk = 2*n + lwork_zgebrd_mn;
                minwrk = 2*n + m;
                if wntqo {
                    maxwrk = max(maxwrk, 2*n + lwork_zunmbr_prc_nn);
                    maxwrk = max(maxwrk, 2*n + lwork_zunmbr_qln_mn);
                    maxwrk += m*n;
                    minwrk += n*n;
                } else if wntqs {
                    maxwrk = max(maxwrk, 2*n + lwork_zunmbr_qln_mn);
                    maxwrk = max(maxwrk, 2*n + lwork_zunmbr_prc_nn);
                } else if wntqa {
                    maxwrk = max(maxwrk, 2*n + lwork_zunmbr_qln_mm);
                    maxwrk = max(maxwrk, 2*n + lwork_zunmbr_prc_nn);
                }
            }
        } else if minmn > 0 {
            lwork_zgebrd_mn = _gesdd_complex_query_gebrd(m, n, RT, MT);
            lwork_zgebrd_mm = _gesdd_complex_query_gebrd(m, m, RT, MT);
            lwork_zgelqf_mn = _gesdd_complex_query_gelqf(m, n, RT, MT);
            lwork_zungbr_p_mn = _gesdd_complex_query_ungbr(#char "p", m, n, m, RT, MT);
            lwork_zungbr_p_nn = _gesdd_complex_query_ungbr(#char "p", n, n, m, RT, MT);
            lwork_zungbr_q_mm = _gesdd_complex_query_ungbr(#char "q", m, m, n, RT, MT);
            lwork_zunglq_mn = _gesdd_complex_query_unglq(m, n, m, RT, MT);
            lwork_zunglq_nn = _gesdd_complex_query_unglq(n, n, m, RT, MT);
            lwork_zunmbr_prc_mm = _gesdd_complex_query_unmbr(#char "p", #char "r", #char "c", m, m, m, RT, MT);
            lwork_zunmbr_prc_mn = _gesdd_complex_query_unmbr(#char "p", #char "r", #char "c", m, n, m, RT, MT);
            lwork_zunmbr_prc_nn = _gesdd_complex_query_unmbr(#char "p", #char "r", #char "c", n, n, m, RT, MT);
            lwork_zunmbr_qln_mm = _gesdd_complex_query_unmbr(#char "q", #char "l", #char "n", m, m, m, RT, MT);

            if n >= mnthr1 {
                if wntqn {
                    maxwrk = m + lwork_zgelqf_mn;
                    maxwrk = max(maxwrk, 2*m + lwork_zgebrd_mm);
                    minwrk = 3*m;
                } else if wntqo {
                    wrkbl = m + lwork_zgelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_zunglq_mn);
                    wrkbl = max(wrkbl, 2*m + lwork_zgebrd_mm);
                    wrkbl = max(wrkbl, 2*m + lwork_zunmbr_qln_mm);
                    wrkbl = max(wrkbl, 2*m + lwork_zunmbr_prc_mm);
                    maxwrk = m*n + m*m + wrkbl;
                    minwrk = 2*m*m + 3*m;
                } else if wntqs {
                    wrkbl = m + lwork_zgelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_zunglq_mn);
                    wrkbl = max(wrkbl, 2*m + lwork_zgebrd_mm);
                    wrkbl = max(wrkbl, 2*m + lwork_zunmbr_qln_mm);
                    wrkbl = max(wrkbl, 2*m + lwork_zunmbr_prc_mm);
                    maxwrk = m*m + wrkbl;
                    minwrk = m*m + 3*m;
                } else {
                    wrkbl = m + lwork_zgelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_zunglq_nn);
                    wrkbl = max(wrkbl, 2*m + lwork_zgebrd_mm);
                    wrkbl = max(wrkbl, 2*m + lwork_zunmbr_qln_mm);
                    wrkbl = max(wrkbl, 2*m + lwork_zunmbr_prc_mm);
                    maxwrk = m*m + wrkbl;
                    minwrk = m*m + max(3*m, m + n);
                }
            } else if n >= mnthr2 {
                maxwrk = 2*m + lwork_zgebrd_mn;
                minwrk = 2*m + n;
                if wntqo {
                    maxwrk = max(maxwrk, 2*m + lwork_zungbr_q_mm);
                    maxwrk = max(maxwrk, 2*m + lwork_zungbr_p_mn);
                    maxwrk += m*n;
                    minwrk += m*m;
                } else if wntqs {
                    maxwrk = max(maxwrk, 2*m + lwork_zungbr_q_mm);
                    maxwrk = max(maxwrk, 2*m + lwork_zungbr_p_mn);
                } else if wntqa {
                    maxwrk = max(maxwrk, 2*m + lwork_zungbr_q_mm);
                    maxwrk = max(maxwrk, 2*m + lwork_zungbr_p_nn);
                }
            } else {
                maxwrk = 2*m + lwork_zgebrd_mn;
                minwrk = 2*m + n;
                if wntqo {
                    maxwrk = max(maxwrk, 2*m + lwork_zunmbr_qln_mm);
                    maxwrk = max(maxwrk, 2*m + lwork_zunmbr_prc_mn);
                    maxwrk += m*n;
                    minwrk += m*m;
                } else if wntqs {
                    maxwrk = max(maxwrk, 2*m + lwork_zunmbr_qln_mm);
                    maxwrk = max(maxwrk, 2*m + lwork_zunmbr_prc_mn);
                } else if wntqa {
                    maxwrk = max(maxwrk, 2*m + lwork_zunmbr_qln_mm);
                    maxwrk = max(maxwrk, 2*m + lwork_zunmbr_prc_nn);
                }
            }
        }

        maxwrk = max(maxwrk, minwrk);
        _gesdd_set_work1(work, _gesdd_roundup_lwork(ncast(RT, maxwrk), RT), RT, MT);
        if lwork < minwrk && !lquery {
            info.* = -12;
        }
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if m == 0 || n == 0 {
        return;
    }
}
