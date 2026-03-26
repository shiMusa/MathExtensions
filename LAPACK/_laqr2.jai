// ? LAPACK *laqr2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqr2/dlaqr2/claqr2/zlaqr2
// *LAQR2 performs aggressive early deflation for the QR algorithm.

slaqr2 :: #bake_arguments _laqr2_real(RT=float32);
dlaqr2 :: #bake_arguments _laqr2_real(RT=float64);
claqr2 :: #bake_arguments _laqr2_complex(RT=float32, MT=C32);
zlaqr2 :: #bake_arguments _laqr2_complex(RT=float64, MT=C64);

#scope_file
_laqr2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_laqr2_storage :: inline (data: *$T, $$count: int, $$inc: int) -> [] T {
    storage_count := ifx count <= 0 then 0 else 1 + (count - 1)*inc;
    return _laqr2_raw_slice(data, storage_count);
}

_laqr2_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _laqr2_raw_slice(x.data + (offset - 1), count);
}

_laqr2_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else lda*(cols - 1) + rows;
    return _laqr2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqr2_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_laqr2_query_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_laqr2_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laqr2_abs1 :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(x);
    } else {
        return abs(real(x)) + abs(imag(x));
    }
}

_laqr2_set_work1_real :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && WORK.data_type == RT;
} { fset(work, 1, x); }

_laqr2_set_work1_complex :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} { fset(work, 1, cpx(x, zero(RT))); }

_laqr2_copy :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        scopy(n, x, incx, y, incy);
    } else if MT == float64 {
        dcopy(n, x, incx, y, incy);
    } else if MT == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_laqr2_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT;
} {
    #if MT == float32 {
        slacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == float64 {
        dlacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_laqr2_laset :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $MT, beta: MT, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT;
} {
    #if MT == float32 {
        slaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == float64 {
        dlaset(uplo, m, n, alpha, beta, a, lda);
    } else if MT == C32 {
        claset(uplo, m, n, alpha, beta, a, lda);
    } else {
        zlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_laqr2_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int, alpha: $MT,
    a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int, beta: MT,
    c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == float32 {
        sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == float64 {
        dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_laqr2_larfg_real :: ($$n: int, alpha: *$RT, x: *$X/ArrayType1D, $$incx: int, tau: *RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        slarfg(n, alpha, x, incx, tau);
    } else {
        dlarfg(n, alpha, x, incx, tau);
    }
}

_laqr2_larfg_complex :: ($$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *MT, $MT: Type) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        clarfg(n, alpha, x, incx, tau);
    } else {
        zlarfg(n, alpha, x, incx, tau);
    }
}

_laqr2_larf1f_real :: (
    $$side: u8, $$m: int, $$n: int, v: $V/ArrayType1D, $$incv: int, $$tau: $RT, c: *$C/ArrayType2D, $$ldc: int,
    work: *$WORK/ArrayType1D, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           V.data_type == RT && C.data_type == RT && WORK.data_type == RT;
} {
    #if RT == float32 {
        slarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else {
        dlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    }
}

_laqr2_larf1f_complex :: (
    $$side: u8, $$m: int, $$n: int, v: $V/ArrayType1D, $$incv: int, tau: $MT, c: *$C/ArrayType2D, $$ldc: int,
    work: *$WORK/ArrayType1D, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           V.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        clarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else {
        zlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    }
}

_laqr2_gehrd_real :: (
    $$n: int, $$ilo: int, $$ihi: int, a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && TAU.data_type == RT && WORK.data_type == RT;
} {
    #if RT == float32 {
        sgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info);
    } else {
        dgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info);
    }
}

_laqr2_gehrd_complex :: (
    $$n: int, $$ilo: int, $$ihi: int, a: *$A/ArrayType2D, $$lda: int, tau: *$TAU/ArrayType1D,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && TAU.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        cgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info);
    } else {
        zgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info);
    }
}

_laqr2_ormhr_real :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$ilo: int, $$ihi: int, a: $A/ArrayType2D, $$lda: int,
    tau: $TAU/ArrayType1D, c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT && TAU.data_type == RT && C.data_type == RT && WORK.data_type == RT;
} {
    #if RT == float32 {
        sormhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        dormhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_laqr2_unmhr_complex :: (
    $$side: u8, $$trans: u8, $$m: int, $$n: int, $$ilo: int, $$ihi: int, a: $A/ArrayType2D, $$lda: int,
    tau: $TAU/ArrayType1D, c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && TAU.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        cunmhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info);
    } else {
        zunmhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info);
    }
}

_laqr2_lahqr_real :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int, h: *$H/ArrayType2D, $$ldh: int,
    wr: *$WR/ArrayType1D, wi: *$WI/ArrayType1D, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT && WR.data_type == RT && WI.data_type == RT && Z.data_type == RT;
} {
    #if RT == float32 {
        slahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info);
    } else {
        dlahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info);
    }
}

_laqr2_lahqr_complex :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int, h: *$H/ArrayType2D, $$ldh: int,
    w: *$W/ArrayType1D, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           H.data_type == MT && W.data_type == MT && Z.data_type == MT;
} {
    #if MT == C32 {
        clahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
    } else {
        zlahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
    }
}

_laqr2_trexc_real :: (
    $$compq: u8, $$n: int, t: *$T/ArrayType2D, $$ldt: int, q: *$Q/ArrayType2D, $$ldq: int,
    ifst: *int, ilst: *int, work: *$WORK/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           T.data_type == RT && Q.data_type == RT && WORK.data_type == RT;
} {
    #if RT == float32 {
        strexc(compq, n, t, ldt, q, ldq, ifst, ilst, work, info);
    } else {
        dtrexc(compq, n, t, ldt, q, ldq, ifst, ilst, work, info);
    }
}

_laqr2_trexc_complex :: (
    $$compq: u8, $$n: int, t: *$T/ArrayType2D, $$ldt: int, q: *$Q/ArrayType2D, $$ldq: int,
    ifst: *int, ilst: *int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           T.data_type == MT && Q.data_type == MT;
} {
    #if MT == C32 {
        ctrexc(compq, n, t, ldt, q, ldq, ifst, ilst, info);
    } else {
        ztrexc(compq, n, t, ldt, q, ldq, ifst, ilst, info);
    }
}

_laqr2_lanv2 :: (
    a: *$RT, b: *RT, c: *RT, d: *RT, rt1r: *RT, rt1i: *RT, rt2r: *RT, rt2i: *RT, cs: *RT, sn: *RT, $RT: Type
) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        slanv2(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn);
    } else {
        dlanv2(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn);
    }
}

_laqr2_real :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ktop: int, $$kbot: int, $$nw: int,
    h: *$H/ArrayType2D, $$ldh: int, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    ns: *int, nd: *int, sr: *$SR/ArrayType1D, si: *$SI/ArrayType1D,
    v: *$V/ArrayType2D, $$ldv: int, $$nh: int, t: *$T/ArrayType2D, $$ldt: int,
    $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int, work: *$WORK/ArrayType1D, $$lwork: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT && Z.data_type == RT && SR.data_type == RT && SI.data_type == RT &&
           V.data_type == RT && T.data_type == RT && WV.data_type == RT && WORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);

    aa, bb, cc, cs, dd, evi, evk, foo, rt1i, rt1r, rt2i, rt2r, s, safmax, safmin, smlnum, sn, tau, ulp : RT = ---;
    i, ifst, ilst, info_local, infqr, j, jw, k, kend, kcol, kln, krow, kwtop, ltop, lwk1, lwk2, lwkopt : int = ---;
    bulge, sorted : bool = ---;

    jw = min(nw, kbot - ktop + 1);
    if jw <= 2 {
        lwkopt = 1;
    } else {
        tau_view := array1d_view(_laqr2_subvector(work.*, 1, work.count));
        work_view := array1d_view(_laqr2_subvector(work.*, 1, work.count));
        _laqr2_gehrd_real(jw, 1, jw - 1, t, ldt, *tau_view, *work_view, -1, info_local, RT);
        lwk1 = _laqr2_query_to_int(fget(work.*, 1), RT, RT);

        _laqr2_ormhr_real(#char "R", #char "N", jw, jw, 1, jw - 1, t.*, ldt, *tau_view, v, ldv, *work_view, -1, info_local, RT);
        lwk2 = _laqr2_query_to_int(fget(work.*, 1), RT, RT);
        lwkopt = jw + max(lwk1, lwk2);
    }

    if lwork == -1 {
        _laqr2_set_work1_real(work, _laqr2_roundup_lwork(lwkopt, RT), RT);
        return;
    }

    ns.* = 0;
    nd.* = 0;
    _laqr2_set_work1_real(work, one_r, RT);
    if ktop > kbot || nw < 1 {
        return;
    }

    safmin = ifx RT == float32 then slamch("Safe minimum".data) else dlamch("Safe minimum".data);
    safmax = one_r / safmin;
    ulp = ifx RT == float32 then slamch("Precision".data) else dlamch("Precision".data);
    smlnum = safmin*(ncast(RT, n) / ulp);
    _ = safmax;

    jw = min(nw, kbot - ktop + 1);
    kwtop = kbot - jw + 1;
    if kwtop == ktop {
        s = zero_r;
    } else {
        s = fget(h.*, kwtop, kwtop - 1);
    }

    if kbot == kwtop {
        fset(sr, kwtop, fget(h.*, kwtop, kwtop));
        fset(si, kwtop, zero_r);
        ns.* = 1;
        nd.* = 0;
        if abs(s) <= max(smlnum, ulp*abs(fget(h.*, kwtop, kwtop))) {
            ns.* = 0;
            nd.* = 1;
            if kwtop > ktop {
                fset(h, kwtop, kwtop - 1, zero_r);
            }
        }
        _laqr2_set_work1_real(work, one_r, RT);
        return;
    }

    hwin_view := array2d_view(_laqr2_matrix_storage(h.*, ldh, kwtop, kwtop, jw, jw), ldh, row_major=false);
    _laqr2_lacpy(#char "U", jw, jw, *hwin_view, ldh, t, ldt, RT);
    _laqr2_copy(jw - 1, _laqr2_storage(h.data + (kwtop - 1)*ldh + kwtop, jw - 1, ldh + 1), ldh + 1,
                _laqr2_storage(t.data + 1, jw - 1, ldt + 1), ldt + 1, RT);

    _laqr2_laset(#char "A", jw, jw, zero_r, one_r, v, ldv, RT);
    sr_win := array1d_view(_laqr2_subvector(sr.*, kwtop, sr.count - kwtop + 1));
    si_win := array1d_view(_laqr2_subvector(si.*, kwtop, si.count - kwtop + 1));
    _laqr2_lahqr_real(true, true, jw, 1, jw, t, ldt, *sr_win, *si_win, 1, jw, v, ldv, infqr, RT);

    for j_loop:1..jw - 3 {
        j = j_loop;
        fset(t, j + 2, j, zero_r);
        fset(t, j + 3, j, zero_r);
    }
    if jw > 2 {
        fset(t, jw, jw - 2, zero_r);
    }

    ns.* = jw;
    ilst = infqr + 1;
    while ilst <= ns.* {
        if ns.* == 1 {
            bulge = false;
        } else {
            bulge = fget(t.*, ns.*, ns.* - 1) != zero_r;
        }

        if !bulge {
            foo = abs(fget(t.*, ns.*, ns.*));
            if foo == zero_r {
                foo = abs(s);
            }
            if abs(s*fget(v.*, 1, ns.*)) <= max(smlnum, ulp*foo) {
                ns.* -= 1;
            } else {
                ifst = ns.*;
                work_vec := array1d_view(_laqr2_subvector(work.*, 1, work.count));
                _laqr2_trexc_real(#char "V", jw, t, ldt, v, ldv, ifst, ilst, *work_vec, info_local, RT);
                ilst += 1;
            }
        } else {
            foo = abs(fget(t.*, ns.*, ns.*)) + sqrt(abs(fget(t.*, ns.*, ns.* - 1)))*sqrt(abs(fget(t.*, ns.* - 1, ns.*)));
            if foo == zero_r {
                foo = abs(s);
            }
            if max(abs(s*fget(v.*, 1, ns.*)), abs(s*fget(v.*, 1, ns.* - 1))) <= max(smlnum, ulp*foo) {
                ns.* -= 2;
            } else {
                ifst = ns.*;
                work_vec := array1d_view(_laqr2_subvector(work.*, 1, work.count));
                _laqr2_trexc_real(#char "V", jw, t, ldt, v, ldv, ifst, ilst, *work_vec, info_local, RT);
                ilst += 2;
            }
        }
    }

    if ns.* == 0 {
        s = zero_r;
    }

    if ns.* < jw {
        sorted = false;
        i = ns.* + 1;
        while !sorted {
            sorted = true;
            kend = i - 1;
            i = infqr + 1;
            if i == ns.* {
                k = i + 1;
            } else if fget(t.*, i + 1, i) == zero_r {
                k = i + 1;
            } else {
                k = i + 2;
            }

            while k <= kend {
                if k == i + 1 {
                    evi = abs(fget(t.*, i, i));
                } else {
                    evi = abs(fget(t.*, i, i)) + sqrt(abs(fget(t.*, i + 1, i)))*sqrt(abs(fget(t.*, i, i + 1)));
                }

                if k == kend {
                    evk = abs(fget(t.*, k, k));
                } else if fget(t.*, k + 1, k) == zero_r {
                    evk = abs(fget(t.*, k, k));
                } else {
                    evk = abs(fget(t.*, k, k)) + sqrt(abs(fget(t.*, k + 1, k)))*sqrt(abs(fget(t.*, k, k + 1)));
                }

                if evi >= evk {
                    i = k;
                } else {
                    sorted = false;
                    ifst = i;
                    ilst = k;
                    work_vec := array1d_view(_laqr2_subvector(work.*, 1, work.count));
                    _laqr2_trexc_real(#char "V", jw, t, ldt, v, ldv, ifst, ilst, *work_vec, info_local, RT);
                    if info_local == 0 {
                        i = ilst;
                    } else {
                        i = k;
                    }
                }

                if i == kend {
                    k = i + 1;
                } else if fget(t.*, i + 1, i) == zero_r {
                    k = i + 1;
                } else {
                    k = i + 2;
                }
            }
        }
    }

    i = jw;
    while i >= infqr + 1 {
        if i == infqr + 1 {
            fset(sr, kwtop + i - 1, fget(t.*, i, i));
            fset(si, kwtop + i - 1, zero_r);
            i -= 1;
        } else if fget(t.*, i, i - 1) == zero_r {
            fset(sr, kwtop + i - 1, fget(t.*, i, i));
            fset(si, kwtop + i - 1, zero_r);
            i -= 1;
        } else {
            aa = fget(t.*, i - 1, i - 1);
            cc = fget(t.*, i, i - 1);
            bb = fget(t.*, i - 1, i);
            dd = fget(t.*, i, i);
            _laqr2_lanv2(aa, bb, cc, dd, rt1r, rt1i, rt2r, rt2i, cs, sn, RT);
            fset(sr, kwtop + i - 2, rt1r);
            fset(si, kwtop + i - 2, rt1i);
            fset(sr, kwtop + i - 1, rt2r);
            fset(si, kwtop + i - 1, rt2i);
            i -= 2;
        }
    }

    if ns.* < jw || s == zero_r {
        if ns.* > 1 && s != zero_r {
            work_vec := array1d_view(_laqr2_subvector(work.*, 1, work.count));
            _laqr2_copy(ns.*, _laqr2_storage(v.data, ns.*, ldv), ldv, *work_vec, 1, RT);
            work_tail := array1d_view(_laqr2_subvector(work.*, 2, max(0, work.count - 1)));
            _laqr2_larfg_real(ns.*, *work.data, *work_tail, 1, tau, RT);

            t_lower_view := array2d_view(_laqr2_matrix_storage(t.*, ldt, 3, 1, jw - 2, jw - 2), ldt, row_major=false);
            _laqr2_laset(#char "L", jw - 2, jw - 2, zero_r, zero_r, *t_lower_view, ldt, RT);

            scratch_view := array1d_view(_laqr2_subvector(work.*, jw + 1, max(0, work.count - jw)));
            _laqr2_larf1f_real(#char "L", ns.*, jw, *work_vec, 1, tau, t, ldt, *scratch_view, RT);
            _laqr2_larf1f_real(#char "R", ns.*, ns.*, *work_vec, 1, tau, t, ldt, *scratch_view, RT);
            _laqr2_larf1f_real(#char "R", jw, ns.*, *work_vec, 1, tau, v, ldv, *scratch_view, RT);

            tau_hess := array1d_view(_laqr2_subvector(work.*, 1, work.count));
            _laqr2_gehrd_real(jw, 1, ns.*, t, ldt, *tau_hess, *scratch_view, lwork - jw, info_local, RT);
        }

        if kwtop > 1 {
            fset(h, kwtop, kwtop - 1, s*fget(v.*, 1, 1));
        }
        _laqr2_lacpy(#char "U", jw, jw, t.*, ldt, *hwin_view, ldh, RT);
        _laqr2_copy(jw - 1, _laqr2_storage(t.data + 1, jw - 1, ldt + 1), ldt + 1,
                    _laqr2_storage(h.data + (kwtop - 1)*ldh + kwtop, jw - 1, ldh + 1), ldh + 1, RT);

        if ns.* > 1 && s != zero_r {
            tau_hess := array1d_view(_laqr2_subvector(work.*, 1, work.count));
            scratch_view := array1d_view(_laqr2_subvector(work.*, jw + 1, max(0, work.count - jw)));
            _laqr2_ormhr_real(#char "R", #char "N", jw, ns.*, 1, ns.*, t.*, ldt, *tau_hess, v, ldv, *scratch_view, lwork - jw, info_local, RT);
        }

        if wantt {
            ltop = 1;
        } else {
            ltop = ktop;
        }
        krow = ltop;
        while krow <= kwtop - 1 {
            kln = min(nv, kwtop - krow);
            h_vert := array2d_view(_laqr2_matrix_storage(h.*, ldh, krow, kwtop, kln, jw), ldh, row_major=false);
            wv_vert := array2d_view(_laqr2_matrix_storage(wv.*, ldwv, 1, 1, kln, jw), ldwv, row_major=false);
            _laqr2_gemm(#char "N", #char "N", kln, jw, jw, one_r, *h_vert, ldh, v.*, ldv, zero_r, *wv_vert, ldwv, RT);
            _laqr2_lacpy(#char "A", kln, jw, *wv_vert, ldwv, *h_vert, ldh, RT);
            krow += nv;
        }

        if wantt {
            kcol = kbot + 1;
            while kcol <= n {
                kln = min(nh, n - kcol + 1);
                h_horz := array2d_view(_laqr2_matrix_storage(h.*, ldh, kwtop, kcol, jw, kln), ldh, row_major=false);
                t_horz := array2d_view(_laqr2_matrix_storage(t.*, ldt, 1, 1, jw, kln), ldt, row_major=false);
                _laqr2_gemm(#char "T", #char "N", jw, kln, jw, one_r, v.*, ldv, *h_horz, ldh, zero_r, *t_horz, ldt, RT);
                _laqr2_lacpy(#char "A", jw, kln, *t_horz, ldt, *h_horz, ldh, RT);
                kcol += nh;
            }
        }

        if wantz {
            krow = iloz;
            while krow <= ihiz {
                kln = min(nv, ihiz - krow + 1);
                z_vert := array2d_view(_laqr2_matrix_storage(z.*, ldz, krow, kwtop, kln, jw), ldz, row_major=false);
                wv_vert := array2d_view(_laqr2_matrix_storage(wv.*, ldwv, 1, 1, kln, jw), ldwv, row_major=false);
                _laqr2_gemm(#char "N", #char "N", kln, jw, jw, one_r, *z_vert, ldz, v.*, ldv, zero_r, *wv_vert, ldwv, RT);
                _laqr2_lacpy(#char "A", kln, jw, *wv_vert, ldwv, *z_vert, ldz, RT);
                krow += nv;
            }
        }
    }

    nd.* = jw - ns.*;
    ns.* = ns.* - infqr;
    _laqr2_set_work1_real(work, _laqr2_roundup_lwork(lwkopt, RT), RT);
}

_laqr2_complex :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ktop: int, $$kbot: int, $$nw: int,
    h: *$H/ArrayType2D, $$ldh: int, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    ns: *int, nd: *int, sh: *$SH/ArrayType1D,
    v: *$V/ArrayType2D, $$ldv: int, $$nh: int, t: *$T/ArrayType2D, $$ldt: int,
    $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int, work: *$WORK/ArrayType1D, $$lwork: int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT && Z.data_type == MT && SH.data_type == MT &&
           V.data_type == MT && T.data_type == MT && WV.data_type == MT && WORK.data_type == MT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    zero_t :: #run zero(MT);
    one_t  :: #run _laqr2_complex_from_real(one_r, RT, MT);

    cdum, s, tau : MT = ---;
    foo, safmax, safmin, smlnum, ulp : RT = ---;
    i, ifst, ilst, info_local, infqr, j, jw, kcol, kln, knt, krow, kwtop, ltop, lwk1, lwk2, lwkopt : int = ---;

    jw = min(nw, kbot - ktop + 1);
    if jw <= 2 {
        lwkopt = 1;
    } else {
        tau_view := array1d_view(_laqr2_subvector(work.*, 1, work.count));
        work_view := array1d_view(_laqr2_subvector(work.*, 1, work.count));
        _laqr2_gehrd_complex(jw, 1, jw - 1, t, ldt, *tau_view, *work_view, -1, info_local, MT);
        lwk1 = _laqr2_query_to_int(fget(work.*, 1), RT, MT);

        _laqr2_unmhr_complex(#char "R", #char "N", jw, jw, 1, jw - 1, t.*, ldt, *tau_view, v, ldv, *work_view, -1, info_local, MT);
        lwk2 = _laqr2_query_to_int(fget(work.*, 1), RT, MT);
        lwkopt = jw + max(lwk1, lwk2);
    }

    if lwork == -1 {
        _laqr2_set_work1_complex(work, ncast(RT, lwkopt), RT, MT);
        return;
    }

    ns.* = 0;
    nd.* = 0;
    _laqr2_set_work1_complex(work, one_r, RT, MT);
    if ktop > kbot || nw < 1 {
        return;
    }

    safmin = ifx RT == float32 then slamch("Safe minimum".data) else dlamch("Safe minimum".data);
    safmax = one_r / safmin;
    ulp = ifx RT == float32 then slamch("Precision".data) else dlamch("Precision".data);
    smlnum = safmin*(ncast(RT, n) / ulp);
    _ = safmax;

    jw = min(nw, kbot - ktop + 1);
    kwtop = kbot - jw + 1;
    if kwtop == ktop {
        s = zero_t;
    } else {
        s = fget(h.*, kwtop, kwtop - 1);
    }

    if kbot == kwtop {
        fset(sh, kwtop, fget(h.*, kwtop, kwtop));
        ns.* = 1;
        nd.* = 0;
        if _laqr2_abs1(s, RT, MT) <= max(smlnum, ulp*_laqr2_abs1(fget(h.*, kwtop, kwtop), RT, MT)) {
            ns.* = 0;
            nd.* = 1;
            if kwtop > ktop {
                fset(h, kwtop, kwtop - 1, zero_t);
            }
        }
        _laqr2_set_work1_complex(work, one_r, RT, MT);
        return;
    }

    hwin_view := array2d_view(_laqr2_matrix_storage(h.*, ldh, kwtop, kwtop, jw, jw), ldh, row_major=false);
    _laqr2_lacpy(#char "U", jw, jw, *hwin_view, ldh, t, ldt, MT);
    _laqr2_copy(jw - 1, _laqr2_storage(h.data + (kwtop - 1)*ldh + kwtop, jw - 1, ldh + 1), ldh + 1,
                _laqr2_storage(t.data + 1, jw - 1, ldt + 1), ldt + 1, MT);

    _laqr2_laset(#char "A", jw, jw, zero_t, one_t, v, ldv, MT);
    sh_win := array1d_view(_laqr2_subvector(sh.*, kwtop, sh.count - kwtop + 1));
    _laqr2_lahqr_complex(true, true, jw, 1, jw, t, ldt, *sh_win, 1, jw, v, ldv, infqr, MT);

    ns.* = jw;
    ilst = infqr + 1;
    for knt:infqr + 1..jw {
        foo = _laqr2_abs1(fget(t.*, ns.*, ns.*), RT, MT);
        if foo == zero_r {
            foo = _laqr2_abs1(s, RT, MT);
        }
        if _laqr2_abs1(s, RT, MT)*_laqr2_abs1(fget(v.*, 1, ns.*), RT, MT) <= max(smlnum, ulp*foo) {
            ns.* -= 1;
        } else {
            ifst = ns.*;
            _laqr2_trexc_complex(#char "V", jw, t, ldt, v, ldv, ifst, ilst, info_local, MT);
            ilst += 1;
        }
    }

    if ns.* == 0 {
        s = zero_t;
    }

    if ns.* < jw {
        for i_loop:infqr + 1..ns.* {
            i = i_loop;
            ifst = i;
            for j_loop:i + 1..ns.* {
                j = j_loop;
                if _laqr2_abs1(fget(t.*, j, j), RT, MT) > _laqr2_abs1(fget(t.*, ifst, ifst), RT, MT) {
                    ifst = j;
                }
            }
            ilst = i;
            if ifst != ilst {
                _laqr2_trexc_complex(#char "V", jw, t, ldt, v, ldv, ifst, ilst, info_local, MT);
            }
        }
    }

    for i_loop:infqr + 1..jw {
        i = i_loop;
        fset(sh, kwtop + i - 1, fget(t.*, i, i));
    }

    if ns.* < jw || s == zero_t {
        if ns.* > 1 && s != zero_t {
            work_vec := array1d_view(_laqr2_subvector(work.*, 1, work.count));
            _laqr2_copy(ns.*, _laqr2_storage(v.data, ns.*, ldv), ldv, *work_vec, 1, MT);
            for i_loop:1..ns.* {
                i = i_loop;
                fset(work, i, conjugate(fget(work.*, i)));
            }
            work_tail := array1d_view(_laqr2_subvector(work.*, 2, max(0, work.count - 1)));
            _laqr2_larfg_complex(ns.*, *work.data, *work_tail, 1, tau, MT);

            t_lower_view := array2d_view(_laqr2_matrix_storage(t.*, ldt, 3, 1, jw - 2, jw - 2), ldt, row_major=false);
            _laqr2_laset(#char "L", jw - 2, jw - 2, zero_t, zero_t, *t_lower_view, ldt, MT);

            scratch_view := array1d_view(_laqr2_subvector(work.*, jw + 1, max(0, work.count - jw)));
            _laqr2_larf1f_complex(#char "L", ns.*, jw, *work_vec, 1, conjugate(tau), t, ldt, *scratch_view, MT);
            _laqr2_larf1f_complex(#char "R", ns.*, ns.*, *work_vec, 1, tau, t, ldt, *scratch_view, MT);
            _laqr2_larf1f_complex(#char "R", jw, ns.*, *work_vec, 1, tau, v, ldv, *scratch_view, MT);

            tau_hess := array1d_view(_laqr2_subvector(work.*, 1, work.count));
            _laqr2_gehrd_complex(jw, 1, ns.*, t, ldt, *tau_hess, *scratch_view, lwork - jw, info_local, MT);
        }

        if kwtop > 1 {
            fset(h, kwtop, kwtop - 1, s*conjugate(fget(v.*, 1, 1)));
        }
        _laqr2_lacpy(#char "U", jw, jw, t.*, ldt, *hwin_view, ldh, MT);
        _laqr2_copy(jw - 1, _laqr2_storage(t.data + 1, jw - 1, ldt + 1), ldt + 1,
                    _laqr2_storage(h.data + (kwtop - 1)*ldh + kwtop, jw - 1, ldh + 1), ldh + 1, MT);

        if ns.* > 1 && s != zero_t {
            tau_hess := array1d_view(_laqr2_subvector(work.*, 1, work.count));
            scratch_view := array1d_view(_laqr2_subvector(work.*, jw + 1, max(0, work.count - jw)));
            _laqr2_unmhr_complex(#char "R", #char "N", jw, ns.*, 1, ns.*, t.*, ldt, *tau_hess, v, ldv, *scratch_view, lwork - jw, info_local, MT);
        }

        if wantt {
            ltop = 1;
        } else {
            ltop = ktop;
        }
        krow = ltop;
        while krow <= kwtop - 1 {
            kln = min(nv, kwtop - krow);
            h_vert := array2d_view(_laqr2_matrix_storage(h.*, ldh, krow, kwtop, kln, jw), ldh, row_major=false);
            wv_vert := array2d_view(_laqr2_matrix_storage(wv.*, ldwv, 1, 1, kln, jw), ldwv, row_major=false);
            _laqr2_gemm(#char "N", #char "N", kln, jw, jw, one_t, *h_vert, ldh, v.*, ldv, zero_t, *wv_vert, ldwv, MT);
            _laqr2_lacpy(#char "A", kln, jw, *wv_vert, ldwv, *h_vert, ldh, MT);
            krow += nv;
        }

        if wantt {
            kcol = kbot + 1;
            while kcol <= n {
                kln = min(nh, n - kcol + 1);
                h_horz := array2d_view(_laqr2_matrix_storage(h.*, ldh, kwtop, kcol, jw, kln), ldh, row_major=false);
                t_horz := array2d_view(_laqr2_matrix_storage(t.*, ldt, 1, 1, jw, kln), ldt, row_major=false);
                _laqr2_gemm(#char "C", #char "N", jw, kln, jw, one_t, v.*, ldv, *h_horz, ldh, zero_t, *t_horz, ldt, MT);
                _laqr2_lacpy(#char "A", jw, kln, *t_horz, ldt, *h_horz, ldh, MT);
                kcol += nh;
            }
        }

        if wantz {
            krow = iloz;
            while krow <= ihiz {
                kln = min(nv, ihiz - krow + 1);
                z_vert := array2d_view(_laqr2_matrix_storage(z.*, ldz, krow, kwtop, kln, jw), ldz, row_major=false);
                wv_vert := array2d_view(_laqr2_matrix_storage(wv.*, ldwv, 1, 1, kln, jw), ldwv, row_major=false);
                _laqr2_gemm(#char "N", #char "N", kln, jw, jw, one_t, *z_vert, ldz, v.*, ldv, zero_t, *wv_vert, ldwv, MT);
                _laqr2_lacpy(#char "A", kln, jw, *wv_vert, ldwv, *z_vert, ldz, MT);
                krow += nv;
            }
        }
    }

    nd.* = jw - ns.*;
    ns.* = ns.* - infqr;
    _laqr2_set_work1_complex(work, ncast(RT, lwkopt), RT, MT);
}
