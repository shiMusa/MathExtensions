// ? LAPACK *laqr0()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqr0/dlaqr0/claqr0/zlaqr0
// *LAQR0 is the multishift QR driver used by *HSEQR.

slaqr0 :: #bake_arguments _laqr0_real(RT=float32);
dlaqr0 :: #bake_arguments _laqr0_real(RT=float64);
claqr0 :: #bake_arguments _laqr0_complex(RT=float32, MT=C32);
zlaqr0 :: #bake_arguments _laqr0_complex(RT=float64, MT=C64);

#scope_file
_laqr0_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_laqr0_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _laqr0_raw_slice(x.data + (offset - 1), count);
}

_laqr0_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else lda*(cols - 1) + rows;
    return _laqr0_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqr0_query_to_int :: inline (x: $MT, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return cast(int) x;
    } else {
        return cast(int) real(x);
    }
}

_laqr0_roundup_lwork :: inline ($$lwork: int, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sroundup_lwork(lwork);
    } else {
        return droundup_lwork(lwork);
    }
}

_laqr0_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laqr0_abs1_complex :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return abs(real(x)) + abs(imag(x));
}

_laqr0_set_work1_real :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && WORK.data_type == RT;
} { fset(work, 1, x); }

_laqr0_set_work1_complex :: inline (work: *$WORK/ArrayType1D, $$x: $RT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && WORK.data_type == MT;
} { fset(work, 1, cpx(x, zero(RT))); }

_laqr0_lacpy :: (
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

_laqr0_lahqr_real :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int, h: *$H/ArrayType2D, $$ldh: int,
    wr: *$WR/ArrayType1D, wi: *$WI/ArrayType1D, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    info: *int, $RT: Type
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

_laqr0_lahqr_complex :: (
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

_laqr0_laqr3_real :: (
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
    #if RT == float32 {
        slaqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
    } else {
        dlaqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
    }
}

_laqr0_laqr3_complex :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ktop: int, $$kbot: int, $$nw: int,
    h: *$H/ArrayType2D, $$ldh: int, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    ns: *int, nd: *int, w: *$W/ArrayType1D,
    v: *$V/ArrayType2D, $$ldv: int, $$nh: int, t: *$T/ArrayType2D, $$ldt: int,
    $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int, work: *$WORK/ArrayType1D, $$lwork: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           H.data_type == MT && Z.data_type == MT && W.data_type == MT &&
           V.data_type == MT && T.data_type == MT && WV.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        claqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, w, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
    } else {
        zlaqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, w, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
    }
}

_laqr0_laqr4_real :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int, h: *$H/ArrayType2D, $$ldh: int,
    wr: *$WR/ArrayType1D, wi: *$WI/ArrayType1D, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT && WR.data_type == RT && WI.data_type == RT && Z.data_type == RT && WORK.data_type == RT;
} {
    #if RT == float32 {
        slaqr4(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
    } else {
        dlaqr4(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
    }
}

_laqr0_laqr4_complex :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int, h: *$H/ArrayType2D, $$ldh: int,
    w: *$W/ArrayType1D, $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           H.data_type == MT && W.data_type == MT && Z.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        claqr4(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
    } else {
        zlaqr4(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
    }
}

_laqr0_laqr5_real :: (
    $$wantt: bool, $$wantz: bool, $$kacc22: int, $$n: int, $$ktop: int, $$kbot: int, $$nshfts: int,
    sr: $SR/ArrayType1D, si: $SI/ArrayType1D, h: *$H/ArrayType2D, $$ldh: int,
    $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int, v: *$V/ArrayType2D, $$ldv: int,
    u: *$U/ArrayType2D, $$ldu: int, $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int,
    $$nh: int, wh: *$WH/ArrayType2D, $$ldwh: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           SR.data_type == RT && SI.data_type == RT && H.data_type == RT && Z.data_type == RT &&
           V.data_type == RT && U.data_type == RT && WV.data_type == RT && WH.data_type == RT;
} {
    #if RT == float32 {
        slaqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
    } else {
        dlaqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
    }
}

_laqr0_laqr5_complex :: (
    $$wantt: bool, $$wantz: bool, $$kacc22: int, $$n: int, $$ktop: int, $$kbot: int, $$nshfts: int,
    s: $S/ArrayType1D, h: *$H/ArrayType2D, $$ldh: int,
    $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int, v: *$V/ArrayType2D, $$ldv: int,
    u: *$U/ArrayType2D, $$ldu: int, $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int,
    $$nh: int, wh: *$WH/ArrayType2D, $$ldwh: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           S.data_type == MT && H.data_type == MT && Z.data_type == MT &&
           V.data_type == MT && U.data_type == MT && WV.data_type == MT && WH.data_type == MT;
} {
    #if MT == C32 {
        claqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
    } else {
        zlaqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
    }
}

_laqr0_lanv2 :: (
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

_laqr0_ilaenv_job :: inline ($$wantt: bool, $$wantz: bool) -> string {
    if wantt {
        if wantz return "SV";
        return "SN";
    }
    if wantz return "EV";
    return "EN";
}

_laqr0_real :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int,
    h: *$H/ArrayType2D, $$ldh: int, wr: *$WR/ArrayType1D, wi: *$WI/ArrayType1D,
    $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT &&
           WR.data_type == RT &&
           WI.data_type == RT &&
           Z.data_type == RT &&
           WORK.data_type == RT;
} {
    NTINY :: 15;
    KEXNW :: 5;
    KEXSH :: 6;
    WILK1 :: #run ncast(RT, 0.75);
    WILK2 :: #run ncast(RT, -0.4375);
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    routine_name :: #run if RT == float32 { "SLAQR0" } else { "DLAQR0" };

    zdum_storage: [1] RT;
    zdum_view := array2d_view(zdum_storage[0..], 1, row_major=false);

    aa, bb, cc, cs, dd, rt1i, rt1r, rt2i, rt2r, sn, ss, swap : RT = ---;
    i, inf, it, itmax, k, kacc22, kbot, kdu, ks, kt, ktop, ku, kv, kwh, kwtop, kwv,
    ld, ls, lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns, nsmax, nsr, nve, nw, nwmax,
    nwr, nwupbd : int = ---;
    sorted : bool = ---;
    jbcmpz := _laqr0_ilaenv_job(wantt, wantz);

    info.* = 0;

    if n == 0 {
        _laqr0_set_work1_real(work, one_r, RT);
        return;
    }

    if n <= NTINY {
        lwkopt = 1;
        if lwork != -1 {
            _laqr0_lahqr_real(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info, RT);
        }
    } else {
        info.* = 0;

        nwr = ilaenv(13, routine_name, jbcmpz, n, ilo, ihi, lwork);
        nwr = max(2, nwr);
        nwr = min(min(ihi - ilo + 1, (n - 1) / 3), nwr);

        nsr = ilaenv(15, routine_name, jbcmpz, n, ilo, ihi, lwork);
        nsr = min(min(nsr, (n - 3) / 6), ihi - ilo);
        nsr = max(2, nsr - mod(nsr, 2));

        _laqr0_laqr3_real(wantt, wantz, n, ilo, ihi, nwr + 1, h, ldh, iloz, ihiz, z, ldz,
                          ls, ld, wr, wi, h, ldh, n, h, ldh, n, h, ldh, work, -1, RT);
        lwkopt = max(3*nsr / 2, _laqr0_query_to_int(fget(work.*, 1), RT, RT));

        if lwork == -1 {
            _laqr0_set_work1_real(work, _laqr0_roundup_lwork(lwkopt, RT), RT);
            return;
        }

        nmin = ilaenv(12, routine_name, jbcmpz, n, ilo, ihi, lwork);
        nmin = max(NTINY, nmin);

        nibble = ilaenv(14, routine_name, jbcmpz, n, ilo, ihi, lwork);
        nibble = max(0, nibble);

        kacc22 = ilaenv(16, routine_name, jbcmpz, n, ilo, ihi, lwork);
        kacc22 = max(0, min(2, kacc22));

        nwmax = min((n - 1) / 3, lwork / 2);
        nw = nwmax;

        nsmax = min((n - 3) / 6, 2*lwork / 3);
        nsmax = nsmax - mod(nsmax, 2);

        ndfl = 1;
        itmax = max(30, 2*KEXSH)*max(10, ihi - ilo + 1);
        kbot = ihi;

        it = 1;
        while it <= itmax {
            if kbot < ilo break;

            k = ilo;
            kk := kbot;
            while kk >= ilo + 1 {
                if fget(h.*, kk, kk - 1) == zero_r {
                    k = kk;
                    break;
                }
                kk -= 1;
            }
            ktop = k;

            nh = kbot - ktop + 1;
            nwupbd = min(nh, nwmax);
            if ndfl < KEXNW {
                nw = min(nwupbd, nwr);
            } else {
                nw = min(nwupbd, 2*nw);
            }
            if nw < nwmax {
                if nw >= nh - 1 {
                    nw = nh;
                } else {
                    kwtop = kbot - nw + 1;
                    if abs(fget(h.*, kwtop, kwtop - 1)) > abs(fget(h.*, kwtop - 1, kwtop - 2)) {
                        nw += 1;
                    }
                }
            }
            if ndfl < KEXNW {
                ndec = -1;
            } else if ndec >= 0 || nw >= nwupbd {
                ndec += 1;
                if nw - ndec < 2 {
                    ndec = 0;
                }
                nw -= ndec;
            }

            kv = n - nw + 1;
            kt = nw + 1;
            nho = (n - nw - 1) - kt + 1;
            kwv = nw + 2;
            nve = (n - nw) - kwv + 1;

            v_view  := array2d_view(_laqr0_matrix_storage(h.*, ldh, kv, 1, n - kv + 1, n), ldh, row_major=false);
            t_view  := array2d_view(_laqr0_matrix_storage(h.*, ldh, kv, kt, n - kv + 1, n - kt + 1), ldh, row_major=false);
            wv_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, kwv, 1, n - kwv + 1, n), ldh, row_major=false);

            _laqr0_laqr3_real(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz,
                              ls, ld, wr, wi, *v_view, ldh, nho, *t_view, ldh, nve, *wv_view, ldh,
                              work, lwork, RT);

            kbot -= ld;
            ks = kbot - ls + 1;

            if ld == 0 || ((100*ld <= nw*nibble) && (kbot - ktop + 1 > min(nmin, nwmax))) {
                ns = min(min(nsmax, nsr), max(2, kbot - ktop));
                ns = ns - mod(ns, 2);

                if mod(ndfl, KEXSH) == 0 {
                    ks = kbot - ns + 1;
                    i = kbot;
                    while i >= max(ks + 1, ktop + 2) {
                        ss = abs(fget(h.*, i, i - 1)) + abs(fget(h.*, i - 1, i - 2));
                        aa = WILK1*ss + fget(h.*, i, i);
                        bb = ss;
                        cc = WILK2*ss;
                        dd = aa;
                        _laqr0_lanv2(aa, bb, cc, dd, rt1r, rt1i, rt2r, rt2i, cs, sn, RT);
                        fset(wr, i - 1, rt1r);
                        fset(wi, i - 1, rt1i);
                        fset(wr, i, rt2r);
                        fset(wi, i, rt2i);
                        i -= 2;
                    }
                    if ks == ktop {
                        fset(wr, ks + 1, fget(h.*, ks + 1, ks + 1));
                        fset(wi, ks + 1, zero_r);
                        fset(wr, ks, fget(wr.*, ks + 1));
                        fset(wi, ks, fget(wi.*, ks + 1));
                    }
                } else {
                    if kbot - ks + 1 <= ns / 2 {
                        ks = kbot - ns + 1;
                        kt = n - ns + 1;
                        h_trail_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, ks, ks, ns, ns), ldh, row_major=false);
                        h_scratch_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, kt, 1, n - kt + 1, n), ldh, row_major=false);
                        _laqr0_lacpy(#char "A", ns, ns, *h_trail_view, ldh, *h_scratch_view, ldh, RT);
                        wr_tail := array1d_view(_laqr0_subvector(wr.*, ks, wr.count - ks + 1));
                        wi_tail := array1d_view(_laqr0_subvector(wi.*, ks, wi.count - ks + 1));
                        if ns > nmin {
                            _laqr0_laqr4_real(false, false, ns, 1, ns, *h_scratch_view, ldh, *wr_tail, *wi_tail, 1, 1, *zdum_view, 1, work, lwork, inf, RT);
                        } else {
                            _laqr0_lahqr_real(false, false, ns, 1, ns, *h_scratch_view, ldh, *wr_tail, *wi_tail, 1, 1, *zdum_view, 1, inf, RT);
                        }
                        ks += inf;

                        if ks >= kbot {
                            aa = fget(h.*, kbot - 1, kbot - 1);
                            cc = fget(h.*, kbot, kbot - 1);
                            bb = fget(h.*, kbot - 1, kbot);
                            dd = fget(h.*, kbot, kbot);
                            _laqr0_lanv2(aa, bb, cc, dd, rt1r, rt1i, rt2r, rt2i, cs, sn, RT);
                            fset(wr, kbot - 1, rt1r);
                            fset(wi, kbot - 1, rt1i);
                            fset(wr, kbot, rt2r);
                            fset(wi, kbot, rt2i);
                            ks = kbot - 1;
                        }
                    }

                    if kbot - ks + 1 > ns {
                        sorted = false;
                        for <k_loop:kbot..ks + 1 {
                            k = k_loop;
                            if sorted break;
                            sorted = true;
                            for i_loop:ks..k - 1 {
                                i = i_loop;
                                if abs(fget(wr.*, i)) + abs(fget(wi.*, i)) < abs(fget(wr.*, i + 1)) + abs(fget(wi.*, i + 1)) {
                                    sorted = false;
                                    swap = fget(wr.*, i);
                                    fset(wr, i, fget(wr.*, i + 1));
                                    fset(wr, i + 1, swap);
                                    swap = fget(wi.*, i);
                                    fset(wi, i, fget(wi.*, i + 1));
                                    fset(wi, i + 1, swap);
                                }
                            }
                        }
                    }

                    i = kbot;
                    while i >= ks + 2 {
                        if fget(wi.*, i) != -fget(wi.*, i - 1) {
                            swap = fget(wr.*, i);
                            fset(wr, i, fget(wr.*, i - 1));
                            fset(wr, i - 1, fget(wr.*, i - 2));
                            fset(wr, i - 2, swap);
                            swap = fget(wi.*, i);
                            fset(wi, i, fget(wi.*, i - 1));
                            fset(wi, i - 1, fget(wi.*, i - 2));
                            fset(wi, i - 2, swap);
                        }
                        i -= 2;
                    }
                }

                if kbot - ks + 1 == 2 && fget(wi.*, kbot) == zero_r {
                    if abs(fget(wr.*, kbot) - fget(h.*, kbot, kbot)) < abs(fget(wr.*, kbot - 1) - fget(h.*, kbot, kbot)) {
                        fset(wr, kbot - 1, fget(wr.*, kbot));
                    } else {
                        fset(wr, kbot, fget(wr.*, kbot - 1));
                    }
                }

                ns = min(ns, kbot - ks + 1);
                ns = ns - mod(ns, 2);
                ks = kbot - ns + 1;

                kdu = 2*ns;
                ku = n - kdu + 1;
                kwh = kdu + 1;
                nho = (n - kdu + 1 - 4) - (kdu + 1) + 1;
                kwv = kdu + 4;
                nve = n - kdu - kwv + 1;

                sr_view := array1d_view(_laqr0_subvector(wr.*, ks, wr.count - ks + 1));
                si_view := array1d_view(_laqr0_subvector(wi.*, ks, wi.count - ks + 1));
                u_view  := array2d_view(_laqr0_matrix_storage(h.*, ldh, ku, 1, n - ku + 1, n), ldh, row_major=false);
                wv2_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, kwv, 1, n - kwv + 1, n), ldh, row_major=false);
                wh_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, ku, kwh, n - ku + 1, n - kwh + 1), ldh, row_major=false);

                v_small_view := array2d_view(_laqr0_raw_slice(work.data, work.count), 3, row_major=false);
                _laqr0_laqr5_real(wantt, wantz, kacc22, n, ktop, kbot, ns, *sr_view, *si_view,
                                  h, ldh, iloz, ihiz, z, ldz, *v_small_view, 3, *u_view, ldh, nve, *wv2_view, ldh,
                                  nho, *wh_view, ldh, RT);
            }

            if ld > 0 {
                ndfl = 1;
            } else {
                ndfl += 1;
            }

            it += 1;
        }

        if it > itmax && kbot >= ilo {
            info.* = kbot;
        }
    }

    _laqr0_set_work1_real(work, _laqr0_roundup_lwork(lwkopt, RT), RT);
}

_laqr0_complex :: (
    $$wantt: bool, $$wantz: bool, $$n: int, $$ilo: int, $$ihi: int,
    h: *$H/ArrayType2D, $$ldh: int, w: *$W/ArrayType1D,
    $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT &&
           W.data_type == MT &&
           Z.data_type == MT &&
           WORK.data_type == MT;
} {
    NTINY :: 15;
    KEXNW :: 5;
    KEXSH :: 6;
    WILK1 :: #run ncast(RT, 0.75);
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_t  :: #run _laqr0_complex_from_real(ncast(RT, 2), RT, MT);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == C32 { "CLAQR0" } else { "ZLAQR0" };

    zdum_storage: [1] MT;
    zdum_view := array2d_view(zdum_storage[0..], 1, row_major=false);

    aa, bb, cc, dd, det, rtdisc, s, swap, tr2 : MT = ---;
    i, inf, it, itmax, k, kacc22, kbot, kdu, ks, kt, ktop, ku, kv, kwh, kwtop, kwv,
    ld, ls, lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns, nsmax, nsr, nve, nw, nwmax,
    nwr, nwupbd : int = ---;
    sorted : bool = ---;
    jbcmpz := _laqr0_ilaenv_job(wantt, wantz);

    info.* = 0;

    if n == 0 {
        _laqr0_set_work1_complex(work, one_r, RT, MT);
        return;
    }

    if n <= NTINY {
        lwkopt = 1;
        if lwork != -1 {
            _laqr0_lahqr_complex(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info, MT);
        }
    } else {
        info.* = 0;

        nwr = ilaenv(13, routine_name, jbcmpz, n, ilo, ihi, lwork);
        nwr = max(2, nwr);
        nwr = min(min(ihi - ilo + 1, (n - 1) / 3), nwr);

        nsr = ilaenv(15, routine_name, jbcmpz, n, ilo, ihi, lwork);
        nsr = min(min(nsr, (n - 3) / 6), ihi - ilo);
        nsr = max(2, nsr - mod(nsr, 2));

        _laqr0_laqr3_complex(wantt, wantz, n, ilo, ihi, nwr + 1, h, ldh, iloz, ihiz, z, ldz,
                             ls, ld, w, h, ldh, n, h, ldh, n, h, ldh, work, -1, MT);
        lwkopt = max(3*nsr / 2, _laqr0_query_to_int(fget(work.*, 1), RT, MT));

        if lwork == -1 {
            _laqr0_set_work1_complex(work, ncast(RT, lwkopt), RT, MT);
            return;
        }

        nmin = ilaenv(12, routine_name, jbcmpz, n, ilo, ihi, lwork);
        nmin = max(NTINY, nmin);

        nibble = ilaenv(14, routine_name, jbcmpz, n, ilo, ihi, lwork);
        nibble = max(0, nibble);

        kacc22 = ilaenv(16, routine_name, jbcmpz, n, ilo, ihi, lwork);
        kacc22 = max(0, min(2, kacc22));

        nwmax = min((n - 1) / 3, lwork / 2);
        nw = nwmax;

        nsmax = min((n - 3) / 6, 2*lwork / 3);
        nsmax = nsmax - mod(nsmax, 2);

        ndfl = 1;
        itmax = max(30, 2*KEXSH)*max(10, ihi - ilo + 1);
        kbot = ihi;

        it = 1;
        while it <= itmax {
            if kbot < ilo break;

            k = ilo;
            kk := kbot;
            while kk >= ilo + 1 {
                if fget(h.*, kk, kk - 1) == zero_t {
                    k = kk;
                    break;
                }
                kk -= 1;
            }
            ktop = k;

            nh = kbot - ktop + 1;
            nwupbd = min(nh, nwmax);
            if ndfl < KEXNW {
                nw = min(nwupbd, nwr);
            } else {
                nw = min(nwupbd, 2*nw);
            }
            if nw < nwmax {
                if nw >= nh - 1 {
                    nw = nh;
                } else {
                    kwtop = kbot - nw + 1;
                    if _laqr0_abs1_complex(fget(h.*, kwtop, kwtop - 1), RT, MT) >
                       _laqr0_abs1_complex(fget(h.*, kwtop - 1, kwtop - 2), RT, MT) {
                        nw += 1;
                    }
                }
            }
            if ndfl < KEXNW {
                ndec = -1;
            } else if ndec >= 0 || nw >= nwupbd {
                ndec += 1;
                if nw - ndec < 2 {
                    ndec = 0;
                }
                nw -= ndec;
            }

            kv = n - nw + 1;
            kt = nw + 1;
            nho = (n - nw - 1) - kt + 1;
            kwv = nw + 2;
            nve = (n - nw) - kwv + 1;

            v_view  := array2d_view(_laqr0_matrix_storage(h.*, ldh, kv, 1, n - kv + 1, n), ldh, row_major=false);
            t_view  := array2d_view(_laqr0_matrix_storage(h.*, ldh, kv, kt, n - kv + 1, n - kt + 1), ldh, row_major=false);
            wv_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, kwv, 1, n - kwv + 1, n), ldh, row_major=false);

            _laqr0_laqr3_complex(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz,
                                 ls, ld, w, *v_view, ldh, nho, *t_view, ldh, nve, *wv_view, ldh,
                                 work, lwork, MT);

            kbot -= ld;
            ks = kbot - ls + 1;

            if ld == 0 || ((100*ld <= nw*nibble) && (kbot - ktop + 1 > min(nmin, nwmax))) {
                ns = min(min(nsmax, nsr), max(2, kbot - ktop));
                ns = ns - mod(ns, 2);

                if mod(ndfl, KEXSH) == 0 {
                    ks = kbot - ns + 1;
                    i = kbot;
                    while i >= ks + 1 {
                        fset(w, i, fget(h.*, i, i) + _laqr0_complex_from_real(WILK1, RT, MT)*_laqr0_complex_from_real(_laqr0_abs1_complex(fget(h.*, i, i - 1), RT, MT), RT, MT));
                        fset(w, i - 1, fget(w.*, i));
                        i -= 2;
                    }
                } else {
                    if kbot - ks + 1 <= ns / 2 {
                        ks = kbot - ns + 1;
                        kt = n - ns + 1;
                        h_trail_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, ks, ks, ns, ns), ldh, row_major=false);
                        h_scratch_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, kt, 1, n - kt + 1, n), ldh, row_major=false);
                        _laqr0_lacpy(#char "A", ns, ns, *h_trail_view, ldh, *h_scratch_view, ldh, MT);
                        w_tail := array1d_view(_laqr0_subvector(w.*, ks, w.count - ks + 1));
                        if ns > nmin {
                            _laqr0_laqr4_complex(false, false, ns, 1, ns, *h_scratch_view, ldh, *w_tail, 1, 1, *zdum_view, 1, work, lwork, inf, MT);
                        } else {
                            _laqr0_lahqr_complex(false, false, ns, 1, ns, *h_scratch_view, ldh, *w_tail, 1, 1, *zdum_view, 1, inf, MT);
                        }
                        ks += inf;

                        if ks >= kbot {
                            s = _laqr0_complex_from_real(_laqr0_abs1_complex(fget(h.*, kbot - 1, kbot - 1), RT, MT) +
                                                         _laqr0_abs1_complex(fget(h.*, kbot, kbot - 1), RT, MT) +
                                                         _laqr0_abs1_complex(fget(h.*, kbot - 1, kbot), RT, MT) +
                                                         _laqr0_abs1_complex(fget(h.*, kbot, kbot), RT, MT), RT, MT);
                            aa = fget(h.*, kbot - 1, kbot - 1) / s;
                            cc = fget(h.*, kbot, kbot - 1) / s;
                            bb = fget(h.*, kbot - 1, kbot) / s;
                            dd = fget(h.*, kbot, kbot) / s;
                            tr2 = (aa + dd) / two_t;
                            det = (aa - tr2)*(dd - tr2) - bb*cc;
                            rtdisc = sqrt(-det);
                            fset(w, kbot - 1, (tr2 + rtdisc)*s);
                            fset(w, kbot, (tr2 - rtdisc)*s);
                            ks = kbot - 1;
                        }
                    }

                    if kbot - ks + 1 > ns {
                        sorted = false;
                        for <k_loop:kbot..ks + 1 {
                            k = k_loop;
                            if sorted break;
                            sorted = true;
                            for i_loop:ks..k - 1 {
                                i = i_loop;
                                if _laqr0_abs1_complex(fget(w.*, i), RT, MT) < _laqr0_abs1_complex(fget(w.*, i + 1), RT, MT) {
                                    sorted = false;
                                    swap = fget(w.*, i);
                                    fset(w, i, fget(w.*, i + 1));
                                    fset(w, i + 1, swap);
                                }
                            }
                        }
                    }
                }

                if kbot - ks + 1 == 2 {
                    if _laqr0_abs1_complex(fget(w.*, kbot) - fget(h.*, kbot, kbot), RT, MT) <
                       _laqr0_abs1_complex(fget(w.*, kbot - 1) - fget(h.*, kbot, kbot), RT, MT) {
                        fset(w, kbot - 1, fget(w.*, kbot));
                    } else {
                        fset(w, kbot, fget(w.*, kbot - 1));
                    }
                }

                ns = min(ns, kbot - ks + 1);
                ns = ns - mod(ns, 2);
                ks = kbot - ns + 1;

                kdu = 2*ns;
                ku = n - kdu + 1;
                kwh = kdu + 1;
                nho = (n - kdu + 1 - 4) - (kdu + 1) + 1;
                kwv = kdu + 4;
                nve = n - kdu - kwv + 1;

                s_view := array1d_view(_laqr0_subvector(w.*, ks, w.count - ks + 1));
                u_view  := array2d_view(_laqr0_matrix_storage(h.*, ldh, ku, 1, n - ku + 1, n), ldh, row_major=false);
                wv2_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, kwv, 1, n - kwv + 1, n), ldh, row_major=false);
                wh_view := array2d_view(_laqr0_matrix_storage(h.*, ldh, ku, kwh, n - ku + 1, n - kwh + 1), ldh, row_major=false);

                v_small_view := array2d_view(_laqr0_raw_slice(work.data, work.count), 3, row_major=false);
                _laqr0_laqr5_complex(wantt, wantz, kacc22, n, ktop, kbot, ns, *s_view,
                                     h, ldh, iloz, ihiz, z, ldz, *v_small_view, 3, *u_view, ldh, nve, *wv2_view, ldh,
                                     nho, *wh_view, ldh, MT);
            }

            if ld > 0 {
                ndfl = 1;
            } else {
                ndfl += 1;
            }

            it += 1;
        }

        if it > itmax && kbot >= ilo {
            info.* = kbot;
        }
    }

    _laqr0_set_work1_complex(work, ncast(RT, lwkopt), RT, MT);
}
