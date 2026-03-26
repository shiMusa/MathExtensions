// ? LAPACK *laqr5()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqr5/dlaqr5/claqr5/zlaqr5
// *LAQR5 carries out a blocked small-bulge multishift QR sweep.

slaqr5 :: #bake_arguments _laqr5_real(RT=float32);
dlaqr5 :: #bake_arguments _laqr5_real(RT=float64);
claqr5 :: #bake_arguments _laqr5_complex(RT=float32, MT=C32);
zlaqr5 :: #bake_arguments _laqr5_complex(RT=float64, MT=C64);

#scope_file
_laqr5_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_laqr5_storage :: inline (data: *$T, $$count: int, $$inc: int) -> [] T {
    storage_count := ifx count <= 0 then 0 else 1 + (count - 1)*inc;
    return _laqr5_raw_slice(data, storage_count);
}

_laqr5_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else lda*(cols - 1) + rows;
    return _laqr5_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqr5_abs1 :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return abs(real(x)) + abs(imag(x));
}

_laqr5_laset :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $MT, beta: MT, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) && A.data_type == MT;
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

_laqr5_lacpy :: (
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

_laqr5_gemm :: (
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

_laqr5_laqr1_real :: (
    $$n: int, h: $H/ArrayType2D, $$ldh: int, $$sr1: $RT, $$si1: RT, $$sr2: RT, $$si2: RT, v: *$V/ArrayType1D, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && H.data_type == RT && V.data_type == RT;
} {
    #if RT == float32 {
        slaqr1(n, h, ldh, sr1, si1, sr2, si2, v);
    } else {
        dlaqr1(n, h, ldh, sr1, si1, sr2, si2, v);
    }
}

_laqr5_laqr1_complex :: (
    $$n: int, h: $H/ArrayType2D, $$ldh: int, s1: $MT, s2: MT, v: *$V/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT && V.data_type == MT;
} {
    #if MT == C32 {
        claqr1(n, h, ldh, s1, s2, v);
    } else {
        zlaqr1(n, h, ldh, s1, s2, v);
    }
}

_laqr5_larfg_real :: ($$n: int, alpha: *$RT, x: *$X/ArrayType1D, $$incx: int, tau: *RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 { slarfg(n, alpha, x, incx, tau); } else { dlarfg(n, alpha, x, incx, tau); }
}

_laqr5_larfg_complex :: ($$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *MT, $MT: Type) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { clarfg(n, alpha, x, incx, tau); } else { zlarfg(n, alpha, x, incx, tau); }
}

_laqr5_real :: (
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
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);

    alpha, beta, h11, h12, h21, h22, refsum, safmax, safmin, scl, smlnum, swap,
    t1, t2, t3, tst1, tst2, ulp : RT = ---;
    i, i2, i4, incol, j, jbot, jcol, jlen, jrow, jtop, k, k1, kdu, kms, krcol,
    m, m22, mbot, mtop, nbmps, ndcol, ns, nu : int = ---;
    accum, bmp22 : bool = ---;
    tau_local : RT = ---;
    vt_storage: [3] RT;
    vt := array1d_view(vt_storage[0..]);

    if nshfts < 2 then return;
    if ktop >= kbot then return;

    i = 1;
    while i <= nshfts - 2 {
        if fget(si, i) != -fget(si, i + 1) {
            swap = fget(sr, i);
            fset(sr, i, fget(sr, i + 1));
            fset(sr, i + 1, fget(sr, i + 2));
            fset(sr, i + 2, swap);

            swap = fget(si, i);
            fset(si, i, fget(si, i + 1));
            fset(si, i + 1, fget(si, i + 2));
            fset(si, i + 2, swap);
        }
        i += 2;
    }

    ns = nshfts - mod(nshfts, 2);
    safmin = ifx RT == float32 then slamch("Safe minimum".data) else dlamch("Safe minimum".data);
    safmax = one_r / safmin;
    ulp = ifx RT == float32 then slamch("Precision".data) else dlamch("Precision".data);
    smlnum = safmin*(ncast(RT, n) / ulp);
    _ = safmax;

    accum = kacc22 == 1 || kacc22 == 2;

    if ktop + 2 <= kbot {
        fset(h, ktop + 2, ktop, zero_r);
    }

    nbmps = ns / 2;
    kdu = 4*nbmps;

    incol = ktop - 2*nbmps + 1;
    while incol <= kbot - 2 {
        if accum {
            jtop = max(ktop, incol);
        } else if wantt {
            jtop = 1;
        } else {
            jtop = ktop;
        }

        ndcol = incol + kdu;
        if accum {
            _laqr5_laset(#char "A", kdu, kdu, zero_r, one_r, u, ldu, RT);
        }

        krcol = incol;
        while krcol <= min(incol + 2*nbmps - 1, kbot - 2) {
            mtop = max(1, (ktop - krcol)/2 + 1);
            mbot = min(nbmps, (kbot - krcol - 1)/2);
            m22 = mbot + 1;
            bmp22 = (mbot < nbmps) && (krcol + 2*(m22 - 1) == (kbot - 2));

            if bmp22 {
                k = krcol + 2*(m22 - 1);
                v2col := array1d_view(_laqr5_raw_slice(v.data + (m22 - 1)*ldv, 2));
                v2tail := array1d_view(_laqr5_raw_slice(v2col.data + 1, 1));

                if k == ktop - 1 {
                    h2_view := array2d_view(_laqr5_matrix_storage(h.*, ldh, k + 1, k + 1, 2, 2), ldh, row_major=false);
                    _laqr5_laqr1_real(2, *h2_view, ldh, fget(sr, 2*m22 - 1), fget(si, 2*m22 - 1), fget(sr, 2*m22), fget(si, 2*m22), *v2col, RT);
                    beta = fget(v.*, 1, m22);
                    _laqr5_larfg_real(2, beta, *v2tail, 1, tau_local, RT);
                    fset(v, 1, m22, tau_local);
                } else {
                    beta = fget(h.*, k + 1, k);
                    fset(v, 2, m22, fget(h.*, k + 2, k));
                    _laqr5_larfg_real(2, beta, *v2tail, 1, tau_local, RT);
                    fset(v, 1, m22, tau_local);
                    fset(h, k + 1, k, beta);
                    fset(h, k + 2, k, zero_r);
                }

                t1 = fget(v.*, 1, m22);
                t2 = t1*fget(v.*, 2, m22);
                for j_loop:jtop..min(kbot, k + 3) {
                    j = j_loop;
                    refsum = fget(h.*, j, k + 1) + fget(v.*, 2, m22)*fget(h.*, j, k + 2);
                    fset(h, j, k + 1, fget(h.*, j, k + 1) - refsum*t1);
                    fset(h, j, k + 2, fget(h.*, j, k + 2) - refsum*t2);
                }

                if accum {
                    jbot = min(ndcol, kbot);
                } else if wantt {
                    jbot = n;
                } else {
                    jbot = kbot;
                }

                t1 = fget(v.*, 1, m22);
                t2 = t1*fget(v.*, 2, m22);
                for j_loop:k + 1..jbot {
                    j = j_loop;
                    refsum = fget(h.*, k + 1, j) + fget(v.*, 2, m22)*fget(h.*, k + 2, j);
                    fset(h, k + 1, j, fget(h.*, k + 1, j) - refsum*t1);
                    fset(h, k + 2, j, fget(h.*, k + 2, j) - refsum*t2);
                }

                if k >= ktop && fget(h.*, k + 1, k) != zero_r {
                    tst1 = abs(fget(h.*, k, k)) + abs(fget(h.*, k + 1, k + 1));
                    if tst1 == zero_r {
                        if k >= ktop + 1 { tst1 += abs(fget(h.*, k, k - 1)); }
                        if k >= ktop + 2 { tst1 += abs(fget(h.*, k, k - 2)); }
                        if k >= ktop + 3 { tst1 += abs(fget(h.*, k, k - 3)); }
                        if k <= kbot - 2 { tst1 += abs(fget(h.*, k + 2, k + 1)); }
                        if k <= kbot - 3 { tst1 += abs(fget(h.*, k + 3, k + 1)); }
                        if k <= kbot - 4 { tst1 += abs(fget(h.*, k + 4, k + 1)); }
                    }
                    if abs(fget(h.*, k + 1, k)) <= max(smlnum, ulp*tst1) {
                        h12 = max(abs(fget(h.*, k + 1, k)), abs(fget(h.*, k, k + 1)));
                        h21 = min(abs(fget(h.*, k + 1, k)), abs(fget(h.*, k, k + 1)));
                        h11 = max(abs(fget(h.*, k + 1, k + 1)), abs(fget(h.*, k, k) - fget(h.*, k + 1, k + 1)));
                        h22 = min(abs(fget(h.*, k + 1, k + 1)), abs(fget(h.*, k, k) - fget(h.*, k + 1, k + 1)));
                        scl = h11 + h12;
                        tst2 = h22*(h11 / scl);
                        if tst2 == zero_r || h21*(h12 / scl) <= max(smlnum, ulp*tst2) {
                            fset(h, k + 1, k, zero_r);
                        }
                    }
                }

                if accum {
                    kms = k - incol;
                    t1 = fget(v.*, 1, m22);
                    t2 = t1*fget(v.*, 2, m22);
                    for j_loop:max(1, ktop - incol)..kdu {
                        j = j_loop;
                        refsum = fget(u.*, j, kms + 1) + fget(v.*, 2, m22)*fget(u.*, j, kms + 2);
                        fset(u, j, kms + 1, fget(u.*, j, kms + 1) - refsum*t1);
                        fset(u, j, kms + 2, fget(u.*, j, kms + 2) - refsum*t2);
                    }
                } else if wantz {
                    t1 = fget(v.*, 1, m22);
                    t2 = t1*fget(v.*, 2, m22);
                    for j_loop:iloz..ihiz {
                        j = j_loop;
                        refsum = fget(z.*, j, k + 1) + fget(v.*, 2, m22)*fget(z.*, j, k + 2);
                        fset(z, j, k + 1, fget(z.*, j, k + 1) - refsum*t1);
                        fset(z, j, k + 2, fget(z.*, j, k + 2) - refsum*t2);
                    }
                }
            }

            for <m_loop:mbot..mtop {
                m = m_loop;
                k = krcol + 2*(m - 1);
                v3col := array1d_view(_laqr5_raw_slice(v.data + (m - 1)*ldv, 3));
                v3tail := array1d_view(_laqr5_raw_slice(v3col.data + 1, 2));

                if k == ktop - 1 {
                    h3_view := array2d_view(_laqr5_matrix_storage(h.*, ldh, ktop, ktop, 3, 3), ldh, row_major=false);
                    _laqr5_laqr1_real(3, *h3_view, ldh, fget(sr, 2*m - 1), fget(si, 2*m - 1), fget(sr, 2*m), fget(si, 2*m), *v3col, RT);
                    alpha = fget(v.*, 1, m);
                    _laqr5_larfg_real(3, alpha, *v3tail, 1, tau_local, RT);
                    fset(v, 1, m, tau_local);
                } else {
                    t1 = fget(v.*, 1, m);
                    t2 = t1*fget(v.*, 2, m);
                    t3 = t1*fget(v.*, 3, m);
                    refsum = fget(v.*, 3, m)*fget(h.*, k + 3, k + 2);
                    fset(h, k + 3, k, -refsum*t1);
                    fset(h, k + 3, k + 1, -refsum*t2);
                    fset(h, k + 3, k + 2, fget(h.*, k + 3, k + 2) - refsum*t3);

                    beta = fget(h.*, k + 1, k);
                    fset(v, 2, m, fget(h.*, k + 2, k));
                    fset(v, 3, m, fget(h.*, k + 3, k));
                    _laqr5_larfg_real(3, beta, *v3tail, 1, tau_local, RT);
                    fset(v, 1, m, tau_local);

                    if fget(h.*, k + 3, k) != zero_r || fget(h.*, k + 3, k + 1) != zero_r || fget(h.*, k + 3, k + 2) == zero_r {
                        fset(h, k + 1, k, beta);
                        fset(h, k + 2, k, zero_r);
                        fset(h, k + 3, k, zero_r);
                    } else {
                        h3_view := array2d_view(_laqr5_matrix_storage(h.*, ldh, k + 1, k + 1, 3, 3), ldh, row_major=false);
                        _laqr5_laqr1_real(3, *h3_view, ldh, fget(sr, 2*m - 1), fget(si, 2*m - 1), fget(sr, 2*m), fget(si, 2*m), *vt, RT);
                        alpha = fget(vt, 1);
                        vt_tail := array1d_view(vt_storage[1..]);
                        _laqr5_larfg_real(3, alpha, *vt_tail, 1, tau_local, RT);
                        fset(vt, 1, tau_local);
                        t1 = fget(vt, 1);
                        t2 = t1*fget(vt, 2);
                        t3 = t2*fget(vt, 3);
                        refsum = fget(h.*, k + 1, k) + fget(vt, 2)*fget(h.*, k + 2, k);

                        if abs(fget(h.*, k + 2, k) - refsum*t2) + abs(refsum*t3) >
                           ulp*(abs(fget(h.*, k, k)) + abs(fget(h.*, k + 1, k + 1)) + abs(fget(h.*, k + 2, k + 2))) {
                            fset(h, k + 1, k, beta);
                            fset(h, k + 2, k, zero_r);
                            fset(h, k + 3, k, zero_r);
                        } else {
                            fset(h, k + 1, k, fget(h.*, k + 1, k) - refsum*t1);
                            fset(h, k + 2, k, zero_r);
                            fset(h, k + 3, k, zero_r);
                            fset(v, 1, m, fget(vt, 1));
                            fset(v, 2, m, fget(vt, 2));
                            fset(v, 3, m, fget(vt, 3));
                        }
                    }
                }

                t1 = fget(v.*, 1, m);
                t2 = t1*fget(v.*, 2, m);
                t3 = t1*fget(v.*, 3, m);
                for j_loop:jtop..min(kbot, k + 3) {
                    j = j_loop;
                    refsum = fget(h.*, j, k + 1) + fget(v.*, 2, m)*fget(h.*, j, k + 2) + fget(v.*, 3, m)*fget(h.*, j, k + 3);
                    fset(h, j, k + 1, fget(h.*, j, k + 1) - refsum*t1);
                    fset(h, j, k + 2, fget(h.*, j, k + 2) - refsum*t2);
                    fset(h, j, k + 3, fget(h.*, j, k + 3) - refsum*t3);
                }

                refsum = fget(h.*, k + 1, k + 1) + fget(v.*, 2, m)*fget(h.*, k + 2, k + 1) + fget(v.*, 3, m)*fget(h.*, k + 3, k + 1);
                fset(h, k + 1, k + 1, fget(h.*, k + 1, k + 1) - refsum*t1);
                fset(h, k + 2, k + 1, fget(h.*, k + 2, k + 1) - refsum*t2);
                fset(h, k + 3, k + 1, fget(h.*, k + 3, k + 1) - refsum*t3);

                if k < ktop {
                    continue;
                }

                if fget(h.*, k + 1, k) != zero_r {
                    tst1 = abs(fget(h.*, k, k)) + abs(fget(h.*, k + 1, k + 1));
                    if tst1 == zero_r {
                        if k >= ktop + 1 { tst1 += abs(fget(h.*, k, k - 1)); }
                        if k >= ktop + 2 { tst1 += abs(fget(h.*, k, k - 2)); }
                        if k >= ktop + 3 { tst1 += abs(fget(h.*, k, k - 3)); }
                        if k <= kbot - 2 { tst1 += abs(fget(h.*, k + 2, k + 1)); }
                        if k <= kbot - 3 { tst1 += abs(fget(h.*, k + 3, k + 1)); }
                        if k <= kbot - 4 { tst1 += abs(fget(h.*, k + 4, k + 1)); }
                    }
                    if abs(fget(h.*, k + 1, k)) <= max(smlnum, ulp*tst1) {
                        h12 = max(abs(fget(h.*, k + 1, k)), abs(fget(h.*, k, k + 1)));
                        h21 = min(abs(fget(h.*, k + 1, k)), abs(fget(h.*, k, k + 1)));
                        h11 = max(abs(fget(h.*, k + 1, k + 1)), abs(fget(h.*, k, k) - fget(h.*, k + 1, k + 1)));
                        h22 = min(abs(fget(h.*, k + 1, k + 1)), abs(fget(h.*, k, k) - fget(h.*, k + 1, k + 1)));
                        scl = h11 + h12;
                        tst2 = h22*(h11 / scl);
                        if tst2 == zero_r || h21*(h12 / scl) <= max(smlnum, ulp*tst2) {
                            fset(h, k + 1, k, zero_r);
                        }
                    }
                }
            }

            if accum {
                jbot = min(ndcol, kbot);
            } else if wantt {
                jbot = n;
            } else {
                jbot = kbot;
            }

            for <m_loop:mbot..mtop {
                m = m_loop;
                k = krcol + 2*(m - 1);
                t1 = fget(v.*, 1, m);
                t2 = t1*fget(v.*, 2, m);
                t3 = t1*fget(v.*, 3, m);
                for j_loop:max(ktop, krcol + 2*m)..jbot {
                    j = j_loop;
                    refsum = fget(h.*, k + 1, j) + fget(v.*, 2, m)*fget(h.*, k + 2, j) + fget(v.*, 3, m)*fget(h.*, k + 3, j);
                    fset(h, k + 1, j, fget(h.*, k + 1, j) - refsum*t1);
                    fset(h, k + 2, j, fget(h.*, k + 2, j) - refsum*t2);
                    fset(h, k + 3, j, fget(h.*, k + 3, j) - refsum*t3);
                }
            }

            if accum {
                for <m_loop:mbot..mtop {
                    m = m_loop;
                    k = krcol + 2*(m - 1);
                    kms = k - incol;
                    i2 = max(1, ktop - incol);
                    i2 = max(i2, kms - (krcol - incol) + 1);
                    i4 = min(kdu, krcol + 2*(mbot - 1) - incol + 5);
                    t1 = fget(v.*, 1, m);
                    t2 = t1*fget(v.*, 2, m);
                    t3 = t1*fget(v.*, 3, m);
                    for j_loop:i2..i4 {
                        j = j_loop;
                        refsum = fget(u.*, j, kms + 1) + fget(v.*, 2, m)*fget(u.*, j, kms + 2) + fget(v.*, 3, m)*fget(u.*, j, kms + 3);
                        fset(u, j, kms + 1, fget(u.*, j, kms + 1) - refsum*t1);
                        fset(u, j, kms + 2, fget(u.*, j, kms + 2) - refsum*t2);
                        fset(u, j, kms + 3, fget(u.*, j, kms + 3) - refsum*t3);
                    }
                }
            } else if wantz {
                for <m_loop:mbot..mtop {
                    m = m_loop;
                    k = krcol + 2*(m - 1);
                    t1 = fget(v.*, 1, m);
                    t2 = t1*fget(v.*, 2, m);
                    t3 = t1*fget(v.*, 3, m);
                    for j_loop:iloz..ihiz {
                        j = j_loop;
                        refsum = fget(z.*, j, k + 1) + fget(v.*, 2, m)*fget(z.*, j, k + 2) + fget(v.*, 3, m)*fget(z.*, j, k + 3);
                        fset(z, j, k + 1, fget(z.*, j, k + 1) - refsum*t1);
                        fset(z, j, k + 2, fget(z.*, j, k + 2) - refsum*t2);
                        fset(z, j, k + 3, fget(z.*, j, k + 3) - refsum*t3);
                    }
                }
            }

            krcol += 1;
        }

        if accum {
            if wantt {
                jtop = 1;
                jbot = n;
            } else {
                jtop = ktop;
                jbot = kbot;
            }
            k1 = max(1, ktop - incol);
            nu = (kdu - max(0, ndcol - kbot)) - k1 + 1;

            jcol = min(ndcol, kbot) + 1;
            while jcol <= jbot {
                jlen = min(nh, jbot - jcol + 1);
                ublk := array2d_view(_laqr5_matrix_storage(u.*, ldu, k1, k1, nu, nu), ldu, row_major=false);
                hblk := array2d_view(_laqr5_matrix_storage(h.*, ldh, incol + k1, jcol, nu, jlen), ldh, row_major=false);
                whblk := array2d_view(_laqr5_matrix_storage(wh.*, ldwh, 1, 1, nu, jlen), ldwh, row_major=false);
                _laqr5_gemm(#char "C", #char "N", nu, jlen, nu, one_r, *ublk, ldu, *hblk, ldh, zero_r, *whblk, ldwh, RT);
                _laqr5_lacpy(#char "A", nu, jlen, *whblk, ldwh, *hblk, ldh, RT);
                jcol += nh;
            }

            jrow = jtop;
            while jrow <= max(ktop, incol) - 1 {
                jlen = min(nv, max(ktop, incol) - jrow);
                hblk := array2d_view(_laqr5_matrix_storage(h.*, ldh, jrow, incol + k1, jlen, nu), ldh, row_major=false);
                ublk := array2d_view(_laqr5_matrix_storage(u.*, ldu, k1, k1, nu, nu), ldu, row_major=false);
                wvblk := array2d_view(_laqr5_matrix_storage(wv.*, ldwv, 1, 1, jlen, nu), ldwv, row_major=false);
                _laqr5_gemm(#char "N", #char "N", jlen, nu, nu, one_r, *hblk, ldh, *ublk, ldu, zero_r, *wvblk, ldwv, RT);
                _laqr5_lacpy(#char "A", jlen, nu, *wvblk, ldwv, *hblk, ldh, RT);
                jrow += nv;
            }

            if wantz {
                jrow = iloz;
                while jrow <= ihiz {
                    jlen = min(nv, ihiz - jrow + 1);
                    zblk := array2d_view(_laqr5_matrix_storage(z.*, ldz, jrow, incol + k1, jlen, nu), ldz, row_major=false);
                    ublk := array2d_view(_laqr5_matrix_storage(u.*, ldu, k1, k1, nu, nu), ldu, row_major=false);
                    wvblk := array2d_view(_laqr5_matrix_storage(wv.*, ldwv, 1, 1, jlen, nu), ldwv, row_major=false);
                    _laqr5_gemm(#char "N", #char "N", jlen, nu, nu, one_r, *zblk, ldz, *ublk, ldu, zero_r, *wvblk, ldwv, RT);
                    _laqr5_lacpy(#char "A", jlen, nu, *wvblk, ldwv, *zblk, ldz, RT);
                    jrow += nv;
                }
            }
        }

        incol += 2*nbmps;
    }
}

_laqr5_complex :: (
    $$wantt: bool, $$wantz: bool, $$kacc22: int, $$n: int, $$ktop: int, $$kbot: int, $$nshfts: int,
    s: $S/ArrayType1D, h: *$H/ArrayType2D, $$ldh: int,
    $$iloz: int, $$ihiz: int, z: *$Z/ArrayType2D, $$ldz: int, v: *$V/ArrayType2D, $$ldv: int,
    u: *$U/ArrayType2D, $$ldu: int, $$nv: int, wv: *$WV/ArrayType2D, $$ldwv: int,
    $$nh: int, wh: *$WH/ArrayType2D, $$ldwh: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           S.data_type == MT && H.data_type == MT && Z.data_type == MT &&
           V.data_type == MT && U.data_type == MT && WV.data_type == MT && WH.data_type == MT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    zero_t :: #run zero(MT);
    one_t  :: #run cpx(one_r, zero_r);

    alpha, beta, refsum, t1, t2, t3 : MT = ---;
    h11, h12, h21, h22, safmax, safmin, scl, smlnum, tst1, tst2, ulp : RT = ---;
    i2, i4, incol, j, jbot, jcol, jlen, jrow, jtop, k, k1, kdu, kms, krcol,
    m, m22, mbot, mtop, nbmps, ndcol, ns, nu : int = ---;
    accum, bmp22 : bool = ---;
    tau_local : MT = ---;
    vt_storage: [3] MT;
    vt := array1d_view(vt_storage[0..]);

    if nshfts < 2 then return;
    if ktop >= kbot then return;

    ns = nshfts - mod(nshfts, 2);
    safmin = ifx RT == float32 then slamch("Safe minimum".data) else dlamch("Safe minimum".data);
    safmax = one_r / safmin;
    ulp = ifx RT == float32 then slamch("Precision".data) else dlamch("Precision".data);
    smlnum = safmin*(ncast(RT, n) / ulp);
    _ = safmax;

    accum = kacc22 == 1 || kacc22 == 2;

    if ktop + 2 <= kbot {
        fset(h, ktop + 2, ktop, zero_t);
    }

    nbmps = ns / 2;
    kdu = 4*nbmps;

    incol = ktop - 2*nbmps + 1;
    while incol <= kbot - 2 {
        if accum {
            jtop = max(ktop, incol);
        } else if wantt {
            jtop = 1;
        } else {
            jtop = ktop;
        }

        ndcol = incol + kdu;
        if accum {
            _laqr5_laset(#char "A", kdu, kdu, zero_t, one_t, u, ldu, MT);
        }

        krcol = incol;
        while krcol <= min(incol + 2*nbmps - 1, kbot - 2) {
            mtop = max(1, (ktop - krcol)/2 + 1);
            mbot = min(nbmps, (kbot - krcol - 1)/2);
            m22 = mbot + 1;
            bmp22 = (mbot < nbmps) && (krcol + 2*(m22 - 1) == (kbot - 2));

            if bmp22 {
                k = krcol + 2*(m22 - 1);
                v2col := array1d_view(_laqr5_raw_slice(v.data + (m22 - 1)*ldv, 2));
                v2tail := array1d_view(_laqr5_raw_slice(v2col.data + 1, 1));

                if k == ktop - 1 {
                    h2_view := array2d_view(_laqr5_matrix_storage(h.*, ldh, k + 1, k + 1, 2, 2), ldh, row_major=false);
                    _laqr5_laqr1_complex(2, *h2_view, ldh, fget(s, 2*m22 - 1), fget(s, 2*m22), *v2col, RT, MT);
                    beta = fget(v.*, 1, m22);
                    _laqr5_larfg_complex(2, beta, *v2tail, 1, tau_local, MT);
                    fset(v, 1, m22, tau_local);
                } else {
                    beta = fget(h.*, k + 1, k);
                    fset(v, 2, m22, fget(h.*, k + 2, k));
                    _laqr5_larfg_complex(2, beta, *v2tail, 1, tau_local, MT);
                    fset(v, 1, m22, tau_local);
                    fset(h, k + 1, k, beta);
                    fset(h, k + 2, k, zero_t);
                }

                t1 = fget(v.*, 1, m22);
                t2 = t1*conjugate(fget(v.*, 2, m22));
                for j_loop:jtop..min(kbot, k + 3) {
                    j = j_loop;
                    refsum = fget(h.*, j, k + 1) + fget(v.*, 2, m22)*fget(h.*, j, k + 2);
                    fset(h, j, k + 1, fget(h.*, j, k + 1) - refsum*t1);
                    fset(h, j, k + 2, fget(h.*, j, k + 2) - refsum*t2);
                }

                if accum {
                    jbot = min(ndcol, kbot);
                } else if wantt {
                    jbot = n;
                } else {
                    jbot = kbot;
                }

                t1 = conjugate(fget(v.*, 1, m22));
                t2 = t1*fget(v.*, 2, m22);
                for j_loop:k + 1..jbot {
                    j = j_loop;
                    refsum = fget(h.*, k + 1, j) + conjugate(fget(v.*, 2, m22))*fget(h.*, k + 2, j);
                    fset(h, k + 1, j, fget(h.*, k + 1, j) - refsum*t1);
                    fset(h, k + 2, j, fget(h.*, k + 2, j) - refsum*t2);
                }

                if k >= ktop && fget(h.*, k + 1, k) != zero_t {
                    tst1 = _laqr5_abs1(fget(h.*, k, k), RT, MT) + _laqr5_abs1(fget(h.*, k + 1, k + 1), RT, MT);
                    if tst1 == zero_r {
                        if k >= ktop + 1 { tst1 += _laqr5_abs1(fget(h.*, k, k - 1), RT, MT); }
                        if k >= ktop + 2 { tst1 += _laqr5_abs1(fget(h.*, k, k - 2), RT, MT); }
                        if k >= ktop + 3 { tst1 += _laqr5_abs1(fget(h.*, k, k - 3), RT, MT); }
                        if k <= kbot - 2 { tst1 += _laqr5_abs1(fget(h.*, k + 2, k + 1), RT, MT); }
                        if k <= kbot - 3 { tst1 += _laqr5_abs1(fget(h.*, k + 3, k + 1), RT, MT); }
                        if k <= kbot - 4 { tst1 += _laqr5_abs1(fget(h.*, k + 4, k + 1), RT, MT); }
                    }
                    if _laqr5_abs1(fget(h.*, k + 1, k), RT, MT) <= max(smlnum, ulp*tst1) {
                        h12 = max(_laqr5_abs1(fget(h.*, k + 1, k), RT, MT), _laqr5_abs1(fget(h.*, k, k + 1), RT, MT));
                        h21 = min(_laqr5_abs1(fget(h.*, k + 1, k), RT, MT), _laqr5_abs1(fget(h.*, k, k + 1), RT, MT));
                        h11 = max(_laqr5_abs1(fget(h.*, k + 1, k + 1), RT, MT), _laqr5_abs1(fget(h.*, k, k) - fget(h.*, k + 1, k + 1), RT, MT));
                        h22 = min(_laqr5_abs1(fget(h.*, k + 1, k + 1), RT, MT), _laqr5_abs1(fget(h.*, k, k) - fget(h.*, k + 1, k + 1), RT, MT));
                        scl = h11 + h12;
                        tst2 = h22*(h11 / scl);
                        if tst2 == zero_r || h21*(h12 / scl) <= max(smlnum, ulp*tst2) {
                            fset(h, k + 1, k, zero_t);
                        }
                    }
                }

                if accum {
                    kms = k - incol;
                    for j_loop:max(1, ktop - incol)..kdu {
                        j = j_loop;
                        refsum = fget(v.*, 1, m22)*(fget(u.*, j, kms + 1) + fget(v.*, 2, m22)*fget(u.*, j, kms + 2));
                        fset(u, j, kms + 1, fget(u.*, j, kms + 1) - refsum);
                        fset(u, j, kms + 2, fget(u.*, j, kms + 2) - refsum*conjugate(fget(v.*, 2, m22)));
                    }
                } else if wantz {
                    for j_loop:iloz..ihiz {
                        j = j_loop;
                        refsum = fget(v.*, 1, m22)*(fget(z.*, j, k + 1) + fget(v.*, 2, m22)*fget(z.*, j, k + 2));
                        fset(z, j, k + 1, fget(z.*, j, k + 1) - refsum);
                        fset(z, j, k + 2, fget(z.*, j, k + 2) - refsum*conjugate(fget(v.*, 2, m22)));
                    }
                }
            }

            for <m_loop:mbot..mtop {
                m = m_loop;
                k = krcol + 2*(m - 1);
                v3col := array1d_view(_laqr5_raw_slice(v.data + (m - 1)*ldv, 3));
                v3tail := array1d_view(_laqr5_raw_slice(v3col.data + 1, 2));

                if k == ktop - 1 {
                    h3_view := array2d_view(_laqr5_matrix_storage(h.*, ldh, ktop, ktop, 3, 3), ldh, row_major=false);
                    _laqr5_laqr1_complex(3, *h3_view, ldh, fget(s, 2*m - 1), fget(s, 2*m), *v3col, RT, MT);
                    alpha = fget(v.*, 1, m);
                    _laqr5_larfg_complex(3, alpha, *v3tail, 1, tau_local, MT);
                    fset(v, 1, m, tau_local);
                } else {
                    t1 = fget(v.*, 1, m);
                    t2 = t1*conjugate(fget(v.*, 2, m));
                    t3 = t1*conjugate(fget(v.*, 3, m));
                    refsum = fget(v.*, 3, m)*fget(h.*, k + 3, k + 2);
                    fset(h, k + 3, k, -refsum*t1);
                    fset(h, k + 3, k + 1, -refsum*t2);
                    fset(h, k + 3, k + 2, fget(h.*, k + 3, k + 2) - refsum*t3);

                    beta = fget(h.*, k + 1, k);
                    fset(v, 2, m, fget(h.*, k + 2, k));
                    fset(v, 3, m, fget(h.*, k + 3, k));
                    _laqr5_larfg_complex(3, beta, *v3tail, 1, tau_local, MT);
                    fset(v, 1, m, tau_local);

                    if fget(h.*, k + 3, k) != zero_t || fget(h.*, k + 3, k + 1) != zero_t || fget(h.*, k + 3, k + 2) == zero_t {
                        fset(h, k + 1, k, beta);
                        fset(h, k + 2, k, zero_t);
                        fset(h, k + 3, k, zero_t);
                    } else {
                        h3_view := array2d_view(_laqr5_matrix_storage(h.*, ldh, k + 1, k + 1, 3, 3), ldh, row_major=false);
                        _laqr5_laqr1_complex(3, *h3_view, ldh, fget(s, 2*m - 1), fget(s, 2*m), *vt, RT, MT);
                        alpha = fget(vt, 1);
                        vt_tail := array1d_view(vt_storage[1..]);
                        _laqr5_larfg_complex(3, alpha, *vt_tail, 1, tau_local, MT);
                        fset(vt, 1, tau_local);
                        t1 = conjugate(fget(vt, 1));
                        t2 = t1*fget(vt, 2);
                        t3 = t1*fget(vt, 3);
                        refsum = fget(h.*, k + 1, k) + conjugate(fget(vt, 2))*fget(h.*, k + 2, k);

                        if _laqr5_abs1(fget(h.*, k + 2, k) - refsum*t2, RT, MT) + _laqr5_abs1(refsum*t3, RT, MT) >
                           ulp*(_laqr5_abs1(fget(h.*, k, k), RT, MT) + _laqr5_abs1(fget(h.*, k + 1, k + 1), RT, MT) + _laqr5_abs1(fget(h.*, k + 2, k + 2), RT, MT)) {
                            fset(h, k + 1, k, beta);
                            fset(h, k + 2, k, zero_t);
                            fset(h, k + 3, k, zero_t);
                        } else {
                            fset(h, k + 1, k, fget(h.*, k + 1, k) - refsum*t1);
                            fset(h, k + 2, k, zero_t);
                            fset(h, k + 3, k, zero_t);
                            fset(v, 1, m, fget(vt, 1));
                            fset(v, 2, m, fget(vt, 2));
                            fset(v, 3, m, fget(vt, 3));
                        }
                    }
                }

                t1 = fget(v.*, 1, m);
                t2 = t1*conjugate(fget(v.*, 2, m));
                t3 = t1*conjugate(fget(v.*, 3, m));
                for j_loop:jtop..min(kbot, k + 3) {
                    j = j_loop;
                    refsum = fget(h.*, j, k + 1) + fget(v.*, 2, m)*fget(h.*, j, k + 2) + fget(v.*, 3, m)*fget(h.*, j, k + 3);
                    fset(h, j, k + 1, fget(h.*, j, k + 1) - refsum*t1);
                    fset(h, j, k + 2, fget(h.*, j, k + 2) - refsum*t2);
                    fset(h, j, k + 3, fget(h.*, j, k + 3) - refsum*t3);
                }

                t1 = conjugate(fget(v.*, 1, m));
                t2 = t1*fget(v.*, 2, m);
                t3 = t1*fget(v.*, 3, m);
                refsum = fget(h.*, k + 1, k + 1) + conjugate(fget(v.*, 2, m))*fget(h.*, k + 2, k + 1) + conjugate(fget(v.*, 3, m))*fget(h.*, k + 3, k + 1);
                fset(h, k + 1, k + 1, fget(h.*, k + 1, k + 1) - refsum*t1);
                fset(h, k + 2, k + 1, fget(h.*, k + 2, k + 1) - refsum*t2);
                fset(h, k + 3, k + 1, fget(h.*, k + 3, k + 1) - refsum*t3);

                if k < ktop {
                    continue;
                }

                if fget(h.*, k + 1, k) != zero_t {
                    tst1 = _laqr5_abs1(fget(h.*, k, k), RT, MT) + _laqr5_abs1(fget(h.*, k + 1, k + 1), RT, MT);
                    if tst1 == zero_r {
                        if k >= ktop + 1 { tst1 += _laqr5_abs1(fget(h.*, k, k - 1), RT, MT); }
                        if k >= ktop + 2 { tst1 += _laqr5_abs1(fget(h.*, k, k - 2), RT, MT); }
                        if k >= ktop + 3 { tst1 += _laqr5_abs1(fget(h.*, k, k - 3), RT, MT); }
                        if k <= kbot - 2 { tst1 += _laqr5_abs1(fget(h.*, k + 2, k + 1), RT, MT); }
                        if k <= kbot - 3 { tst1 += _laqr5_abs1(fget(h.*, k + 3, k + 1), RT, MT); }
                        if k <= kbot - 4 { tst1 += _laqr5_abs1(fget(h.*, k + 4, k + 1), RT, MT); }
                    }
                    if _laqr5_abs1(fget(h.*, k + 1, k), RT, MT) <= max(smlnum, ulp*tst1) {
                        h12 = max(_laqr5_abs1(fget(h.*, k + 1, k), RT, MT), _laqr5_abs1(fget(h.*, k, k + 1), RT, MT));
                        h21 = min(_laqr5_abs1(fget(h.*, k + 1, k), RT, MT), _laqr5_abs1(fget(h.*, k, k + 1), RT, MT));
                        h11 = max(_laqr5_abs1(fget(h.*, k + 1, k + 1), RT, MT), _laqr5_abs1(fget(h.*, k, k) - fget(h.*, k + 1, k + 1), RT, MT));
                        h22 = min(_laqr5_abs1(fget(h.*, k + 1, k + 1), RT, MT), _laqr5_abs1(fget(h.*, k, k) - fget(h.*, k + 1, k + 1), RT, MT));
                        scl = h11 + h12;
                        tst2 = h22*(h11 / scl);
                        if tst2 == zero_r || h21*(h12 / scl) <= max(smlnum, ulp*tst2) {
                            fset(h, k + 1, k, zero_t);
                        }
                    }
                }
            }

            if accum {
                jbot = min(ndcol, kbot);
            } else if wantt {
                jbot = n;
            } else {
                jbot = kbot;
            }

            for <m_loop:mbot..mtop {
                m = m_loop;
                k = krcol + 2*(m - 1);
                t1 = conjugate(fget(v.*, 1, m));
                t2 = t1*fget(v.*, 2, m);
                t3 = t1*fget(v.*, 3, m);
                for j_loop:max(ktop, krcol + 2*m)..jbot {
                    j = j_loop;
                    refsum = fget(h.*, k + 1, j) + conjugate(fget(v.*, 2, m))*fget(h.*, k + 2, j) + conjugate(fget(v.*, 3, m))*fget(h.*, k + 3, j);
                    fset(h, k + 1, j, fget(h.*, k + 1, j) - refsum*t1);
                    fset(h, k + 2, j, fget(h.*, k + 2, j) - refsum*t2);
                    fset(h, k + 3, j, fget(h.*, k + 3, j) - refsum*t3);
                }
            }

            if accum {
                for <m_loop:mbot..mtop {
                    m = m_loop;
                    k = krcol + 2*(m - 1);
                    kms = k - incol;
                    i2 = max(1, ktop - incol);
                    i2 = max(i2, kms - (krcol - incol) + 1);
                    i4 = min(kdu, krcol + 2*(mbot - 1) - incol + 5);
                    t1 = fget(v.*, 1, m);
                    t2 = t1*conjugate(fget(v.*, 2, m));
                    t3 = t1*conjugate(fget(v.*, 3, m));
                    for j_loop:i2..i4 {
                        j = j_loop;
                        refsum = fget(u.*, j, kms + 1) + fget(v.*, 2, m)*fget(u.*, j, kms + 2) + fget(v.*, 3, m)*fget(u.*, j, kms + 3);
                        fset(u, j, kms + 1, fget(u.*, j, kms + 1) - refsum*t1);
                        fset(u, j, kms + 2, fget(u.*, j, kms + 2) - refsum*t2);
                        fset(u, j, kms + 3, fget(u.*, j, kms + 3) - refsum*t3);
                    }
                }
            } else if wantz {
                for <m_loop:mbot..mtop {
                    m = m_loop;
                    k = krcol + 2*(m - 1);
                    t1 = fget(v.*, 1, m);
                    t2 = t1*conjugate(fget(v.*, 2, m));
                    t3 = t1*conjugate(fget(v.*, 3, m));
                    for j_loop:iloz..ihiz {
                        j = j_loop;
                        refsum = fget(z.*, j, k + 1) + fget(v.*, 2, m)*fget(z.*, j, k + 2) + fget(v.*, 3, m)*fget(z.*, j, k + 3);
                        fset(z, j, k + 1, fget(z.*, j, k + 1) - refsum*t1);
                        fset(z, j, k + 2, fget(z.*, j, k + 2) - refsum*t2);
                        fset(z, j, k + 3, fget(z.*, j, k + 3) - refsum*t3);
                    }
                }
            }

            krcol += 1;
        }

        if accum {
            if wantt {
                jtop = 1;
                jbot = n;
            } else {
                jtop = ktop;
                jbot = kbot;
            }
            k1 = max(1, ktop - incol);
            nu = (kdu - max(0, ndcol - kbot)) - k1 + 1;

            jcol = min(ndcol, kbot) + 1;
            while jcol <= jbot {
                jlen = min(nh, jbot - jcol + 1);
                ublk := array2d_view(_laqr5_matrix_storage(u.*, ldu, k1, k1, nu, nu), ldu, row_major=false);
                hblk := array2d_view(_laqr5_matrix_storage(h.*, ldh, incol + k1, jcol, nu, jlen), ldh, row_major=false);
                whblk := array2d_view(_laqr5_matrix_storage(wh.*, ldwh, 1, 1, nu, jlen), ldwh, row_major=false);
                _laqr5_gemm(#char "C", #char "N", nu, jlen, nu, one_t, *ublk, ldu, *hblk, ldh, zero_t, *whblk, ldwh, MT);
                _laqr5_lacpy(#char "A", nu, jlen, *whblk, ldwh, *hblk, ldh, MT);
                jcol += nh;
            }

            jrow = jtop;
            while jrow <= max(ktop, incol) - 1 {
                jlen = min(nv, max(ktop, incol) - jrow);
                hblk := array2d_view(_laqr5_matrix_storage(h.*, ldh, jrow, incol + k1, jlen, nu), ldh, row_major=false);
                ublk := array2d_view(_laqr5_matrix_storage(u.*, ldu, k1, k1, nu, nu), ldu, row_major=false);
                wvblk := array2d_view(_laqr5_matrix_storage(wv.*, ldwv, 1, 1, jlen, nu), ldwv, row_major=false);
                _laqr5_gemm(#char "N", #char "N", jlen, nu, nu, one_t, *hblk, ldh, *ublk, ldu, zero_t, *wvblk, ldwv, MT);
                _laqr5_lacpy(#char "A", jlen, nu, *wvblk, ldwv, *hblk, ldh, MT);
                jrow += nv;
            }

            if wantz {
                jrow = iloz;
                while jrow <= ihiz {
                    jlen = min(nv, ihiz - jrow + 1);
                    zblk := array2d_view(_laqr5_matrix_storage(z.*, ldz, jrow, incol + k1, jlen, nu), ldz, row_major=false);
                    ublk := array2d_view(_laqr5_matrix_storage(u.*, ldu, k1, k1, nu, nu), ldu, row_major=false);
                    wvblk := array2d_view(_laqr5_matrix_storage(wv.*, ldwv, 1, 1, jlen, nu), ldwv, row_major=false);
                    _laqr5_gemm(#char "N", #char "N", jlen, nu, nu, one_t, *zblk, ldz, *ublk, ldu, zero_t, *wvblk, ldwv, MT);
                    _laqr5_lacpy(#char "A", jlen, nu, *wvblk, ldwv, *zblk, ldz, MT);
                    jrow += nv;
                }
            }
        }

        incol += 2*nbmps;
    }
}
