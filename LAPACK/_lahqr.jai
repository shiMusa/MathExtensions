// ? LAPACK *lahqr()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slahqr/dlahqr/clahqr/zlahqr
// *LAHQR computes the eigenvalues and Schur factorization of an upper
// Hessenberg matrix.

slahqr :: #bake_arguments _lahqr_real(RT=float32);
dlahqr :: #bake_arguments _lahqr_real(RT=float64);
clahqr :: #bake_arguments _lahqr_complex(RT=float32, MT=C32);
zlahqr :: #bake_arguments _lahqr_complex(RT=float64, MT=C64);

#scope_file
_lahqr_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_lahqr_cabs1 :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return abs(real(x)) + abs(imag(x));
}

_lahqr_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_lahqr_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_lahqr_larfg_real :: ($$n: int, alpha: *$RT, x: *$X/ArrayType1D, $$incx: int, tau: *RT, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        slarfg(n, alpha, x, incx, tau);
    } else {
        dlarfg(n, alpha, x, incx, tau);
    }
}

_lahqr_larfg_complex :: ($$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *$MT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        clarfg(n, alpha, x, incx, tau);
    } else {
        zlarfg(n, alpha, x, incx, tau);
    }
}

_lahqr_lanv2 :: (
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

_lahqr_ladiv :: ($$x: $MT, $$y: MT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        return cladiv(x, y);
    } else {
        return zladiv(x, y);
    }
}

_lahqr_real :: (
    $$wantt: bool,
    $$wantz: bool,
    $$n: int,
    $$ilo: int,
    $$ihi: int,
    h: *$H/ArrayType2D,
    $$ldh: int,
    wr: *$WR/ArrayType1D,
    wi: *$WI/ArrayType1D,
    $$iloz: int,
    $$ihiz: int,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           H.data_type == RT &&
           WR.data_type == RT &&
           WI.data_type == RT &&
           Z.data_type == RT;
} {
    DAT1 :: #run ncast(RT, 3.0 / 4.0);
    DAT2 :: #run ncast(RT, -0.4375);
    KEXSH :: 10;
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);

    aa, ab, ba, bb, cs, det, h11, h12, h21, h21s, h22, rt1i, rt1r, rt2i, rt2r, rtdisc, s, safmax, safmin, smlnum,
    sn, sum, t1, t2, t3, tr, tst, ulp, v2, v3 : RT = ---;
    i, i1, i2, its, itmax, j, k, l, m, nh, nr, nz, kdefl : int = ---;
    split : bool = ---;
    v_storage: [3] RT;
    v := array1d_view(v_storage[0..]);
    larfg_tail := array1d_view(v_storage[1..]);

    info.* = 0;
    if n == 0 then return;

    if ilo == ihi {
        fset(wr, ilo, fget(h.*, ilo, ilo));
        fset(wi, ilo, zero_r);
        return;
    }

    for j:ilo..ihi - 3 {
        fset(h, j + 2, j, zero_r);
        fset(h, j + 3, j, zero_r);
    }
    if ilo <= ihi - 2 {
        fset(h, ihi, ihi - 2, zero_r);
    }

    nh = ihi - ilo + 1;
    nz = ihiz - iloz + 1;

    safmin = _lahqr_lamch_safe_minimum(RT);
    safmax = one_r / safmin;
    ulp = _lahqr_lamch_precision(RT);
    smlnum = safmin*(ncast(RT, nh) / ulp);
    _ = safmax;

    if wantt {
        i1 = 1;
        i2 = n;
    } else {
        i1 = ilo;
        i2 = ihi;
    }

    itmax = 30*max(10, nh);
    kdefl = 0;
    i = ihi;

    while true {
        l = ilo;
        if i < ilo break;

        split = false;
        for its:0..itmax {
            k = l;
            for <kk:i..l + 1 {
                k = kk;
                if abs(fget(h.*, k, k - 1)) <= smlnum break;

                tst = abs(fget(h.*, k - 1, k - 1)) + abs(fget(h.*, k, k));
                if tst == zero_r {
                    if k - 2 >= ilo {
                        tst += abs(fget(h.*, k - 1, k - 2));
                    }
                    if k + 1 <= ihi {
                        tst += abs(fget(h.*, k + 1, k));
                    }
                }

                if abs(fget(h.*, k, k - 1)) <= ulp*tst {
                    ab = max(abs(fget(h.*, k, k - 1)), abs(fget(h.*, k - 1, k)));
                    ba = min(abs(fget(h.*, k, k - 1)), abs(fget(h.*, k - 1, k)));
                    aa = max(abs(fget(h.*, k, k)),
                              abs(fget(h.*, k - 1, k - 1) - fget(h.*, k, k)));
                    bb = min(abs(fget(h.*, k, k)),
                              abs(fget(h.*, k - 1, k - 1) - fget(h.*, k, k)));
                    s = aa + ab;
                    if ba*(ab / s) <= max(smlnum, ulp*(bb*(aa / s))) break;
                }
            }

            l = k;
            if l > ilo {
                fset(h, l, l - 1, zero_r);
            }

            if l >= i - 1 {
                split = true;
                break;
            }

            kdefl += 1;

            if !wantt {
                i1 = l;
                i2 = i;
            }

            if mod(kdefl, 2*KEXSH) == 0 {
                s = abs(fget(h.*, i, i - 1)) + abs(fget(h.*, i - 1, i - 2));
                h11 = DAT1*s + fget(h.*, i, i);
                h12 = DAT2*s;
                h21 = s;
                h22 = h11;
            } else if mod(kdefl, KEXSH) == 0 {
                s = abs(fget(h.*, l + 1, l)) + abs(fget(h.*, l + 2, l + 1));
                h11 = DAT1*s + fget(h.*, l, l);
                h12 = DAT2*s;
                h21 = s;
                h22 = h11;
            } else {
                h11 = fget(h.*, i - 1, i - 1);
                h21 = fget(h.*, i, i - 1);
                h12 = fget(h.*, i - 1, i);
                h22 = fget(h.*, i, i);
            }

            s = abs(h11) + abs(h12) + abs(h21) + abs(h22);
            if s == zero_r {
                rt1r = zero_r;
                rt1i = zero_r;
                rt2r = zero_r;
                rt2i = zero_r;
            } else {
                h11 /= s;
                h21 /= s;
                h12 /= s;
                h22 /= s;
                tr = (h11 + h22) / two_r;
                det = (h11 - tr)*(h22 - tr) - h12*h21;
                rtdisc = sqrt(abs(det));
                if det >= zero_r {
                    rt1r = tr*s;
                    rt2r = rt1r;
                    rt1i = rtdisc*s;
                    rt2i = -rt1i;
                } else {
                    rt1r = tr + rtdisc;
                    rt2r = tr - rtdisc;
                    if abs(rt1r - h22) <= abs(rt2r - h22) {
                        rt1r *= s;
                        rt2r = rt1r;
                    } else {
                        rt2r *= s;
                        rt1r = rt2r;
                    }
                    rt1i = zero_r;
                    rt2i = zero_r;
                }
            }

            m = l;
            for <mm:i - 2..l {
                m = mm;
                h21s = fget(h.*, m + 1, m);
                s = abs(fget(h.*, m, m) - rt2r) + abs(rt2i) + abs(h21s);
                h21s = fget(h.*, m + 1, m) / s;
                fset(v, 1, h21s*fget(h.*, m, m + 1) +
                           (fget(h.*, m, m) - rt1r)*((fget(h.*, m, m) - rt2r) / s) -
                           rt1i*(rt2i / s));
                fset(v, 2, h21s*(fget(h.*, m, m) + fget(h.*, m + 1, m + 1) - rt1r - rt2r));
                fset(v, 3, h21s*fget(h.*, m + 2, m + 1));
                s = abs(fget(v, 1)) + abs(fget(v, 2)) + abs(fget(v, 3));
                fset(v, 1, fget(v, 1) / s);
                fset(v, 2, fget(v, 2) / s);
                fset(v, 3, fget(v, 3) / s);
                if m == l break;
                if abs(fget(h.*, m, m - 1))*(abs(fget(v, 2)) + abs(fget(v, 3))) <=
                   ulp*abs(fget(v, 1))*(abs(fget(h.*, m - 1, m - 1)) + abs(fget(h.*, m, m)) + abs(fget(h.*, m + 1, m + 1))) {
                    break;
                }
            }

            for k:m..i - 1 {
                nr = min(3, i - k + 1);
                if k > m {
                    for jj:1..nr {
                        fset(v, jj, fget(h.*, k + jj - 1, k - 1));
                    }
                }

                _lahqr_larfg_real(nr, v.data, *larfg_tail, 1, *t1, RT);

                if k > m {
                    fset(h, k, k - 1, fget(v, 1));
                    fset(h, k + 1, k - 1, zero_r);
                    if k < i - 1 {
                        fset(h, k + 2, k - 1, zero_r);
                    }
                } else if m > l {
                    fset(h, k, k - 1, fget(h.*, k, k - 1)*(one_r - t1));
                }

                v2 = fget(v, 2);
                t2 = t1*v2;

                if nr == 3 {
                    v3 = fget(v, 3);
                    t3 = t1*v3;

                    for j:k..i2 {
                        sum = fget(h.*, k, j) + v2*fget(h.*, k + 1, j) + v3*fget(h.*, k + 2, j);
                        fset(h, k, j, fget(h.*, k, j) - sum*t1);
                        fset(h, k + 1, j, fget(h.*, k + 1, j) - sum*t2);
                        fset(h, k + 2, j, fget(h.*, k + 2, j) - sum*t3);
                    }

                    for j:i1..min(k + 3, i) {
                        sum = fget(h.*, j, k) + v2*fget(h.*, j, k + 1) + v3*fget(h.*, j, k + 2);
                        fset(h, j, k, fget(h.*, j, k) - sum*t1);
                        fset(h, j, k + 1, fget(h.*, j, k + 1) - sum*t2);
                        fset(h, j, k + 2, fget(h.*, j, k + 2) - sum*t3);
                    }

                    if wantz {
                        for j:iloz..ihiz {
                            sum = fget(z.*, j, k) + v2*fget(z.*, j, k + 1) + v3*fget(z.*, j, k + 2);
                            fset(z, j, k, fget(z.*, j, k) - sum*t1);
                            fset(z, j, k + 1, fget(z.*, j, k + 1) - sum*t2);
                            fset(z, j, k + 2, fget(z.*, j, k + 2) - sum*t3);
                        }
                    }
                } else if nr == 2 {
                    for j:k..i2 {
                        sum = fget(h.*, k, j) + v2*fget(h.*, k + 1, j);
                        fset(h, k, j, fget(h.*, k, j) - sum*t1);
                        fset(h, k + 1, j, fget(h.*, k + 1, j) - sum*t2);
                    }

                    for j:i1..i {
                        sum = fget(h.*, j, k) + v2*fget(h.*, j, k + 1);
                        fset(h, j, k, fget(h.*, j, k) - sum*t1);
                        fset(h, j, k + 1, fget(h.*, j, k + 1) - sum*t2);
                    }

                    if wantz {
                        for j:iloz..ihiz {
                            sum = fget(z.*, j, k) + v2*fget(z.*, j, k + 1);
                            fset(z, j, k, fget(z.*, j, k) - sum*t1);
                            fset(z, j, k + 1, fget(z.*, j, k + 1) - sum*t2);
                        }
                    }
                }
            }
        }

        if !split {
            info.* = i;
            return;
        }

        if l == i {
            fset(wr, i, fget(h.*, i, i));
            fset(wi, i, zero_r);
        } else if l == i - 1 {
            _lahqr_lanv2(
                h.data + (i - 2)*ldh + (i - 2),
                h.data + (i - 1)*ldh + (i - 2),
                h.data + (i - 2)*ldh + (i - 1),
                h.data + (i - 1)*ldh + (i - 1),
                wr.data + (i - 2), wi.data + (i - 2), wr.data + (i - 1), wi.data + (i - 1),
                *cs, *sn, RT
            );

            if wantt {
                if i2 > i {
                    for j:i + 1..i2 {
                        aa = fget(h.*, i - 1, j);
                        bb = fget(h.*, i, j);
                        fset(h, i - 1, j, cs*aa + sn*bb);
                        fset(h, i, j, cs*bb - sn*aa);
                    }
                }
                for j:i1..i - 2 {
                    aa = fget(h.*, j, i - 1);
                    bb = fget(h.*, j, i);
                    fset(h, j, i - 1, cs*aa + sn*bb);
                    fset(h, j, i, cs*bb - sn*aa);
                }
            }

            if wantz {
                for j:iloz..ihiz {
                    aa = fget(z.*, j, i - 1);
                    bb = fget(z.*, j, i);
                    fset(z, j, i - 1, cs*aa + sn*bb);
                    fset(z, j, i, cs*bb - sn*aa);
                }
            }
        }

        kdefl = 0;
        i = l - 1;
    }
}

_lahqr_complex :: (
    $$wantt: bool,
    $$wantz: bool,
    $$n: int,
    $$ilo: int,
    $$ihi: int,
    h: *$H/ArrayType2D,
    $$ldh: int,
    w: *$W/ArrayType1D,
    $$iloz: int,
    $$ihiz: int,
    z: *$Z/ArrayType2D,
    $$ldz: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           H.data_type == MT &&
           W.data_type == MT &&
           Z.data_type == MT;
} {
    DAT1 :: #run ncast(RT, 3.0 / 4.0);
    KEXSH :: 10;
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    half_r :: #run ncast(RT, 0.5);
    zero_t :: #run zero(MT);
    one_t  :: #run _lahqr_complex_from_real(one_r, RT, MT);

    cdum, h11, h11s, h22, sc, sum, t, t1, temp, u, v2, x, y : MT = ---;
    aa, ab, ba, bb, h10, h21, rtemp, s, safmax, safmin, smlnum, sx, t2, tst, ulp : RT = ---;
    i, i1, i2, its, itmax, j, jhi, jlo, k, l, m, nh, nz, kdefl : int = ---;
    split : bool = ---;
    v_storage: [2] MT;
    v := array1d_view(v_storage[0..]);
    larfg_tail := array1d_view(v_storage[1..]);

    _ = cdum;
    info.* = 0;

    if n == 0 then return;
    if ilo == ihi {
        fset(w, ilo, fget(h.*, ilo, ilo));
        return;
    }

    for j:ilo..ihi - 3 {
        fset(h, j + 2, j, zero_t);
        fset(h, j + 3, j, zero_t);
    }
    if ilo <= ihi - 2 {
        fset(h, ihi, ihi - 2, zero_t);
    }

    if wantt {
        jlo = 1;
        jhi = n;
    } else {
        jlo = ilo;
        jhi = ihi;
    }

    for i:ilo + 1..ihi {
        if imag(fget(h.*, i, i - 1)) != zero_r {
            sc = fget(h.*, i, i - 1) / _lahqr_complex_from_real(_lahqr_cabs1(fget(h.*, i, i - 1), RT, MT), RT, MT);
            sc = conjugate(sc) / _lahqr_complex_from_real(abs(sc), RT, MT);
            fset(h, i, i - 1, _lahqr_complex_from_real(abs(fget(h.*, i, i - 1)), RT, MT));

            for j:i..jhi {
                fset(h, i, j, sc*fget(h.*, i, j));
            }
            for j:jlo..min(jhi, i + 1) {
                fset(h, j, i, conjugate(sc)*fget(h.*, j, i));
            }
            if wantz {
                for j:iloz..ihiz {
                    fset(z, j, i, conjugate(sc)*fget(z.*, j, i));
                }
            }
        }
    }

    nh = ihi - ilo + 1;
    nz = ihiz - iloz + 1;
    safmin = _lahqr_lamch_safe_minimum(RT);
    safmax = one_r / safmin;
    ulp = _lahqr_lamch_precision(RT);
    smlnum = safmin*(ncast(RT, nh) / ulp);
    _ = safmax;

    if wantt {
        i1 = 1;
        i2 = n;
    } else {
        i1 = ilo;
        i2 = ihi;
    }

    itmax = 30*max(10, nh);
    kdefl = 0;
    i = ihi;

    while true {
        if i < ilo break;

        l = ilo;
        split = false;
        for its:0..itmax {
            k = l;
            for <kk:i..l + 1 {
                k = kk;
                if _lahqr_cabs1(fget(h.*, k, k - 1), RT, MT) <= smlnum break;

                tst = _lahqr_cabs1(fget(h.*, k - 1, k - 1), RT, MT) + _lahqr_cabs1(fget(h.*, k, k), RT, MT);
                if tst == zero_r {
                    if k - 2 >= ilo {
                        tst += abs(real(fget(h.*, k - 1, k - 2)));
                    }
                    if k + 1 <= ihi {
                        tst += abs(real(fget(h.*, k + 1, k)));
                    }
                }

                if abs(real(fget(h.*, k, k - 1))) <= ulp*tst {
                    ab = max(_lahqr_cabs1(fget(h.*, k, k - 1), RT, MT), _lahqr_cabs1(fget(h.*, k - 1, k), RT, MT));
                    ba = min(_lahqr_cabs1(fget(h.*, k, k - 1), RT, MT), _lahqr_cabs1(fget(h.*, k - 1, k), RT, MT));
                    aa = max(_lahqr_cabs1(fget(h.*, k, k), RT, MT),
                             _lahqr_cabs1(fget(h.*, k - 1, k - 1) - fget(h.*, k, k), RT, MT));
                    bb = min(_lahqr_cabs1(fget(h.*, k, k), RT, MT),
                             _lahqr_cabs1(fget(h.*, k - 1, k - 1) - fget(h.*, k, k), RT, MT));
                    s = aa + ab;
                    if ba*(ab / s) <= max(smlnum, ulp*(bb*(aa / s))) break;
                }
            }

            l = k;
            if l > ilo {
                fset(h, l, l - 1, zero_t);
            }

            if l >= i {
                split = true;
                break;
            }

            kdefl += 1;
            if !wantt {
                i1 = l;
                i2 = i;
            }

            if mod(kdefl, 2*KEXSH) == 0 {
                s = DAT1*abs(real(fget(h.*, i, i - 1)));
                t = _lahqr_complex_from_real(s, RT, MT) + fget(h.*, i, i);
            } else if mod(kdefl, KEXSH) == 0 {
                s = DAT1*abs(real(fget(h.*, l + 1, l)));
                t = _lahqr_complex_from_real(s, RT, MT) + fget(h.*, l, l);
            } else {
                t = fget(h.*, i, i);
                u = sqrt(fget(h.*, i - 1, i))*sqrt(fget(h.*, i, i - 1));
                s = _lahqr_cabs1(u, RT, MT);
                if s != zero_r {
                    x = _lahqr_complex_from_real(half_r, RT, MT)*(fget(h.*, i - 1, i - 1) - t);
                    sx = _lahqr_cabs1(x, RT, MT);
                    s = max(s, _lahqr_cabs1(x, RT, MT));
                    y = _lahqr_complex_from_real(s, RT, MT)*sqrt((x / _lahqr_complex_from_real(s, RT, MT))*(x / _lahqr_complex_from_real(s, RT, MT)) +
                                                                 (u / _lahqr_complex_from_real(s, RT, MT))*(u / _lahqr_complex_from_real(s, RT, MT)));
                    if sx > zero_r {
                        if real(x / _lahqr_complex_from_real(sx, RT, MT))*real(y) +
                           imag(x / _lahqr_complex_from_real(sx, RT, MT))*imag(y) < zero_r {
                            y = -y;
                        }
                    }
                    t = t - u*_lahqr_ladiv(u, x + y, RT, MT);
                }
            }

            m = l;
            for <mm:i - 1..l + 1 {
                m = mm;
                h11 = fget(h.*, m, m);
                h22 = fget(h.*, m + 1, m + 1);
                h11s = h11 - t;
                h21 = real(fget(h.*, m + 1, m));
                s = _lahqr_cabs1(h11s, RT, MT) + abs(h21);
                h11s /= _lahqr_complex_from_real(s, RT, MT);
                h21 /= s;
                fset(v, 1, h11s);
                fset(v, 2, _lahqr_complex_from_real(h21, RT, MT));
                h10 = real(fget(h.*, m, m - 1));
                if abs(h10)*abs(h21) <= ulp*(_lahqr_cabs1(h11s, RT, MT)*(_lahqr_cabs1(h11, RT, MT) + _lahqr_cabs1(h22, RT, MT))) {
                    break;
                }
            }

            if m == l {
                h11 = fget(h.*, l, l);
                h22 = fget(h.*, l + 1, l + 1);
                h11s = h11 - t;
                h21 = real(fget(h.*, l + 1, l));
                s = _lahqr_cabs1(h11s, RT, MT) + abs(h21);
                h11s /= _lahqr_complex_from_real(s, RT, MT);
                h21 /= s;
                fset(v, 1, h11s);
                fset(v, 2, _lahqr_complex_from_real(h21, RT, MT));
            }

            for k:m..i - 1 {
                if k > m {
                    fset(v, 1, fget(h.*, k, k - 1));
                    fset(v, 2, fget(h.*, k + 1, k - 1));
                }

                _lahqr_larfg_complex(2, v.data, *larfg_tail, 1, *t1, RT, MT);

                if k > m {
                    fset(h, k, k - 1, fget(v, 1));
                    fset(h, k + 1, k - 1, zero_t);
                }

                v2 = fget(v, 2);
                t2 = real(t1*v2);

                for j:k..i2 {
                    sum = conjugate(t1)*fget(h.*, k, j) + _lahqr_complex_from_real(t2, RT, MT)*fget(h.*, k + 1, j);
                    fset(h, k, j, fget(h.*, k, j) - sum);
                    fset(h, k + 1, j, fget(h.*, k + 1, j) - sum*v2);
                }

                for j:i1..min(k + 2, i) {
                    sum = t1*fget(h.*, j, k) + _lahqr_complex_from_real(t2, RT, MT)*fget(h.*, j, k + 1);
                    fset(h, j, k, fget(h.*, j, k) - sum);
                    fset(h, j, k + 1, fget(h.*, j, k + 1) - sum*conjugate(v2));
                }

                if wantz {
                    for j:iloz..ihiz {
                        sum = t1*fget(z.*, j, k) + _lahqr_complex_from_real(t2, RT, MT)*fget(z.*, j, k + 1);
                        fset(z, j, k, fget(z.*, j, k) - sum);
                        fset(z, j, k + 1, fget(z.*, j, k + 1) - sum*conjugate(v2));
                    }
                }

                if k == m && m > l {
                    temp = one_t - t1;
                    temp /= _lahqr_complex_from_real(abs(temp), RT, MT);
                    fset(h, m + 1, m, fget(h.*, m + 1, m)*conjugate(temp));
                    if m + 2 <= i {
                        fset(h, m + 2, m + 1, fget(h.*, m + 2, m + 1)*temp);
                    }
                    for j:m..i {
                        if j != m + 1 {
                            for jj:j + 1..i2 {
                                fset(h, j, jj, temp*fget(h.*, j, jj));
                            }
                            for jj:i1..j - 1 {
                                fset(h, jj, j, conjugate(temp)*fget(h.*, jj, j));
                            }
                            if wantz {
                                for jj:iloz..ihiz {
                                    fset(z, jj, j, conjugate(temp)*fget(z.*, jj, j));
                                }
                            }
                        }
                    }
                }
            }

            temp = fget(h.*, i, i - 1);
            if imag(temp) != zero_r {
                rtemp = abs(temp);
                fset(h, i, i - 1, _lahqr_complex_from_real(rtemp, RT, MT));
                temp /= _lahqr_complex_from_real(rtemp, RT, MT);

                for j:i + 1..i2 {
                    fset(h, i, j, conjugate(temp)*fget(h.*, i, j));
                }
                for j:i1..i - 1 {
                    fset(h, j, i, temp*fget(h.*, j, i));
                }
                if wantz {
                    for j:iloz..ihiz {
                        fset(z, j, i, temp*fget(z.*, j, i));
                    }
                }
            }
        }

        if !split {
            info.* = i;
            return;
        }

        fset(w, i, fget(h.*, i, i));
        kdefl = 0;
        i = l - 1;
    }
}
