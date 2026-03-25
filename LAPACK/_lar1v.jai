// ? LAPACK *lar1v()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slar1v/dlar1v/clar1v/zlar1v
// *LAR1V computes an FP vector for a tridiagonal eigensolver.

slar1v :: #bake_arguments _lar1v(RT=float32, MT=float32);
dlar1v :: #bake_arguments _lar1v(RT=float64, MT=float64);
clar1v :: #bake_arguments _lar1v(RT=float32, MT=C32);
zlar1v :: #bake_arguments _lar1v(RT=float64, MT=C64);

#scope_file
_lar1v_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "P");
    } else {
        return dlamch(#char "P");
    }
}

_lar1v_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_lar1v_absz :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == RT {
        return abs(z);
    } else {
        zr := real(z);
        zi := imag(z);
        return sqrt(zr*zr + zi*zi);
    }
}

_lar1v_ztz_term :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == RT {
        return z*z;
    } else {
        return real(z*z);
    }
}

_lar1v :: (
    $$n: int, $$b1: int, $$bn: int, lambda: $RT, d: $D/ArrayType1D, l: $L/ArrayType1D, ld: $LD/ArrayType1D, lld: $LLD/ArrayType1D,
    $$pivmin: RT, $$gaptol: RT, z: *$Z/ArrayType1D, $$wantnc: bool, negcnt: *int, ztz: *RT, mingma: *RT,
    r: *int, isuppz: *$ISUPPZ/ArrayType1D, nrminv: *RT, resid: *RT, rqcorr: *RT, work: *$WORK/ArrayType1D,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           D.data_type == RT && L.data_type == RT && LD.data_type == RT && LLD.data_type == RT &&
           Z.data_type == MT && ISUPPZ.data_type == int && WORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    one_t  :: #run one(MT);
    zero_t :: #run zero(MT);

    i, indlpl, indp, inds, indumn, neg1, neg2, r1, r2 : int = ---;
    dminus, dplus, eps, s, tmp : RT = ---;
    sawnan1, sawnan2, stop_up, stop_dn : bool = ---;

    eps = _lar1v_lamch_precision(RT);

    if r.* == 0 {
        r1 = b1;
        r2 = bn;
    } else {
        r1 = r.*;
        r2 = r.*;
    }

    indlpl = 0;
    indumn = n;
    inds   = 2*n + 1;
    indp   = 3*n + 1;

    if b1 == 1 {
        fset(work, inds, zero_r);
    } else {
        fset(work, inds + b1 - 1, fget(lld, b1 - 1));
    }

    sawnan1 = false;
    neg1 = 0;
    s = fget(work, inds + b1 - 1) - lambda;
    for i:b1..r1 - 1 {
        dplus = fget(d, i) + s;
        fset(work, indlpl + i, fget(ld, i) / dplus);
        if dplus < zero_r then neg1 += 1;
        fset(work, inds + i, s*fget(work, indlpl + i)*fget(l, i));
        s = fget(work, inds + i) - lambda;
    }
    sawnan1 = _lar1v_isnan(s, RT);

    if !sawnan1 {
        for i:r1..r2 - 1 {
            dplus = fget(d, i) + s;
            fset(work, indlpl + i, fget(ld, i) / dplus);
            fset(work, inds + i, s*fget(work, indlpl + i)*fget(l, i));
            s = fget(work, inds + i) - lambda;
        }
        sawnan1 = _lar1v_isnan(s, RT);
    }

    if sawnan1 {
        neg1 = 0;
        s = fget(work, inds + b1 - 1) - lambda;
        for i:b1..r1 - 1 {
            dplus = fget(d, i) + s;
            if abs(dplus) < pivmin then dplus = -pivmin;
            fset(work, indlpl + i, fget(ld, i) / dplus);
            if dplus < zero_r then neg1 += 1;
            fset(work, inds + i, s*fget(work, indlpl + i)*fget(l, i));
            if fget(work, indlpl + i) == zero_r {
                fset(work, inds + i, fget(lld, i));
            }
            s = fget(work, inds + i) - lambda;
        }
        for i:r1..r2 - 1 {
            dplus = fget(d, i) + s;
            if abs(dplus) < pivmin then dplus = -pivmin;
            fset(work, indlpl + i, fget(ld, i) / dplus);
            fset(work, inds + i, s*fget(work, indlpl + i)*fget(l, i));
            if fget(work, indlpl + i) == zero_r {
                fset(work, inds + i, fget(lld, i));
            }
            s = fget(work, inds + i) - lambda;
        }
    }

    sawnan2 = false;
    neg2 = 0;
    fset(work, indp + bn - 1, fget(d, bn) - lambda);
    for <i:r1..bn - 1 {
        dminus = fget(lld, i) + fget(work, indp + i);
        tmp = fget(d, i) / dminus;
        if dminus < zero_r then neg2 += 1;
        fset(work, indumn + i, fget(l, i)*tmp);
        fset(work, indp + i - 1, fget(work, indp + i)*tmp - lambda);
    }
    tmp = fget(work, indp + r1 - 1);
    sawnan2 = _lar1v_isnan(tmp, RT);

    if sawnan2 {
        neg2 = 0;
        for <i:r1..bn - 1 {
            dminus = fget(lld, i) + fget(work, indp + i);
            if abs(dminus) < pivmin then dminus = -pivmin;
            tmp = fget(d, i) / dminus;
            if dminus < zero_r then neg2 += 1;
            fset(work, indumn + i, fget(l, i)*tmp);
            fset(work, indp + i - 1, fget(work, indp + i)*tmp - lambda);
            if tmp == zero_r {
                fset(work, indp + i - 1, fget(d, i) - lambda);
            }
        }
    }

    mingma.* = fget(work, inds + r1 - 1) + fget(work, indp + r1 - 1);
    if mingma.* < zero_r then neg1 += 1;
    if wantnc {
        negcnt.* = neg1 + neg2;
    } else {
        negcnt.* = -1;
    }
    if abs(mingma.*) == zero_r {
        mingma.* = eps*fget(work, inds + r1 - 1);
    }

    r.* = r1;
    for i:r1..r2 - 1 {
        tmp = fget(work, inds + i) + fget(work, indp + i);
        if tmp == zero_r {
            tmp = eps*fget(work, inds + i);
        }
        if abs(tmp) <= abs(mingma.*) {
            mingma.* = tmp;
            r.* = i + 1;
        }
    }

    fset(isuppz, 1, b1);
    fset(isuppz, 2, bn);
    fset(z, r.*, one_t);
    ztz.* = one_r;

    if !sawnan1 && !sawnan2 {
        stop_up = false;
        for <i:b1..r.* - 1 {
            fset(z, i, -(fget(work, indlpl + i)*fget(z.*, i + 1)));
            if (_lar1v_absz(fget(z.*, i), RT, MT) + _lar1v_absz(fget(z.*, i + 1), RT, MT))*abs(fget(ld, i)) < gaptol {
                fset(z, i, zero_t);
                fset(isuppz, 1, i + 1);
                stop_up = true;
            } else {
                ztz.* = ztz.* + _lar1v_ztz_term(fget(z.*, i), RT, MT);
            }
            if stop_up then break;
        }
    } else {
        stop_up = false;
        for <i:b1..r.* - 1 {
            if fget(z.*, i + 1) == zero_t {
                fset(z, i, -((fget(ld, i + 1) / fget(ld, i))*fget(z.*, i + 2)));
            } else {
                fset(z, i, -(fget(work, indlpl + i)*fget(z.*, i + 1)));
            }
            if (_lar1v_absz(fget(z.*, i), RT, MT) + _lar1v_absz(fget(z.*, i + 1), RT, MT))*abs(fget(ld, i)) < gaptol {
                fset(z, i, zero_t);
                fset(isuppz, 1, i + 1);
                stop_up = true;
            } else {
                ztz.* = ztz.* + _lar1v_ztz_term(fget(z.*, i), RT, MT);
            }
            if stop_up then break;
        }
    }

    if !sawnan1 && !sawnan2 {
        stop_dn = false;
        for i:r.*..bn - 1 {
            fset(z, i + 1, -(fget(work, indumn + i)*fget(z.*, i)));
            if (_lar1v_absz(fget(z.*, i), RT, MT) + _lar1v_absz(fget(z.*, i + 1), RT, MT))*abs(fget(ld, i)) < gaptol {
                fset(z, i + 1, zero_t);
                fset(isuppz, 2, i);
                stop_dn = true;
            } else {
                ztz.* = ztz.* + _lar1v_ztz_term(fget(z.*, i + 1), RT, MT);
            }
            if stop_dn then break;
        }
    } else {
        stop_dn = false;
        for i:r.*..bn - 1 {
            if fget(z.*, i) == zero_t {
                fset(z, i + 1, -((fget(ld, i - 1) / fget(ld, i))*fget(z.*, i - 1)));
            } else {
                fset(z, i + 1, -(fget(work, indumn + i)*fget(z.*, i)));
            }
            if (_lar1v_absz(fget(z.*, i), RT, MT) + _lar1v_absz(fget(z.*, i + 1), RT, MT))*abs(fget(ld, i)) < gaptol {
                fset(z, i + 1, zero_t);
                fset(isuppz, 2, i);
                stop_dn = true;
            } else {
                ztz.* = ztz.* + _lar1v_ztz_term(fget(z.*, i + 1), RT, MT);
            }
            if stop_dn then break;
        }
    }

    tmp = one_r / ztz.*;
    nrminv.* = sqrt(tmp);
    resid.* = abs(mingma.*)*nrminv.*;
    rqcorr.* = mingma.* * tmp;
}
