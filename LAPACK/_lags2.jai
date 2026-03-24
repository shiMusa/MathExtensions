// ? LAPACK *lags2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slags2/clags2/dlags2/zlags2
// *LAGS2 computes 2-by-2 orthogonal/unitary matrices U, V, and Q,
// and applies them so the rows of the transformed A and B are parallel.

slags2 :: #bake_arguments _lags2_real(RT=float32);
dlags2 :: #bake_arguments _lags2_real(RT=float64);
clags2 :: #bake_arguments _lags2_complex(RT=float32, MT=C32);
zlags2 :: #bake_arguments _lags2_complex(RT=float64, MT=C64);

#scope_file
_lags2_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_lags2_lasv2 :: (
    $$f: $RT, $$g: RT, $$h: RT,
    ssmin: *RT, ssmax: *RT, snr: *RT, csr: *RT, snl: *RT, csl: *RT,
    $RT: Type
) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        slasv2(f, g, h, ssmin, ssmax, snr, csr, snl, csl);
    } else {
        dlasv2(f, g, h, ssmin, ssmax, snr, csr, snl, csl);
    }
}

_lags2_lartg_real :: ($$f: $RT, $$g: RT, c: *RT, s: *RT, r: *RT, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        slartg(f, g, c, s, r);
    } else {
        dlartg(f, g, c, s, r);
    }
}

_lags2_lartg_complex :: ($$f: $MT, $$g: MT, c: *$RT, s: *MT, r: *MT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 {
        clartg(f, g, c, s, r);
    } else {
        zlartg(f, g, c, s, r);
    }
}

_lags2_real :: (
    $$upper: bool,
    $$a1: $RT,
    $$a2: RT,
    $$a3: RT,
    $$b1: RT,
    $$b2: RT,
    $$b3: RT,
    csu: *RT,
    snu: *RT,
    csv: *RT,
    snv: *RT,
    csq: *RT,
    snq: *RT,
    $RT: Type
) #modify {
    return RT == float32 || RT == float64;
} {
    zero_r :: #run zero(RT);

    a, aua11, aua12, aua21, aua22, avb11, avb12, avb21, avb22, csl, csr, d, s1, s2, snl, snr, ua11r, ua22r, vb11r, vb22r, b, c, r, ua11, ua12, ua21, ua22, vb11, vb12, vb21, vb22 : RT = ---;

    if upper {
        a = a1*b3;
        d = a3*b1;
        b = a2*b1 - a1*b2;

        _lags2_lasv2(a, b, d, *s1, *s2, *snr, *csr, *snl, *csl, RT);

        if abs(csl) >= abs(snl) || abs(csr) >= abs(snr) {
            ua11r = csl*a1;
            ua12 = csl*a2 + snl*a3;

            vb11r = csr*b1;
            vb12 = csr*b2 + snr*b3;

            aua12 = abs(csl)*abs(a2) + abs(snl)*abs(a3);
            avb12 = abs(csr)*abs(b2) + abs(snr)*abs(b3);

            if (abs(ua11r) + abs(ua12)) != zero_r {
                if aua12 / (abs(ua11r) + abs(ua12)) <= avb12 / (abs(vb11r) + abs(vb12)) {
                    _lags2_lartg_real(-ua11r, ua12, csq, snq, *r, RT);
                } else {
                    _lags2_lartg_real(-vb11r, vb12, csq, snq, *r, RT);
                }
            } else {
                _lags2_lartg_real(-vb11r, vb12, csq, snq, *r, RT);
            }

            csu.* = csl;
            snu.* = -snl;
            csv.* = csr;
            snv.* = -snr;
        } else {
            ua21 = -snl*a1;
            ua22 = -snl*a2 + csl*a3;

            vb21 = -snr*b1;
            vb22 = -snr*b2 + csr*b3;

            aua22 = abs(snl)*abs(a2) + abs(csl)*abs(a3);
            avb22 = abs(snr)*abs(b2) + abs(csr)*abs(b3);

            if (abs(ua21) + abs(ua22)) != zero_r {
                if aua22 / (abs(ua21) + abs(ua22)) <= avb22 / (abs(vb21) + abs(vb22)) {
                    _lags2_lartg_real(-ua21, ua22, csq, snq, *r, RT);
                } else {
                    _lags2_lartg_real(-vb21, vb22, csq, snq, *r, RT);
                }
            } else {
                _lags2_lartg_real(-vb21, vb22, csq, snq, *r, RT);
            }

            csu.* = snl;
            snu.* = csl;
            csv.* = snr;
            snv.* = csr;
        }
    } else {
        a = a1*b3;
        d = a3*b1;
        c = a2*b3 - a3*b2;

        _lags2_lasv2(a, c, d, *s1, *s2, *snr, *csr, *snl, *csl, RT);

        if abs(csr) >= abs(snr) || abs(csl) >= abs(snl) {
            ua21 = -snr*a1 + csr*a2;
            ua22r = csr*a3;

            vb21 = -snl*b1 + csl*b2;
            vb22r = csl*b3;

            aua21 = abs(snr)*abs(a1) + abs(csr)*abs(a2);
            avb21 = abs(snl)*abs(b1) + abs(csl)*abs(b2);

            if (abs(ua21) + abs(ua22r)) != zero_r {
                if aua21 / (abs(ua21) + abs(ua22r)) <= avb21 / (abs(vb21) + abs(vb22r)) {
                    _lags2_lartg_real(ua22r, ua21, csq, snq, *r, RT);
                } else {
                    _lags2_lartg_real(vb22r, vb21, csq, snq, *r, RT);
                }
            } else {
                _lags2_lartg_real(vb22r, vb21, csq, snq, *r, RT);
            }

            csu.* = csr;
            snu.* = -snr;
            csv.* = csl;
            snv.* = -snl;
        } else {
            ua11 = csr*a1 + snr*a2;
            ua12 = snr*a3;

            vb11 = csl*b1 + snl*b2;
            vb12 = snl*b3;

            aua11 = abs(csr)*abs(a1) + abs(snr)*abs(a2);
            avb11 = abs(csl)*abs(b1) + abs(snl)*abs(b2);

            if (abs(ua11) + abs(ua12)) != zero_r {
                if aua11 / (abs(ua11) + abs(ua12)) <= avb11 / (abs(vb11) + abs(vb12)) {
                    _lags2_lartg_real(ua12, ua11, csq, snq, *r, RT);
                } else {
                    _lags2_lartg_real(vb12, vb11, csq, snq, *r, RT);
                }
            } else {
                _lags2_lartg_real(vb12, vb11, csq, snq, *r, RT);
            }

            csu.* = snr;
            snu.* = csr;
            csv.* = snl;
            snv.* = csl;
        }
    }
}

_lags2_complex :: (
    $$upper: bool,
    $$a1: $RT,
    $$a2: $MT,
    $$a3: RT,
    $$b1: RT,
    $$b2: MT,
    $$b3: RT,
    csu: *RT,
    snu: *MT,
    csv: *RT,
    snv: *MT,
    csq: *RT,
    snq: *MT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    one_t  :: #run _lags2_complex_from_real(one_r, RT, MT);

    a, aua11, aua12, aua21, aua22, avb12, avb11, avb21, avb22, csl, csr, d, fb, fc, s1, s2, snl, snr, ua11r, ua22r, vb11r, vb22r : RT = ---;
    b, c, d1, r, t, ua11, ua12, ua21, ua22, vb11, vb12, vb21, vb22 : MT = ---;

    if upper {
        a = a1*b3;
        d = a3*b1;
        b = a2*b1 - _lags2_complex_from_real(a1, RT, MT)*b2;
        fb = abs(b);

        d1 = one_t;
        if fb != zero_r {
            d1 = b / fb;
        }

        _lags2_lasv2(a, fb, d, *s1, *s2, *snr, *csr, *snl, *csl, RT);

        if abs(csl) >= abs(snl) || abs(csr) >= abs(snr) {
            ua11r = csl*a1;
            ua12 = csl*a2 + d1*(snl*a3);

            vb11r = csr*b1;
            vb12 = csr*b2 + d1*(snr*b3);

            aua12 = abs(csl)*_abs1(a2) + abs(snl)*abs(a3);
            avb12 = abs(csr)*_abs1(b2) + abs(snr)*abs(b3);

            if (abs(ua11r) + _abs1(ua12)) == zero_r {
                _lags2_lartg_complex(-_lags2_complex_from_real(vb11r, RT, MT), conjugate(vb12), csq, snq, *r, RT, MT);
            } else if (abs(vb11r) + _abs1(vb12)) == zero_r {
                _lags2_lartg_complex(-_lags2_complex_from_real(ua11r, RT, MT), conjugate(ua12), csq, snq, *r, RT, MT);
            } else if aua12 / (abs(ua11r) + _abs1(ua12)) <= avb12 / (abs(vb11r) + _abs1(vb12)) {
                _lags2_lartg_complex(-_lags2_complex_from_real(ua11r, RT, MT), conjugate(ua12), csq, snq, *r, RT, MT);
            } else {
                _lags2_lartg_complex(-_lags2_complex_from_real(vb11r, RT, MT), conjugate(vb12), csq, snq, *r, RT, MT);
            }

            csu.* = csl;
            snu.* = -d1*snl;
            csv.* = csr;
            snv.* = -d1*snr;
        } else {
            ua21 = -conjugate(d1)*(snl*a1);
            ua22 = -conjugate(d1)*(snl*a2) + _lags2_complex_from_real(csl*a3, RT, MT);

            vb21 = -conjugate(d1)*(snr*b1);
            vb22 = -conjugate(d1)*(snr*b2) + _lags2_complex_from_real(csr*b3, RT, MT);

            aua22 = abs(snl)*_abs1(a2) + abs(csl)*abs(a3);
            avb22 = abs(snr)*_abs1(b2) + abs(csr)*abs(b3);

            if (_abs1(ua21) + _abs1(ua22)) == zero_r {
                _lags2_lartg_complex(-conjugate(vb21), conjugate(vb22), csq, snq, *r, RT, MT);
            } else if (_abs1(vb21) + _abs1(vb22)) == zero_r {
                _lags2_lartg_complex(-conjugate(ua21), conjugate(ua22), csq, snq, *r, RT, MT);
            } else if aua22 / (_abs1(ua21) + _abs1(ua22)) <= avb22 / (_abs1(vb21) + _abs1(vb22)) {
                _lags2_lartg_complex(-conjugate(ua21), conjugate(ua22), csq, snq, *r, RT, MT);
            } else {
                _lags2_lartg_complex(-conjugate(vb21), conjugate(vb22), csq, snq, *r, RT, MT);
            }

            csu.* = snl;
            snu.* = d1*csl;
            csv.* = snr;
            snv.* = d1*csr;
        }
    } else {
        a = a1*b3;
        d = a3*b1;
        c = a2*b3 - _lags2_complex_from_real(a3, RT, MT)*b2;
        fc = abs(c);

        d1 = one_t;
        if fc != zero_r {
            d1 = c / fc;
        }

        _lags2_lasv2(a, fc, d, *s1, *s2, *snr, *csr, *snl, *csl, RT);

        if abs(csr) >= abs(snr) || abs(csl) >= abs(snl) {
            ua21 = -d1*(snr*a1) + csr*a2;
            ua22r = csr*a3;

            vb21 = -d1*(snl*b1) + csl*b2;
            vb22r = csl*b3;

            aua21 = abs(snr)*abs(a1) + abs(csr)*_abs1(a2);
            avb21 = abs(snl)*abs(b1) + abs(csl)*_abs1(b2);

            if (_abs1(ua21) + abs(ua22r)) == zero_r {
                _lags2_lartg_complex(_lags2_complex_from_real(vb22r, RT, MT), vb21, csq, snq, *r, RT, MT);
            } else if (_abs1(vb21) + abs(vb22r)) == zero_r {
                _lags2_lartg_complex(_lags2_complex_from_real(ua22r, RT, MT), ua21, csq, snq, *r, RT, MT);
            } else if aua21 / (_abs1(ua21) + abs(ua22r)) <= avb21 / (_abs1(vb21) + abs(vb22r)) {
                _lags2_lartg_complex(_lags2_complex_from_real(ua22r, RT, MT), ua21, csq, snq, *r, RT, MT);
            } else {
                _lags2_lartg_complex(_lags2_complex_from_real(vb22r, RT, MT), vb21, csq, snq, *r, RT, MT);
            }

            csu.* = csr;
            snu.* = -conjugate(d1)*snr;
            csv.* = csl;
            snv.* = -conjugate(d1)*snl;
        } else {
            ua11 = csr*a1 + conjugate(d1)*(snr*a2);
            ua12 = conjugate(d1)*(snr*a3);

            vb11 = csl*b1 + conjugate(d1)*(snl*b2);
            vb12 = conjugate(d1)*(snl*b3);

            aua11 = abs(csr)*abs(a1) + abs(snr)*_abs1(a2);
            avb11 = abs(csl)*abs(b1) + abs(snl)*_abs1(b2);

            if (_abs1(ua11) + _abs1(ua12)) == zero_r {
                _lags2_lartg_complex(vb12, vb11, csq, snq, *r, RT, MT);
            } else if (_abs1(vb11) + _abs1(vb12)) == zero_r {
                _lags2_lartg_complex(ua12, ua11, csq, snq, *r, RT, MT);
            } else if aua11 / (_abs1(ua11) + _abs1(ua12)) <= avb11 / (_abs1(vb11) + _abs1(vb12)) {
                _lags2_lartg_complex(ua12, ua11, csq, snq, *r, RT, MT);
            } else {
                _lags2_lartg_complex(vb12, vb11, csq, snq, *r, RT, MT);
            }

            csu.* = snr;
            snu.* = conjugate(d1)*csr;
            csv.* = snl;
            snv.* = conjugate(d1)*csl;
        }
    }
}
