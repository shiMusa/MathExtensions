// ? LAPACK *lartg()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slartg/dlartg/clartg/zlartg
// *LARTG generates a plane rotation.

slartg :: #bake_arguments _lartg_real(RT=float32);
dlartg :: #bake_arguments _lartg_real(RT=float64);
clartg :: #bake_arguments _lartg_complex(RT=float32, MT=C32);
zlartg :: #bake_arguments _lartg_complex(RT=float64, MT=C64);

#scope_file
_lartg_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_lartg_abssq :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return real(z)*real(z) + imag(z)*imag(z);
}

_lartg_real :: (f: $RT, g: RT, cs: *RT, sn: *RT, r: *RT, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);

    #if RT == float64 {
        safmin :: SAFMIN_64;
        safmax :: SAFMAX_64;
        rtmin  :: RTMIN_64;
    } else {
        safmin :: SAFMIN_32;
        safmax :: SAFMAX_32;
        rtmin  :: RTMIN_32;
    }

    d, f1, fs, g1, gs, u, rtmax : RT = ---;

    rtmax = sqrt(safmax / two_r);
    f1 = abs(f);
    g1 = abs(g);

    if g == zero_r {
        cs.* = one_r;
        sn.* = zero_r;
        r.* = f;
    } else if f == zero_r {
        cs.* = zero_r;
        sn.* = sign(one_r, g);
        r.* = g1;
    } else if f1 > rtmin && f1 < rtmax &&
              g1 > rtmin && g1 < rtmax {
        d = sqrt(f*f + g*g);
        cs.* = f1 / d;
        r.* = sign(d, f);
        sn.* = g / r.*;
    } else {
        u = min(safmax, max(safmin, f1, g1));
        fs = f / u;
        gs = g / u;
        d = sqrt(fs*fs + gs*gs);
        cs.* = abs(fs) / d;
        r.* = sign(d, f);
        sn.* = gs / r.*;
        r.* = r.* * u;
    }
}

_lartg_complex :: (f: $MT, g: MT, cs: *$RT, sn: *MT, r: *MT, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);
    four_r :: #run ncast(RT, 4);
    zero_t :: #run zero(MT);

    #if RT == float64 {
        safmin :: SAFMIN_64;
        safmax :: SAFMAX_64;
        rtmin  :: RTMIN_64;
    } else {
        safmin :: SAFMIN_32;
        safmax :: SAFMAX_32;
        rtmin  :: RTMIN_32;
    }

    d, f1, f2, g1, g2, h2, u, v, w, rtmax : RT = ---;
    fs, gs : MT = ---;

    if g == zero_t {
        cs.* = one_r;
        sn.* = zero_t;
        r.* = f;
    } else if f == zero_t {
        cs.* = zero_r;
        if real(g) == zero_r {
            r.* = _lartg_complex_from_real(abs(imag(g)), RT, MT);
            sn.* = conjugate(g) / real(r.*);
        } else if imag(g) == zero_r {
            r.* = _lartg_complex_from_real(abs(real(g)), RT, MT);
            sn.* = conjugate(g) / real(r.*);
        } else {
            g1 = max(abs(real(g)), abs(imag(g)));
            rtmax = sqrt(safmax / two_r);
            if g1 > rtmin && g1 < rtmax {
                g2 = _lartg_abssq(g, RT, MT);
                d = sqrt(g2);
                sn.* = conjugate(g) / d;
                r.* = _lartg_complex_from_real(d, RT, MT);
            } else {
                u = min(safmax, max(safmin, g1));
                gs = g / u;
                g2 = _lartg_abssq(gs, RT, MT);
                d = sqrt(g2);
                sn.* = conjugate(gs) / d;
                r.* = _lartg_complex_from_real(d*u, RT, MT);
            }
        }
    } else {
        f1 = max(abs(real(f)), abs(imag(f)));
        g1 = max(abs(real(g)), abs(imag(g)));
        rtmax = sqrt(safmax / four_r);
        if f1 > rtmin && f1 < rtmax &&
           g1 > rtmin && g1 < rtmax {
            f2 = _lartg_abssq(f, RT, MT);
            g2 = _lartg_abssq(g, RT, MT);
            h2 = f2 + g2;
            if f2 >= h2*safmin {
                cs.* = sqrt(f2 / h2);
                r.* = f / cs.*;
                rtmax = rtmax * two_r;
                if f2 > rtmin && h2 < rtmax {
                    sn.* = conjugate(g) * (f / sqrt(f2*h2));
                } else {
                    sn.* = conjugate(g) * (r.* / h2);
                }
            } else {
                d = sqrt(f2*h2);
                cs.* = f2 / d;
                if cs.* >= safmin {
                    r.* = f / cs.*;
                } else {
                    r.* = f * (h2 / d);
                }
                sn.* = conjugate(g) * (f / d);
            }
        } else {
            u = min(safmax, max(safmin, f1, g1));
            gs = g / u;
            g2 = _lartg_abssq(gs, RT, MT);
            if f1 / u < rtmin {
                v = min(safmax, max(safmin, f1));
                w = v / u;
                fs = f / v;
                f2 = _lartg_abssq(fs, RT, MT);
                h2 = f2*w*w + g2;
            } else {
                w = one_r;
                fs = f / u;
                f2 = _lartg_abssq(fs, RT, MT);
                h2 = f2 + g2;
            }

            if f2 >= h2*safmin {
                cs.* = sqrt(f2 / h2);
                r.* = fs / cs.*;
                rtmax = rtmax * two_r;
                if f2 > rtmin && h2 < rtmax {
                    sn.* = conjugate(gs) * (fs / sqrt(f2*h2));
                } else {
                    sn.* = conjugate(gs) * (r.* / h2);
                }
            } else {
                d = sqrt(f2*h2);
                cs.* = f2 / d;
                if cs.* >= safmin {
                    r.* = fs / cs.*;
                } else {
                    r.* = fs * (h2 / d);
                }
                sn.* = conjugate(gs) * (fs / d);
            }

            cs.* = cs.* * w;
            r.* = r.* * u;
        }
    }
}
