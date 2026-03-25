// ? LAPACK *largv()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slargv/dlargv/clargv/zlargv
// *LARGV generates a vector of plane rotations.

slargv :: #bake_arguments _largv_real(RT=float32);
dlargv :: #bake_arguments _largv_real(RT=float64);
clargv :: #bake_arguments _largv_complex(RT=float32, MT=C32);
zlargv :: #bake_arguments _largv_complex(RT=float64, MT=C64);

#scope_file
_largv_lamch :: ($$cmach: u8, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(cmach); } else { return dlamch(cmach); }
}

_largv_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return max(abs(real(z)), abs(imag(z)));
}

_largv_abssq :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return real(z)*real(z) + imag(z)*imag(z);
}

_largv_lapy2 :: inline ($$x: $RT, $$y: RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slapy2(x, y); } else { return dlapy2(x, y); }
}

_largv_real :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int,
    c: *$C/ArrayType1D, $$incc: int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT && Y.data_type == RT && C.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);

    f, g, t, tt : RT = ---;
    i, ic, ix, iy : int = ---;

    ix = 1;
    iy = 1;
    ic = 1;
    for i:1..n {
        f = fget(x.*, ix);
        g = fget(y.*, iy);
        if g == zero_r {
            fset(c, ic, one_r);
        } else if f == zero_r {
            fset(c, ic, zero_r);
            fset(y, iy, one_r);
            fset(x, ix, g);
        } else if abs(f) > abs(g) {
            t = g / f;
            tt = sqrt(one_r + t*t);
            fset(c, ic, one_r / tt);
            fset(y, iy, t*fget(c.*, ic));
            fset(x, ix, f*tt);
        } else {
            t = f / g;
            tt = sqrt(one_r + t*t);
            fset(y, iy, one_r / tt);
            fset(c, ic, t*fget(y.*, iy));
            fset(x, ix, g*tt);
        }

        ic += incc;
        iy += incy;
        ix += incx;
    }
}

_largv_complex :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int,
    c: *$C/ArrayType1D, $$incc: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT && Y.data_type == MT && C.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);
    zero_t :: #run zero(MT);

    count, i, ic, ix, iy, j : int = ---;
    cs, d, di, dr, eps, f2, f2s, g2, g2s, safmin, safmn2, safmx2, scale : RT = ---;
    f, ff, fs, g, gs, r, sn : MT = ---;

    safmin = _largv_lamch(#char "S", RT);
    eps = _largv_lamch(#char "E", RT);
    safmn2 = _largv_lamch(#char "B", RT) ** cast(int) (log(safmin / eps) / log(_largv_lamch(#char "B", RT)) / two_r);
    safmx2 = one_r / safmn2;

    ix = 1;
    iy = 1;
    ic = 1;
    for i:1..n {
        f = fget(x.*, ix);
        g = fget(y.*, iy);

        scale = max(_largv_abs1(f, RT, MT), _largv_abs1(g, RT, MT));
        fs = f;
        gs = g;
        count = 0;

        if scale >= safmx2 {
            while scale >= safmx2 && count < 20 {
                count += 1;
                fs = fs*safmn2;
                gs = gs*safmn2;
                scale = scale*safmn2;
            }
        } else if scale <= safmn2 {
            if g == zero_t {
                cs = one_r;
                sn = zero_t;
                r = f;
                fset(c, ic, cs);
                fset(y, iy, sn);
                fset(x, ix, r);
                ic += incc;
                iy += incy;
                ix += incx;
                continue;
            }
            while scale <= safmn2 {
                count -= 1;
                fs = fs*safmx2;
                gs = gs*safmx2;
                scale = scale*safmx2;
            }
        }

        f2 = _largv_abssq(fs, RT, MT);
        g2 = _largv_abssq(gs, RT, MT);
        if f2 <= max(g2, one_r)*safmin {
            if f == zero_t {
                cs = zero_r;
                r = cpx(_largv_lapy2(real(g), imag(g), RT), zero_r);
                d = _largv_lapy2(real(gs), imag(gs), RT);
                sn = cpx(real(gs) / d, -imag(gs) / d);
            } else {
                f2s = _largv_lapy2(real(fs), imag(fs), RT);
                g2s = sqrt(g2);
                cs = f2s / g2s;
                if _largv_abs1(f, RT, MT) > one_r {
                    d = _largv_lapy2(real(f), imag(f), RT);
                    ff = cpx(real(f) / d, imag(f) / d);
                } else {
                    dr = safmx2*real(f);
                    di = safmx2*imag(f);
                    d = _largv_lapy2(dr, di, RT);
                    ff = cpx(dr / d, di / d);
                }
                sn = ff * cpx(real(gs) / g2s, -imag(gs) / g2s);
                r = cs*f + sn*g;
            }
        } else {
            f2s = sqrt(one_r + g2 / f2);
            r = cpx(f2s*real(fs), f2s*imag(fs));
            cs = one_r / f2s;
            d = f2 + g2;
            sn = cpx(real(r) / d, imag(r) / d);
            sn = sn*conjugate(gs);
            if count != 0 {
                if count > 0 {
                    for j:1..count {
                        r = r*safmx2;
                    }
                } else {
                    for j:1..-count {
                        r = r*safmn2;
                    }
                }
            }
        }

        fset(c, ic, cs);
        fset(y, iy, sn);
        fset(x, ix, r);
        ic += incc;
        iy += incy;
        ix += incx;
    }
}
