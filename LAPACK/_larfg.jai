// ? LAPACK *larfg()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarfg/dlarfg/clarfg/zlarfg
// *LARFG generates a Householder reflector.

slarfg :: #bake_arguments _larfg_real(RT=float32);
dlarfg :: #bake_arguments _larfg_real(RT=float64);
clarfg :: #bake_arguments _larfg_complex(RT=float32, MT=C32);
zlarfg :: #bake_arguments _larfg_complex(RT=float64, MT=C64);

#scope_file
_larfg_nrm2 :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && (MT == RT || MT == complex_with_inner_type(RT)) && X.data_type == MT;
} {
    #if MT == float32 { return snrm2(n, x, incx); }
    #if MT == float64 { return dnrm2(n, x, incx); }
    #if MT == C32 { return scnrm2(n, x, incx); }
    return dznrm2(n, x, incx);
}

_larfg_lamch :: ($$cmach: u8, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(cmach); } else { return dlamch(cmach); }
}

_larfg_scal_real :: ($$n: int, $$alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 { sscal(n, alpha, x, incx); } else { dscal(n, alpha, x, incx); }
}

_larfg_scal_complex_real :: ($$n: int, $$alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_larfg_scal_complex :: ($$n: int, alpha: $MT, x: *$X/ArrayType1D, $$incx: int, $MT: Type) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { cscal(n, alpha, x, incx); } else { zscal(n, alpha, x, incx); }
}

_larfg_real :: (
    $$n: int, alpha: *$RT, x: *$X/ArrayType1D, $$incx: int, tau: *RT, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    j, knt : int = ---;
    beta, rsafmn, safmin, xnorm : RT = ---;

    if n <= 1 {
        tau.* = zero_r;
        return;
    }

    xnorm = _larfg_nrm2(n - 1, x.*, incx, RT, RT);
    if xnorm == zero_r {
        tau.* = zero_r;
    } else {
        beta = -sign(if RT == float32 { slapy2(alpha.*, xnorm); } else { dlapy2(alpha.*, xnorm); }, alpha.*);
        safmin = _larfg_lamch(#char "S", RT) / _larfg_lamch(#char "E", RT);
        knt = 0;
        if abs(beta) < safmin {
            rsafmn = one_r / safmin;
            while abs(beta) < safmin && knt < 20 {
                knt = knt + 1;
                _larfg_scal_real(n - 1, rsafmn, x, incx, RT);
                beta = beta*rsafmn;
                alpha.* = alpha.* * rsafmn;
            }
            xnorm = _larfg_nrm2(n - 1, x.*, incx, RT, RT);
            beta = -sign(if RT == float32 { slapy2(alpha.*, xnorm); } else { dlapy2(alpha.*, xnorm); }, alpha.*);
        }
        tau.* = (beta - alpha.*) / beta;
        _larfg_scal_real(n - 1, one_r / (alpha.* - beta), x, incx, RT);
        for j:1..knt {
            beta = beta*safmin;
        }
        alpha.* = beta;
    }
}

_larfg_complex :: (
    $$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *MT, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    one_t  :: #run cpx(one_r, zero_r);
    zero_t :: #run cpx(zero_r, zero_r);

    j, knt : int = ---;
    alphi, alphr, beta, rsafmn, safmin, xnorm : RT = ---;

    if n <= 0 {
        tau.* = zero_t;
        return;
    }

    xnorm = _larfg_nrm2(n - 1, x.*, incx, RT, MT);
    alphr = real(alpha.*);
    alphi = imag(alpha.*);
    if xnorm == zero_r && alphi == zero_r {
        tau.* = zero_t;
    } else {
        #if RT == float32 {
            beta = -sign(slapy3(alphr, alphi, xnorm), alphr);
        } else {
            beta = -sign(dlapy3(alphr, alphi, xnorm), alphr);
        }
        safmin = _larfg_lamch(#char "S", RT) / _larfg_lamch(#char "E", RT);
        rsafmn = one_r / safmin;
        knt = 0;
        if abs(beta) < safmin {
            while abs(beta) < safmin && knt < 20 {
                knt = knt + 1;
                _larfg_scal_complex_real(n - 1, rsafmn, x, incx, RT, MT);
                beta = beta*rsafmn;
                alphi = alphi*rsafmn;
                alphr = alphr*rsafmn;
            }
            xnorm = _larfg_nrm2(n - 1, x.*, incx, RT, MT);
            beta = -sign(if RT == float32 { slapy3(alphr, alphi, xnorm); } else { dlapy3(alphr, alphi, xnorm); }, alphr);
        }
        tau.* = cpx((beta - alphr) / beta, -alphi / beta);
        alpha.* = one_t / (alpha.* - beta);
        _larfg_scal_complex(n - 1, alpha.*, x, incx, MT);
        for j:1..knt {
            beta = beta*safmin;
        }
        alpha.* = cpx(beta, zero_r);
    }
}
