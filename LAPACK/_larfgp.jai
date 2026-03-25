// ? LAPACK *larfgp()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarfgp/dlarfgp/clarfgp/zlarfgp
// *LARFGP generates a Householder reflector with nonnegative beta.

slarfgp :: #bake_arguments _larfgp_real(RT=float32);
dlarfgp :: #bake_arguments _larfgp_real(RT=float64);
clarfgp :: #bake_arguments _larfgp_complex(RT=float32, MT=C32);
zlarfgp :: #bake_arguments _larfgp_complex(RT=float64, MT=C64);

#scope_file
_larfgp_nrm2 :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 { return snrm2(n, x, incx); }
    #if MT == float64 { return dnrm2(n, x, incx); }
    #if MT == C32 { return scnrm2(n, x, incx); }
    return dznrm2(n, x, incx);
}

_larfgp_lamch :: ($$cmach: u8, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(cmach); } else { return dlamch(cmach); }
}

_larfgp_scal_real :: ($$n: int, $$alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 { sscal(n, alpha, x, incx); } else { dscal(n, alpha, x, incx); }
}

_larfgp_scal_complex_real :: ($$n: int, $$alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_larfgp_scal_complex :: ($$n: int, alpha: $MT, x: *$X/ArrayType1D, $$incx: int, $MT: Type) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { cscal(n, alpha, x, incx); } else { zscal(n, alpha, x, incx); }
}

_larfgp_real :: (
    $$n: int, alpha: *$RT, x: *$X/ArrayType1D, $$incx: int, tau: *RT, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    two_r  :: #run ncast(RT, 2);
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);

    beta, bignum, eps, savealpha, smlnum, xnorm : RT = ---;
    j, knt : int = ---;

    if n <= 0 {
        tau.* = zero_r;
        return;
    }

    eps = _larfgp_lamch(#char "P", RT);
    xnorm = _larfgp_nrm2(n - 1, x.*, incx, RT, RT);

    if xnorm <= eps*abs(alpha.*) {
        if alpha.* >= zero_r {
            tau.* = zero_r;
        } else {
            tau.* = two_r;
            for j:1..n - 1 {
                fset(x, 1 + (j - 1)*incx, zero_r);
            }
            alpha.* = -alpha.*;
        }
    } else {
        beta = sign(if RT == float32 { slapy2(alpha.*, xnorm); } else { dlapy2(alpha.*, xnorm); }, alpha.*);
        smlnum = _larfgp_lamch(#char "S", RT) / _larfgp_lamch(#char "E", RT);
        knt = 0;
        if abs(beta) < smlnum {
            bignum = one_r / smlnum;
            while abs(beta) < smlnum && knt < 20 {
                knt += 1;
                _larfgp_scal_real(n - 1, bignum, x, incx, RT);
                beta = beta*bignum;
                alpha.* = alpha.* * bignum;
            }

            xnorm = _larfgp_nrm2(n - 1, x.*, incx, RT, RT);
            beta = sign(if RT == float32 { slapy2(alpha.*, xnorm); } else { dlapy2(alpha.*, xnorm); }, alpha.*);
        }

        savealpha = alpha.*;
        alpha.* = alpha.* + beta;
        if beta < zero_r {
            beta = -beta;
            tau.* = -alpha.* / beta;
        } else {
            alpha.* = xnorm * (xnorm / alpha.*);
            tau.* = alpha.* / beta;
            alpha.* = -alpha.*;
        }

        if abs(tau.*) <= smlnum {
            if savealpha >= zero_r {
                tau.* = zero_r;
            } else {
                tau.* = two_r;
                for j:1..n - 1 {
                    fset(x, 1 + (j - 1)*incx, zero_r);
                }
                beta = -savealpha;
            }
        } else {
            _larfgp_scal_real(n - 1, one_r / alpha.*, x, incx, RT);
        }

        for j:1..knt {
            beta = beta*smlnum;
        }
        alpha.* = beta;
    }
}

_larfgp_complex :: (
    $$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *MT, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    two_r  :: #run ncast(RT, 2);
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    zero_t :: #run zero(MT);

    alphi, alphr, beta, bignum, eps, smlnum, xnorm : RT = ---;
    savealpha : MT = ---;
    j, knt : int = ---;

    if n <= 0 {
        tau.* = zero_t;
        return;
    }

    eps = _larfgp_lamch(#char "P", RT);
    xnorm = _larfgp_nrm2(n - 1, x.*, incx, RT, MT);
    alphr = real(alpha.*);
    alphi = imag(alpha.*);

    if xnorm <= eps*abs(alpha.*) && alphi == zero_r {
        if alphr >= zero_r {
            tau.* = zero_t;
        } else {
            tau.* = cpx(two_r, zero_r);
            for j:1..n - 1 {
                fset(x, 1 + (j - 1)*incx, zero_t);
            }
            alpha.* = -alpha.*;
        }
    } else {
        #if RT == float32 {
            beta = sign(slapy3(alphr, alphi, xnorm), alphr);
        } else {
            beta = sign(dlapy3(alphr, alphi, xnorm), alphr);
        }

        smlnum = _larfgp_lamch(#char "S", RT) / _larfgp_lamch(#char "E", RT);
        bignum = one_r / smlnum;
        knt = 0;

        if abs(beta) < smlnum {
            while abs(beta) < smlnum && knt < 20 {
                knt += 1;
                _larfgp_scal_complex_real(n - 1, bignum, x, incx, RT, MT);
                beta = beta*bignum;
                alphi = alphi*bignum;
                alphr = alphr*bignum;
            }

            xnorm = _larfgp_nrm2(n - 1, x.*, incx, RT, MT);
            alpha.* = cpx(alphr, alphi);
            beta = sign(if RT == float32 { slapy3(alphr, alphi, xnorm); } else { dlapy3(alphr, alphi, xnorm); }, alphr);
        }

        savealpha = alpha.*;
        alpha.* = alpha.* + beta;
        if beta < zero_r {
            beta = -beta;
            tau.* = -alpha.* / beta;
        } else {
            alphr = alphi * (alphi / real(alpha.*));
            alphr = alphr + xnorm * (xnorm / real(alpha.*));
            tau.* = cpx(alphr / beta, -alphi / beta);
            alpha.* = cpx(-alphr, alphi);
        }
        #if MT == C32 {
            alpha.* = cladiv(cpx(one_r, zero_r), alpha.*);
        } else {
            alpha.* = zladiv(cpx(one_r, zero_r), alpha.*);
        }

        if abs(tau.*) <= smlnum {
            alphr = real(savealpha);
            alphi = imag(savealpha);
            if alphi == zero_r {
                if alphr >= zero_r {
                    tau.* = zero_t;
                } else {
                    tau.* = cpx(two_r, zero_r);
                    for j:1..n - 1 {
                        fset(x, 1 + (j - 1)*incx, zero_t);
                    }
                    beta = real(-savealpha);
                }
            } else {
                #if RT == float32 {
                    xnorm = slapy2(alphr, alphi);
                } else {
                    xnorm = dlapy2(alphr, alphi);
                }
                tau.* = cpx(one_r - alphr / xnorm, -alphi / xnorm);
                for j:1..n - 1 {
                    fset(x, 1 + (j - 1)*incx, zero_t);
                }
                beta = xnorm;
            }
        } else {
            _larfgp_scal_complex(n - 1, alpha.*, x, incx, MT);
        }

        for j:1..knt {
            beta = beta*smlnum;
        }
        alpha.* = cpx(beta, zero_r);
    }
}
