// ? LAPACK *laic1()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaic1/dlaic1/claic1/zlaic1
// *LAIC1 applies one step of incremental condition estimation.

slaic1 :: #bake_arguments _laic1_real(RT=float32);
dlaic1 :: #bake_arguments _laic1_real(RT=float64);
claic1 :: #bake_arguments _laic1_complex(RT=float32, MT=C32);
zlaic1 :: #bake_arguments _laic1_complex(RT=float64, MT=C64);

#scope_file
_laic1_real_dot :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    $RT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           Y.data_type == RT;
} {
    #if RT == float32 {
        return sdot(n, x, incx, y, incy);
    } else {
        return ddot(n, x, incx, y, incy);
    }
}

_laic1_complex_dot :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: $Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) -> MT #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        return cdotc(n, x, incx, y, incy);
    } else {
        return zdotc(n, x, incx, y, incy);
    }
}

_laic1_lamch :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "E");
    } else {
        return dlamch(#char "E");
    }
}

_laic1_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_laic1_complex_norm :: inline (s: $MT, c: MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return sqrt(real(s*conjugate(s) + c*conjugate(c)));
}

_laic1_real :: (
    $$job: int,
    $$j: int,
    x: $X/ArrayType1D,
    sest: $RT,
    w: $W/ArrayType1D,
    gamma: RT,
    sestpr: *RT,
    s: *RT,
    c: *RT,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           W.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2.0);
    half_r :: #run ncast(RT, 0.5);
    four_r :: #run ncast(RT, 4.0);

    absalp, absest, absgam, alpha, cosine, norma, s1, s2, sine, t, test, tmp, zeta1, zeta2 : RT = ---;
    term_b, term_c, eps : RT = ---;

    eps = _laic1_lamch(RT);
    alpha = _laic1_real_dot(j, x, 1, w, 1, RT);

    absalp = abs(alpha);
    absgam = abs(gamma);
    absest = abs(sest);

    if job == 1 {
        if sest == zero_r {
            s1 = max(absgam, absalp);
            if s1 == zero_r {
                s.* = zero_r;
                c.* = one_r;
                sestpr.* = zero_r;
            } else {
                s.* = alpha / s1;
                c.* = gamma / s1;
                tmp = sqrt(s.* * s.* + c.* * c.*);
                s.* = s.* / tmp;
                c.* = c.* / tmp;
                sestpr.* = s1 * tmp;
            }
            return;
        } else if absgam <= eps*absest {
            s.* = one_r;
            c.* = zero_r;
            tmp = max(absest, absalp);
            s1 = absest / tmp;
            s2 = absalp / tmp;
            sestpr.* = tmp * sqrt(s1*s1 + s2*s2);
            return;
        } else if absalp <= eps*absest {
            s1 = absgam;
            s2 = absest;
            if s1 <= s2 {
                s.* = one_r;
                c.* = zero_r;
                sestpr.* = s2;
            } else {
                s.* = zero_r;
                c.* = one_r;
                sestpr.* = s1;
            }
            return;
        } else if absest <= eps*absalp || absest <= eps*absgam {
            s1 = absgam;
            s2 = absalp;
            if s1 <= s2 {
                tmp = s1 / s2;
                s.* = sqrt(one_r + tmp*tmp);
                sestpr.* = s2 * s.*;
                c.* = (gamma / s2) / s.*;
                s.* = sign(one_r, alpha) / s.*;
            } else {
                tmp = s2 / s1;
                c.* = sqrt(one_r + tmp*tmp);
                sestpr.* = s1 * c.*;
                s.* = (alpha / s1) / c.*;
                c.* = sign(one_r, gamma) / c.*;
            }
            return;
        } else {
            zeta1 = alpha / absest;
            zeta2 = gamma / absest;

            term_b = (one_r - zeta1*zeta1 - zeta2*zeta2) * half_r;
            term_c = zeta1 * zeta1;
            if term_b > zero_r {
                t = term_c / (term_b + sqrt(term_b*term_b + term_c));
            } else {
                t = sqrt(term_b*term_b + term_c) - term_b;
            }

            sine = -zeta1 / t;
            cosine = -zeta2 / (one_r + t);
            tmp = sqrt(sine*sine + cosine*cosine);
            s.* = sine / tmp;
            c.* = cosine / tmp;
            sestpr.* = sqrt(t + one_r) * absest;
            return;
        }
    } else if job == 2 {
        if sest == zero_r {
            sestpr.* = zero_r;
            if max(absgam, absalp) == zero_r {
                sine = one_r;
                cosine = zero_r;
            } else {
                sine = -gamma;
                cosine = alpha;
            }
            s1 = max(abs(sine), abs(cosine));
            s.* = sine / s1;
            c.* = cosine / s1;
            tmp = sqrt(s.* * s.* + c.* * c.*);
            s.* = s.* / tmp;
            c.* = c.* / tmp;
            return;
        } else if absgam <= eps*absest {
            s.* = zero_r;
            c.* = one_r;
            sestpr.* = absgam;
            return;
        } else if absalp <= eps*absest {
            s1 = absgam;
            s2 = absest;
            if s1 <= s2 {
                s.* = zero_r;
                c.* = one_r;
                sestpr.* = s1;
            } else {
                s.* = one_r;
                c.* = zero_r;
                sestpr.* = s2;
            }
            return;
        } else if absest <= eps*absalp || absest <= eps*absgam {
            s1 = absgam;
            s2 = absalp;
            if s1 <= s2 {
                tmp = s1 / s2;
                c.* = sqrt(one_r + tmp*tmp);
                sestpr.* = absest * (tmp / c.*);
                s.* = -(gamma / s2) / c.*;
                c.* = sign(one_r, alpha) / c.*;
            } else {
                tmp = s2 / s1;
                s.* = sqrt(one_r + tmp*tmp);
                sestpr.* = absest / s.*;
                c.* = (alpha / s1) / s.*;
                s.* = -sign(one_r, gamma) / s.*;
            }
            return;
        } else {
            zeta1 = alpha / absest;
            zeta2 = gamma / absest;

            norma = max(one_r + zeta1*zeta1 + abs(zeta1*zeta2),
                        abs(zeta1*zeta2) + zeta2*zeta2);

            test = one_r + two_r*(zeta1 - zeta2)*(zeta1 + zeta2);
            if test >= zero_r {
                term_b = (zeta1*zeta1 + zeta2*zeta2 + one_r) * half_r;
                term_c = zeta2 * zeta2;
                t = term_c / (term_b + sqrt(abs(term_b*term_b - term_c)));
                sine = zeta1 / (one_r - t);
                cosine = -zeta2 / t;
                sestpr.* = sqrt(t + four_r*eps*eps*norma) * absest;
            } else {
                term_b = (zeta2*zeta2 + zeta1*zeta1 - one_r) * half_r;
                term_c = zeta1 * zeta1;
                if term_b >= zero_r {
                    t = -term_c / (term_b + sqrt(term_b*term_b + term_c));
                } else {
                    t = term_b - sqrt(term_b*term_b + term_c);
                }
                sine = -zeta1 / t;
                cosine = -zeta2 / (one_r + t);
                sestpr.* = sqrt(one_r + t + four_r*eps*eps*norma) * absest;
            }
            tmp = sqrt(sine*sine + cosine*cosine);
            s.* = sine / tmp;
            c.* = cosine / tmp;
            return;
        }
    }
}

_laic1_complex :: (
    $$job: int,
    $$j: int,
    x: $X/ArrayType1D,
    sest: $RT,
    w: $W/ArrayType1D,
    gamma: $MT,
    sestpr: *RT,
    s: *MT,
    c: *MT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT &&
           W.data_type == MT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2.0);
    half_r :: #run ncast(RT, 0.5);
    four_r :: #run ncast(RT, 4.0);
    zero_t :: #run _laic1_complex_from_real(zero_r, RT, MT);
    one_t  :: #run _laic1_complex_from_real(one_r, RT, MT);

    absalp, absest, absgam, norma, s1, s2, scl, t, test, tmp, zeta1, zeta2 : RT = ---;
    term_b, term_c, eps : RT = ---;
    alpha, cosine, sine : MT = ---;

    eps = _laic1_lamch(RT);
    alpha = _laic1_complex_dot(j, x, 1, w, 1, MT);

    absalp = abs(alpha);
    absgam = abs(gamma);
    absest = abs(sest);

    if job == 1 {
        if sest == zero_r {
            s1 = max(absgam, absalp);
            if s1 == zero_r {
                s.* = zero_t;
                c.* = one_t;
                sestpr.* = zero_r;
            } else {
                s.* = alpha / s1;
                c.* = gamma / s1;
                tmp = _laic1_complex_norm(s.*, c.*, RT, MT);
                s.* = s.* / tmp;
                c.* = c.* / tmp;
                sestpr.* = s1 * tmp;
            }
            return;
        } else if absgam <= eps*absest {
            s.* = one_t;
            c.* = zero_t;
            tmp = max(absest, absalp);
            s1 = absest / tmp;
            s2 = absalp / tmp;
            sestpr.* = tmp * sqrt(s1*s1 + s2*s2);
            return;
        } else if absalp <= eps*absest {
            s1 = absgam;
            s2 = absest;
            if s1 <= s2 {
                s.* = one_t;
                c.* = zero_t;
                sestpr.* = s2;
            } else {
                s.* = zero_t;
                c.* = one_t;
                sestpr.* = s1;
            }
            return;
        } else if absest <= eps*absalp || absest <= eps*absgam {
            s1 = absgam;
            s2 = absalp;
            if s1 <= s2 {
                tmp = s1 / s2;
                scl = sqrt(one_r + tmp*tmp);
                sestpr.* = s2 * scl;
                s.* = (alpha / s2) / scl;
                c.* = (gamma / s2) / scl;
            } else {
                tmp = s2 / s1;
                scl = sqrt(one_r + tmp*tmp);
                sestpr.* = s1 * scl;
                s.* = (alpha / s1) / scl;
                c.* = (gamma / s1) / scl;
            }
            return;
        } else {
            zeta1 = absalp / absest;
            zeta2 = absgam / absest;

            term_b = (one_r - zeta1*zeta1 - zeta2*zeta2) * half_r;
            term_c = zeta1 * zeta1;
            if term_b > zero_r {
                t = term_c / (term_b + sqrt(term_b*term_b + term_c));
            } else {
                t = sqrt(term_b*term_b + term_c) - term_b;
            }

            sine = -(alpha / absest) / t;
            cosine = -(gamma / absest) / (one_r + t);
            tmp = _laic1_complex_norm(sine, cosine, RT, MT);
            s.* = sine / tmp;
            c.* = cosine / tmp;
            sestpr.* = sqrt(t + one_r) * absest;
            return;
        }
    } else if job == 2 {
        if sest == zero_r {
            sestpr.* = zero_r;
            if max(absgam, absalp) == zero_r {
                sine = one_t;
                cosine = zero_t;
            } else {
                sine = -conjugate(gamma);
                cosine = conjugate(alpha);
            }
            s1 = max(abs(sine), abs(cosine));
            s.* = sine / s1;
            c.* = cosine / s1;
            tmp = _laic1_complex_norm(s.*, c.*, RT, MT);
            s.* = s.* / tmp;
            c.* = c.* / tmp;
            return;
        } else if absgam <= eps*absest {
            s.* = zero_t;
            c.* = one_t;
            sestpr.* = absgam;
            return;
        } else if absalp <= eps*absest {
            s1 = absgam;
            s2 = absest;
            if s1 <= s2 {
                s.* = zero_t;
                c.* = one_t;
                sestpr.* = s1;
            } else {
                s.* = one_t;
                c.* = zero_t;
                sestpr.* = s2;
            }
            return;
        } else if absest <= eps*absalp || absest <= eps*absgam {
            s1 = absgam;
            s2 = absalp;
            if s1 <= s2 {
                tmp = s1 / s2;
                scl = sqrt(one_r + tmp*tmp);
                sestpr.* = absest * (tmp / scl);
                s.* = -(conjugate(gamma) / s2) / scl;
                c.* = (conjugate(alpha) / s2) / scl;
            } else {
                tmp = s2 / s1;
                scl = sqrt(one_r + tmp*tmp);
                sestpr.* = absest / scl;
                s.* = -(conjugate(gamma) / s1) / scl;
                c.* = (conjugate(alpha) / s1) / scl;
            }
            return;
        } else {
            zeta1 = absalp / absest;
            zeta2 = absgam / absest;

            norma = max(one_r + zeta1*zeta1 + zeta1*zeta2,
                        zeta1*zeta2 + zeta2*zeta2);

            test = one_r + two_r*(zeta1 - zeta2)*(zeta1 + zeta2);
            if test >= zero_r {
                term_b = (zeta1*zeta1 + zeta2*zeta2 + one_r) * half_r;
                term_c = zeta2 * zeta2;
                t = term_c / (term_b + sqrt(abs(term_b*term_b - term_c)));
                sine = (alpha / absest) / (one_r - t);
                cosine = -(gamma / absest) / t;
                sestpr.* = sqrt(t + four_r*eps*eps*norma) * absest;
            } else {
                term_b = (zeta2*zeta2 + zeta1*zeta1 - one_r) * half_r;
                term_c = zeta1 * zeta1;
                if term_b >= zero_r {
                    t = -term_c / (term_b + sqrt(term_b*term_b + term_c));
                } else {
                    t = term_b - sqrt(term_b*term_b + term_c);
                }
                sine = -(alpha / absest) / t;
                cosine = -(gamma / absest) / (one_r + t);
                sestpr.* = sqrt(one_r + t + four_r*eps*eps*norma) * absest;
            }
            tmp = _laic1_complex_norm(sine, cosine, RT, MT);
            s.* = sine / tmp;
            c.* = cosine / tmp;
            return;
        }
    }
}
