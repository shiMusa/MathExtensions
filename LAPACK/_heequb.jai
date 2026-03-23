// ? LAPACK *heequb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cheequb/zheequb
// *HEEQUB computes Hermitian equilibration factors using the BIN
// algorithm, with scale factors restricted to radix powers.

cheequb :: #bake_arguments _heequb(RT=float32, MT=C32);
zheequb :: #bake_arguments _heequb(RT=float64, MT=C64);

#scope_file
_heequb_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_heequb_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_heequb_lamch_safemin :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "S");
    } else {
        return dlamch(#char "S");
    }
}

_heequb_lamch_base :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "B");
    } else {
        return dlamch(#char "B");
    }
}

_heequb_lassq :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    scale: *$RT,
    sumsq: *RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        classq(n, x, incx, scale, sumsq);
    } else {
        zlassq(n, x, incx, scale, sumsq);
    }
}

_heequb :: (
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    s: *$S/ArrayType1D,
    scond: *RT,
    amax: *RT,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           S.data_type == RT &&
           WORK.data_type == MT;
} {
    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    zero_t :: #run zero(MT);
    max_iter :: 100;
    routine_name :: #run if MT == C32 { "CHEEQUB" } else { "ZHEEQUB" };

    up : bool = ---;
    avg, std, tol, c0, c1, c2, t, u, si, dval, base, smin, smax, smlnum,
    bignum, scale, sumsq : RT = ---;
    i, j, iter : int = ---;

    info.* = 0;
    if !(lsame(uplo, #char "U") || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    up = lsame(uplo, #char "U");
    amax.* = zero_r;

    if n == 0 {
        scond.* = one_r;
        return;
    }

    for i:1..n {
        fset(s, i, zero_r);
    }

    if up {
        for j:1..n {
            for i:1..j - 1 {
                t = _heequb_abs1(fget(a, i, j), RT, MT);
                fset(s, i, max(fget(s.*, i), t));
                fset(s, j, max(fget(s.*, j), t));
                amax.* = max(amax.*, t);
            }
            t = _heequb_abs1(fget(a, j, j), RT, MT);
            fset(s, j, max(fget(s.*, j), t));
            amax.* = max(amax.*, t);
        }
    } else {
        for j:1..n {
            t = _heequb_abs1(fget(a, j, j), RT, MT);
            fset(s, j, max(fget(s.*, j), t));
            amax.* = max(amax.*, t);
            for i:j + 1..n {
                t = _heequb_abs1(fget(a, i, j), RT, MT);
                fset(s, i, max(fget(s.*, i), t));
                fset(s, j, max(fget(s.*, j), t));
                amax.* = max(amax.*, t);
            }
        }
    }

    for j:1..n {
        fset(s, j, one_r / fget(s.*, j));
    }

    tol = one_r / sqrt(cast(RT) (2*n));

    for iter:1..max_iter {
        scale = zero_r;
        sumsq = zero_r;

        for i:1..n {
            fset(work, i, zero_t);
        }

        if up {
            for j:1..n {
                for i:1..j - 1 {
                    t = _heequb_abs1(fget(a, i, j), RT, MT);
                    fset(work, i, fget(work.*, i) + _heequb_complex_from_real(t * fget(s.*, j), RT, MT));
                    fset(work, j, fget(work.*, j) + _heequb_complex_from_real(t * fget(s.*, i), RT, MT));
                }
                t = _heequb_abs1(fget(a, j, j), RT, MT);
                fset(work, j, fget(work.*, j) + _heequb_complex_from_real(t * fget(s.*, j), RT, MT));
            }
        } else {
            for j:1..n {
                t = _heequb_abs1(fget(a, j, j), RT, MT);
                fset(work, j, fget(work.*, j) + _heequb_complex_from_real(t * fget(s.*, j), RT, MT));
                for i:j + 1..n {
                    t = _heequb_abs1(fget(a, i, j), RT, MT);
                    fset(work, i, fget(work.*, i) + _heequb_complex_from_real(t * fget(s.*, j), RT, MT));
                    fset(work, j, fget(work.*, j) + _heequb_complex_from_real(t * fget(s.*, i), RT, MT));
                }
            }
        }

        avg = zero_r;
        for i:1..n {
            avg += fget(s.*, i) * real(fget(work.*, i));
        }
        avg /= cast(RT) n;

        for i:n + 1..2*n {
            fset(work, i, _heequb_complex_from_real(fget(s.*, i - n) *
                                                    real(fget(work.*, i - n)) - avg,
                                                    RT, MT));
        }

        work_tail : [] MT;
        work_tail.data = work.data + n;
        work_tail.count = n;
        _heequb_lassq(n, *work_tail, 1, *scale, *sumsq, RT, MT);
        std = scale * sqrt(sumsq / cast(RT) n);

        if std < tol * avg then break;

        for i:1..n {
            t = _heequb_abs1(fget(a, i, i), RT, MT);
            si = fget(s.*, i);
            c2 = cast(RT) (n - 1) * t;
            c1 = cast(RT) (n - 2) * (real(fget(work.*, i)) - t*si);
            c0 = -(t*si)*si + cast(RT) 2 * real(fget(work.*, i)) * si - cast(RT) n * avg;
            dval = c1*c1 - cast(RT) 4 * c0*c2;

            if dval <= zero_r {
                info.* = -1;
                return;
            }
            si = -cast(RT) 2 * c0 / (c1 + sqrt(dval));

            dval = si - fget(s.*, i);
            u = zero_r;
            if up {
                for j:1..i {
                    t = _heequb_abs1(fget(a, j, i), RT, MT);
                    u += fget(s.*, j) * t;
                    fset(work, j, fget(work.*, j) + _heequb_complex_from_real(dval * t, RT, MT));
                }
                for j:i + 1..n {
                    t = _heequb_abs1(fget(a, i, j), RT, MT);
                    u += fget(s.*, j) * t;
                    fset(work, j, fget(work.*, j) + _heequb_complex_from_real(dval * t, RT, MT));
                }
            } else {
                for j:1..i {
                    t = _heequb_abs1(fget(a, i, j), RT, MT);
                    u += fget(s.*, j) * t;
                    fset(work, j, fget(work.*, j) + _heequb_complex_from_real(dval * t, RT, MT));
                }
                for j:i + 1..n {
                    t = _heequb_abs1(fget(a, j, i), RT, MT);
                    u += fget(s.*, j) * t;
                    fset(work, j, fget(work.*, j) + _heequb_complex_from_real(dval * t, RT, MT));
                }
            }

            avg += (u + real(fget(work.*, i))) * dval / cast(RT) n;
            fset(s, i, si);
        }
    }

    smlnum = _heequb_lamch_safemin(RT);
    bignum = one_r / smlnum;
    smin = bignum;
    smax = zero_r;
    t = one_r / sqrt(avg);
    base = _heequb_lamch_base(RT);
    u = one_r / log(base);
    for i:1..n {
        fset(s, i, base ** cast(int) (u * log(fget(s.*, i) * t)));
        smin = min(smin, fget(s.*, i));
        smax = max(smax, fget(s.*, i));
    }
    scond.* = max(smin, smlnum) / min(smax, bignum);
}
