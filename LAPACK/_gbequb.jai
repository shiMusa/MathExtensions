// ? LAPACK *gbequb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbequb/cgbequb/dgbequb/zgbequb
// *GBEQUB computes row and column scalings intended to equilibrate a
// general band matrix, with scale factors restricted to radix powers.

sgbequb :: #bake_arguments _gbequb(RT=float32, MT=float32);
dgbequb :: #bake_arguments _gbequb(RT=float64, MT=float64);
cgbequb :: #bake_arguments _gbequb(RT=float32, MT=C32);
zgbequb :: #bake_arguments _gbequb(RT=float64, MT=C64);

#scope_file
_gbequb_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gbequb_lamch_radix :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "b");
    } else {
        return dlamch(#char "b");
    }
}

_gbequb_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_gbequb_round_to_radix_power :: inline (x: $RT, radix: RT, logrdx: RT) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    if x > zero(RT) {
        return radix ** cast(int) (log(x) / logrdx);
    }
    return zero(RT);
}

_gbequb :: (
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    r: *$R/ArrayType1D,
    c: *$C/ArrayType1D,
    rowcnd: *RT,
    colcnd: *RT,
    amax: *RT,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if CHECKS {
        #run {
            assert(AB.data_type == MT);
            assert(R.data_type  == RT);
            assert(C.data_type  == RT);
        }
    }

    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == float32 {
        "SGBEQUB"
    } else if MT == float64 {
        "DGBEQUB"
    } else if MT == C32 {
        "CGBEQUB"
    } else {
        "ZGBEQUB"
    };

    bignum, logrdx, radix, rcmax, rcmin, smlnum : RT = ---;
    kd, i, j : int = ---;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if kl < 0 {
        info.* = -3;
    } else if ku < 0 {
        info.* = -4;
    } else if ldab < kl + ku + 1 {
        info.* = -6;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if m == 0 || n == 0 {
        rowcnd.* = one_r;
        colcnd.* = one_r;
        amax.*   = zero_r;
        return;
    }

    smlnum = _gbequb_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;
    radix  = _gbequb_lamch_radix(RT);
    logrdx = log(radix);

    for i:1..m {
        fset(r, i, zero_r);
    }

    kd = ku + 1;
    for j:1..n {
        for i:max(j-ku, 1)..min(j+kl, m) {
            fset(r, i, max(fget(r.*, i), _gbequb_abs1(fget(ab, kd+i-j, j), RT, MT)));
        }
    }

    for i:1..m {
        if fget(r.*, i) > zero_r {
            fset(r, i, _gbequb_round_to_radix_power(fget(r.*, i), radix, logrdx));
        }
    }

    rcmin = bignum;
    rcmax = zero_r;
    for i:1..m {
        rcmax = max(rcmax, fget(r.*, i));
        rcmin = min(rcmin, fget(r.*, i));
    }
    amax.* = rcmax;

    if rcmin == zero_r {
        for i:1..m {
            if fget(r.*, i) == zero_r {
                info.* = i;
                return;
            }
        }
    } else {
        for i:1..m {
            fset(r, i, one_r / min(max(fget(r.*, i), smlnum), bignum));
        }
        rowcnd.* = max(rcmin, smlnum) / min(rcmax, bignum);
    }

    for j:1..n {
        fset(c, j, zero_r);
    }

    kd = ku + 1;
    for j:1..n {
        for i:max(j-ku, 1)..min(j+kl, m) {
            fset(c, j, max(fget(c.*, j), _gbequb_abs1(fget(ab, kd+i-j, j), RT, MT) * fget(r.*, i)));
        }
        if fget(c.*, j) > zero_r {
            fset(c, j, _gbequb_round_to_radix_power(fget(c.*, j), radix, logrdx));
        }
    }

    rcmin = bignum;
    rcmax = zero_r;
    for j:1..n {
        rcmin = min(rcmin, fget(c.*, j));
        rcmax = max(rcmax, fget(c.*, j));
    }

    if rcmin == zero_r {
        for j:1..n {
            if fget(c.*, j) == zero_r {
                info.* = m + j;
                return;
            }
        }
    } else {
        for j:1..n {
            fset(c, j, one_r / min(max(fget(c.*, j), smlnum), bignum));
        }
        colcnd.* = max(rcmin, smlnum) / min(rcmax, bignum);
    }
}
