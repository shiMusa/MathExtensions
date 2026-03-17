// ? LAPACK *geequ()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgeequ/cgeequ/dgeequ/zgeequ
// *GEEQU computes row and column scalings intended to equilibrate a
// general matrix.

sgeequ :: #bake_arguments _geequ(RT=float32, MT=float32);
dgeequ :: #bake_arguments _geequ(RT=float64, MT=float64);
cgeequ :: #bake_arguments _geequ(RT=float32, MT=C32);
zgeequ :: #bake_arguments _geequ(RT=float64, MT=C64);

#scope_file
_geequ_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_geequ_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_geequ :: (
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
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
            assert(A.data_type == MT);
            assert(R.data_type == RT);
            assert(C.data_type == RT);
        }
    }

    one_r  :: #run one(RT);
    zero_r :: #run zero(RT);
    routine_name :: #run if MT == float32 {
        "SGEEQU"
    } else if MT == float64 {
        "DGEEQU"
    } else if MT == C32 {
        "CGEEQU"
    } else {
        "ZGEEQU"
    };

    bignum, rcmax, rcmin, smlnum : RT = ---;
    i, j : int = ---;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
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

    smlnum = _geequ_lamch_safe_minimum(RT);
    bignum = one_r / smlnum;

    for i:1..m {
        fset(r, i, zero_r);
    }

    for j:1..n {
        for i:1..m {
            fset(r, i, max(fget(r.*, i), _geequ_abs1(fget(a, i, j), RT, MT)));
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

    for j:1..n {
        for i:1..m {
            fset(c, j, max(fget(c.*, j), _geequ_abs1(fget(a, i, j), RT, MT) * fget(r.*, i)));
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
