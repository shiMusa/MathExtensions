// ? LAPACK *gesc2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgesc2/cgesc2/dgesc2/zgesc2
// *GESC2 solves a system of linear equations using the LU
// factorization with complete pivoting computed by *GETC2.

sgesc2 :: #bake_arguments _gesc2(RT=float32, MT=float32);
dgesc2 :: #bake_arguments _gesc2(RT=float64, MT=float64);
cgesc2 :: #bake_arguments _gesc2(RT=float32, MT=C32);
zgesc2 :: #bake_arguments _gesc2(RT=float64, MT=C64);

#scope_file
_gesc2_lamch_precision :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "p");
    } else {
        return dlamch(#char "p");
    }
}

_gesc2_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gesc2_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $MT: Type
) -> int #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return isamax(n, x, 1) + 1;
    } else if MT == float64 {
        return idamax(n, x, 1) + 1;
    } else if MT == C32 {
        return icamax(n, x, 1) + 1;
    } else {
        return izamax(n, x, 1) + 1;
    }
}

_gesc2_abs :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    return abs(x);
}

_gesc2_scale_vector :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    for i:1..n {
        fset(x, i, fget(x.*, i) * ncast(MT, alpha));
    }
}

_gesc2 :: (
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    rhs: *$RHS/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    jpiv: $JPIV/ArrayType1D,
    scale: *$RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           RHS.data_type == MT &&
           IPIV.data_type == int &&
           JPIV.data_type == int;
} {
    one_r :: #run one(RT);
    two_r :: #run ncast(RT, 2);

    bignum, eps, smlnum, temp_r : RT = ---;
    temp_t : MT = ---;
    i, j : int = ---;

    eps = _gesc2_lamch_precision(RT);
    smlnum = _gesc2_lamch_safe_minimum(RT) / eps;
    bignum = one_r / smlnum;
    _ = bignum;

    for i:1..n - 1 {
        ip := fget(ipiv, i);
        if ip != i {
            temp_t = fget(rhs.*, i);
            fset(rhs, i, fget(rhs.*, ip));
            fset(rhs, ip, temp_t);
        }
    }

    for i:1..n - 1 {
        for j:i + 1..n {
            fset(rhs, j, fget(rhs.*, j) - fget(a, j, i)*fget(rhs.*, i));
        }
    }

    scale.* = one_r;
    i = _gesc2_iamax_1based(n, rhs.*, MT);
    if two_r*smlnum*_gesc2_abs(fget(rhs.*, i), RT, MT) > _gesc2_abs(fget(a, n, n), RT, MT) {
        temp_r = (one_r / two_r) / _gesc2_abs(fget(rhs.*, i), RT, MT);
        _gesc2_scale_vector(n, temp_r, rhs, RT, MT);
        scale.* *= temp_r;
    }

    for <i:1..n {
        temp_t = ncast(MT, one_r) / fget(a, i, i);
        fset(rhs, i, fget(rhs.*, i) * temp_t);
        for j:i + 1..n {
            fset(rhs, i, fget(rhs.*, i) - fget(rhs.*, j)*(fget(a, i, j)*temp_t));
        }
    }

    for <j:1..n - 1 {
        jp := fget(jpiv, j);
        if jp != j {
            temp_t = fget(rhs.*, j);
            fset(rhs, j, fget(rhs.*, jp));
            fset(rhs, jp, temp_t);
        }
    }
}
