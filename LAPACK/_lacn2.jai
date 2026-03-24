// ? LAPACK *lacn2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slacn2/dlacn2/clacn2/zlacn2
// *LACN2 estimates the 1-norm of a square matrix using reverse communication.

slacn2 :: #bake_arguments _lacn2_real(RT=float32);
dlacn2 :: #bake_arguments _lacn2_real(RT=float64);
clacn2 :: #bake_arguments _lacn2_complex(RT=float32, MT=C32);
zlacn2 :: #bake_arguments _lacn2_complex(RT=float64, MT=C64);

#scope_file
_lacn2_real_asum :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        return sasum(n, x, incx);
    } else {
        return dasum(n, x, incx);
    }
}

_lacn2_real_copy :: inline (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT &&
           Y.data_type == RT;
} {
    #if RT == float32 {
        scopy(n, x, incx, y, incy);
    } else {
        dcopy(n, x, incx, y, incy);
    }
}

_lacn2_real_imax :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           X.data_type == RT;
} {
    #if RT == float32 {
        return isamax(n, x, incx);
    } else {
        return idamax(n, x, incx);
    }
}

_lacn2_complex_asum :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        return scsum1(n, x, incx);
    } else {
        return dzsum1(n, x, incx);
    }
}

_lacn2_complex_copy :: inline (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_lacn2_complex_imax :: inline ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        return icmax1(n, x, incx);
    } else {
        return izmax1(n, x, incx);
    }
}

_lacn2_safe_minimum :: inline ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch("Safe minimum".data);
    } else {
        return dlamch("Safe minimum".data);
    }
}

_lacn2_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_lacn2_real :: (
    $$n: int,
    v: *$V/ArrayType1D,
    x: *$X/ArrayType1D,
    isgn: *$ISGN/ArrayType1D,
    est: *$RT,
    kase: *int,
    isave: *$ISAVE/ArrayType1D,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           V.data_type == RT &&
           X.data_type == RT &&
           ISGN.data_type == int &&
           ISAVE.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);
    itmax :: 5;

    i, jlast : int = ---;
    altsgn, estold, temp, xs : RT = ---;

    if kase.* == 0 {
        for i:1..n {
            fset(x, i, one_r / ncast(RT, n));
        }
        kase.* = 1;
        fset(isave, 1, 1);
        return;
    }

    jump := fget(isave.*, 1);
    if jump == 1 {
        if n == 1 {
            fset(v, 1, fget(x.*, 1));
            est.* = abs(fget(v.*, 1));
            kase.* = 0;
            return;
        }

        est.* = _lacn2_real_asum(n, x.*, 1, RT);
        for i:1..n {
            if fget(x.*, i) >= zero_r {
                xs = one_r;
            } else {
                xs = -one_r;
            }
            fset(x, i, xs);
            fset(isgn, i, cast(int) xs);
        }
        kase.* = 2;
        fset(isave, 1, 2);
        return;
    } else if jump == 2 {
        fset(isave, 2, _lacn2_real_imax(n, x.*, 1, RT));
        fset(isave, 3, 2);
    } else if jump == 3 {
        _lacn2_real_copy(n, x.*, 1, v, 1, RT);
        estold = est.*;
        est.* = _lacn2_real_asum(n, v.*, 1, RT);
        for i:1..n {
            if fget(x.*, i) >= zero_r {
                xs = one_r;
            } else {
                xs = -one_r;
            }
            if cast(int) xs != fget(isgn.*, i) {
                goto after_sign_repeat_check_real;
            }
        }
        goto final_stage_real;
after_sign_repeat_check_real:
        if est.* <= estold {
            goto final_stage_real;
        }

        for i:1..n {
            if fget(x.*, i) >= zero_r {
                xs = one_r;
            } else {
                xs = -one_r;
            }
            fset(x, i, xs);
            fset(isgn, i, cast(int) xs);
        }
        kase.* = 2;
        fset(isave, 1, 4);
        return;
    } else if jump == 4 {
        jlast = fget(isave.*, 2);
        fset(isave, 2, _lacn2_real_imax(n, x.*, 1, RT));
        if fget(x.*, jlast) != abs(fget(x.*, fget(isave.*, 2))) &&
           fget(isave.*, 3) < itmax {
            fset(isave, 3, fget(isave.*, 3) + 1);
            for i:1..n {
                fset(x, i, zero_r);
            }
            fset(x, fget(isave.*, 2), one_r);
            kase.* = 1;
            fset(isave, 1, 3);
            return;
        }

final_stage_real:
        altsgn = one_r;
        for i:1..n {
            fset(x, i, altsgn*(one_r + ncast(RT, i - 1) / ncast(RT, n - 1)));
            altsgn = -altsgn;
        }
        kase.* = 1;
        fset(isave, 1, 5);
        return;
    } else {
        temp = two_r*(_lacn2_real_asum(n, x.*, 1, RT) / ncast(RT, 3*n));
        if temp > est.* {
            _lacn2_real_copy(n, x.*, 1, v, 1, RT);
            est.* = temp;
        }
        kase.* = 0;
        return;
    }

    for i:1..n {
        fset(x, i, zero_r);
    }
    fset(x, fget(isave.*, 2), one_r);
    kase.* = 1;
    fset(isave, 1, 3);
}

_lacn2_complex :: (
    $$n: int,
    v: *$V/ArrayType1D,
    x: *$X/ArrayType1D,
    est: *$RT,
    kase: *int,
    isave: *$ISAVE/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           V.data_type == MT &&
           X.data_type == MT &&
           ISAVE.data_type == int;
} {
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);
    czero  :: #run zero(MT);
    cone   :: #run _lacn2_complex_from_real(one(RT), RT, MT);
    itmax :: 5;

    i, jlast : int = ---;
    absxi, altsgn, estold, safmin, temp : RT = ---;

    safmin = _lacn2_safe_minimum(RT);

    if kase.* == 0 {
        for i:1..n {
            fset(x, i, _lacn2_complex_from_real(one_r / ncast(RT, n), RT, MT));
        }
        kase.* = 1;
        fset(isave, 1, 1);
        return;
    }

    jump := fget(isave.*, 1);
    if jump == 1 {
        if n == 1 {
            fset(v, 1, fget(x.*, 1));
            est.* = abs(fget(v.*, 1));
            kase.* = 0;
            return;
        }

        est.* = _lacn2_complex_asum(n, x.*, 1, RT, MT);
        for i:1..n {
            absxi = abs(fget(x.*, i));
            if absxi > safmin {
                fset(x, i, cpx(real(fget(x.*, i)) / absxi, imag(fget(x.*, i)) / absxi));
            } else {
                fset(x, i, cone);
            }
        }
        kase.* = 2;
        fset(isave, 1, 2);
        return;
    } else if jump == 2 {
        fset(isave, 2, _lacn2_complex_imax(n, x.*, 1, RT, MT));
        fset(isave, 3, 2);
    } else if jump == 3 {
        _lacn2_complex_copy(n, x.*, 1, v, 1, RT, MT);
        estold = est.*;
        est.* = _lacn2_complex_asum(n, v.*, 1, RT, MT);
        if est.* <= estold {
            goto final_stage_complex;
        }

        for i:1..n {
            absxi = abs(fget(x.*, i));
            if absxi > safmin {
                fset(x, i, cpx(real(fget(x.*, i)) / absxi, imag(fget(x.*, i)) / absxi));
            } else {
                fset(x, i, cone);
            }
        }
        kase.* = 2;
        fset(isave, 1, 4);
        return;
    } else if jump == 4 {
        jlast = fget(isave.*, 2);
        fset(isave, 2, _lacn2_complex_imax(n, x.*, 1, RT, MT));
        if abs(fget(x.*, jlast)) != abs(fget(x.*, fget(isave.*, 2))) &&
           fget(isave.*, 3) < itmax {
            fset(isave, 3, fget(isave.*, 3) + 1);
            for i:1..n {
                fset(x, i, czero);
            }
            fset(x, fget(isave.*, 2), cone);
            kase.* = 1;
            fset(isave, 1, 3);
            return;
        }

final_stage_complex:
        altsgn = one_r;
        for i:1..n {
            fset(x, i, _lacn2_complex_from_real(altsgn*(one_r + ncast(RT, i - 1) / ncast(RT, n - 1)), RT, MT));
            altsgn = -altsgn;
        }
        kase.* = 1;
        fset(isave, 1, 5);
        return;
    } else {
        temp = two_r*(_lacn2_complex_asum(n, x.*, 1, RT, MT) / ncast(RT, 3*n));
        if temp > est.* {
            _lacn2_complex_copy(n, x.*, 1, v, 1, RT, MT);
            est.* = temp;
        }
        kase.* = 0;
        return;
    }

    for i:1..n {
        fset(x, i, czero);
    }
    fset(x, fget(isave.*, 2), cone);
    kase.* = 1;
    fset(isave, 1, 3);
}
