// ? LAPACK *lacon()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slacon/dlacon/clacon/zlacon
// *LACON estimates the 1-norm of a square matrix using reverse communication.

slacon :: #bake_arguments _lacon_real(RT=float32);
dlacon :: #bake_arguments _lacon_real(RT=float64);
clacon :: #bake_arguments _lacon_complex(RT=float32, MT=C32);
zlacon :: #bake_arguments _lacon_complex(RT=float64, MT=C64);

#scope_file
_slacon_jump: int;
_slacon_j: int;
_slacon_iter: int;

_dlacon_jump: int;
_dlacon_j: int;
_dlacon_iter: int;

_clacon_jump: int;
_clacon_j: int;
_clacon_iter: int;

_zlacon_jump: int;
_zlacon_j: int;
_zlacon_iter: int;

_lacon_real_get_jump :: inline ($RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return _slacon_jump; } else { return _dlacon_jump; }
}

_lacon_real_set_jump :: inline (value: int, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { _slacon_jump = value; } else { _dlacon_jump = value; }
}

_lacon_real_get_j :: inline ($RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return _slacon_j; } else { return _dlacon_j; }
}

_lacon_real_set_j :: inline (value: int, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { _slacon_j = value; } else { _dlacon_j = value; }
}

_lacon_real_get_iter :: inline ($RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return _slacon_iter; } else { return _dlacon_iter; }
}

_lacon_real_set_iter :: inline (value: int, $RT: Type) #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { _slacon_iter = value; } else { _dlacon_iter = value; }
}

_lacon_complex_get_jump :: inline ($RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 { return _clacon_jump; } else { return _zlacon_jump; }
}

_lacon_complex_set_jump :: inline (value: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 { _clacon_jump = value; } else { _zlacon_jump = value; }
}

_lacon_complex_get_j :: inline ($RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 { return _clacon_j; } else { return _zlacon_j; }
}

_lacon_complex_set_j :: inline (value: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 { _clacon_j = value; } else { _zlacon_j = value; }
}

_lacon_complex_get_iter :: inline ($RT: Type, $MT: Type) -> int #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 { return _clacon_iter; } else { return _zlacon_iter; }
}

_lacon_complex_set_iter :: inline (value: int, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    #if MT == C32 { _clacon_iter = value; } else { _zlacon_iter = value; }
}

_lacon_real :: (
    $$n: int,
    v: *$V/ArrayType1D,
    x: *$X/ArrayType1D,
    isgn: *$ISGN/ArrayType1D,
    est: *$RT,
    kase: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           V.data_type == RT &&
           X.data_type == RT &&
           ISGN.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);
    itmax :: 5;

    i, jlast, jump, jstate, iter_state : int = ---;
    altsgn, estold, temp, xs : RT = ---;

    jump = _lacon_real_get_jump(RT);
    jstate = _lacon_real_get_j(RT);
    iter_state = _lacon_real_get_iter(RT);

    if kase.* == 0 {
        for i:1..n {
            fset(x, i, one_r / ncast(RT, n));
        }
        kase.* = 1;
        _lacon_real_set_jump(1, RT);
        return;
    }

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
        _lacon_real_set_jump(2, RT);
        return;
    } else if jump == 2 {
        jstate = _lacn2_real_imax(n, x.*, 1, RT);
        iter_state = 2;
        _lacon_real_set_j(jstate, RT);
        _lacon_real_set_iter(iter_state, RT);
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
                goto after_sign_repeat_check_lacon_real;
            }
        }
        goto final_stage_lacon_real;
after_sign_repeat_check_lacon_real:
        if est.* <= estold {
            goto final_stage_lacon_real;
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
        _lacon_real_set_jump(4, RT);
        return;
    } else if jump == 4 {
        jlast = jstate;
        jstate = _lacn2_real_imax(n, x.*, 1, RT);
        _lacon_real_set_j(jstate, RT);
        if fget(x.*, jlast) != abs(fget(x.*, jstate)) &&
           iter_state < itmax {
            iter_state = iter_state + 1;
            _lacon_real_set_iter(iter_state, RT);
            for i:1..n {
                fset(x, i, zero_r);
            }
            fset(x, jstate, one_r);
            kase.* = 1;
            _lacon_real_set_jump(3, RT);
            return;
        }

final_stage_lacon_real:
        altsgn = one_r;
        for i:1..n {
            fset(x, i, altsgn*(one_r + ncast(RT, i - 1) / ncast(RT, n - 1)));
            altsgn = -altsgn;
        }
        kase.* = 1;
        _lacon_real_set_jump(5, RT);
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
    fset(x, jstate, one_r);
    kase.* = 1;
    _lacon_real_set_jump(3, RT);
}

_lacon_complex :: (
    $$n: int,
    v: *$V/ArrayType1D,
    x: *$X/ArrayType1D,
    est: *$RT,
    kase: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           V.data_type == MT &&
           X.data_type == MT;
} {
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);
    czero  :: #run zero(MT);
    cone   :: #run _lacn2_complex_from_real(one(RT), RT, MT);
    itmax :: 5;

    i, jlast, jump, jstate, iter_state : int = ---;
    absxi, altsgn, estold, safmin, temp : RT = ---;

    safmin = _lacn2_safe_minimum(RT);
    jump = _lacon_complex_get_jump(RT, MT);
    jstate = _lacon_complex_get_j(RT, MT);
    iter_state = _lacon_complex_get_iter(RT, MT);

    if kase.* == 0 {
        for i:1..n {
            fset(x, i, _lacn2_complex_from_real(one_r / ncast(RT, n), RT, MT));
        }
        kase.* = 1;
        _lacon_complex_set_jump(1, RT, MT);
        return;
    }

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
        _lacon_complex_set_jump(2, RT, MT);
        return;
    } else if jump == 2 {
        jstate = _lacn2_complex_imax(n, x.*, 1, RT, MT);
        iter_state = 2;
        _lacon_complex_set_j(jstate, RT, MT);
        _lacon_complex_set_iter(iter_state, RT, MT);
    } else if jump == 3 {
        _lacn2_complex_copy(n, x.*, 1, v, 1, RT, MT);
        estold = est.*;
        est.* = _lacn2_complex_asum(n, v.*, 1, RT, MT);
        if est.* <= estold {
            goto final_stage_lacon_complex;
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
        _lacon_complex_set_jump(4, RT, MT);
        return;
    } else if jump == 4 {
        jlast = jstate;
        jstate = _lacn2_complex_imax(n, x.*, 1, RT, MT);
        _lacon_complex_set_j(jstate, RT, MT);
        if abs(fget(x.*, jlast)) != abs(fget(x.*, jstate)) &&
           iter_state < itmax {
            iter_state = iter_state + 1;
            _lacon_complex_set_iter(iter_state, RT, MT);
            for i:1..n {
                fset(x, i, czero);
            }
            fset(x, jstate, cone);
            kase.* = 1;
            _lacon_complex_set_jump(3, RT, MT);
            return;
        }

final_stage_lacon_complex:
        altsgn = one_r;
        for i:1..n {
            fset(x, i, _lacn2_complex_from_real(altsgn*(one_r + ncast(RT, i - 1) / ncast(RT, n - 1)), RT, MT));
            altsgn = -altsgn;
        }
        kase.* = 1;
        _lacon_complex_set_jump(5, RT, MT);
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
    fset(x, jstate, cone);
    kase.* = 1;
    _lacon_complex_set_jump(3, RT, MT);
}
