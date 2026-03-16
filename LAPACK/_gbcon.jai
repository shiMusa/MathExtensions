// ? LAPACK *gbcon()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbcon/cgbcon/dgbcon/zgbcon
// *GBCON estimates the reciprocal of the condition number of a
// general band matrix using the LU factorization computed by *GBTRF.

sgbcon :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    ipiv: $IPIV/ArrayType1D,
    anorm: float32,
    rcond: *float32,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, iwork, info, RT=float32, MT=float32);
}

dgbcon :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    ipiv: $IPIV/ArrayType1D,
    anorm: float64,
    rcond: *float64,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, iwork, info, RT=float64, MT=float64);
}

cgbcon :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    ipiv: $IPIV/ArrayType1D,
    anorm: float32,
    rcond: *float32,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, rwork, info, RT=float32, MT=C32);
}

zgbcon :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    ipiv: $IPIV/ArrayType1D,
    anorm: float64,
    rcond: *float64,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gbcon_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gbcon_lamch_safe_minimum :: () -> $RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gbcon_lacn2 :: (
    $$n: int,
    v: *$V/ArrayType1D,
    x: *$X/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    est: *$RT,
    kase: *int,
    isave: *$ISAVE/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        slacn2(n, v, x, aux, est, kase, isave);
    } else if MT == float64 {
        dlacn2(n, v, x, aux, est, kase, isave);
    } else if MT == C32 {
        clacn2(n, v, x, est, kase, isave);
    } else {
        zlacn2(n, v, x, est, kase, isave);
    }
}

_gbcon_latbs :: (
    $$trans: u8,
    $$normin: u8,
    $$n: int,
    $$kd: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    x: *$X/ArrayType1D,
    scale: *$RT,
    cnorm: *$CNORM/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT && X.data_type == MT && CNORM.data_type == RT;
} {
    #if MT == float32 {
        slatbs(#char "u", trans, #char "n", normin, n, kd, ab, ldab, x, scale, cnorm, info);
    } else if MT == float64 {
        dlatbs(#char "u", trans, #char "n", normin, n, kd, ab, ldab, x, scale, cnorm, info);
    } else if MT == C32 {
        clatbs(#char "u", trans, #char "n", normin, n, kd, ab, ldab, x, scale, cnorm, info);
    } else {
        zlatbs(#char "u", trans, #char "n", normin, n, kd, ab, ldab, x, scale, cnorm, info);
    }
}

_gbcon_rscl :: (
    $$n: int,
    sa: $RT,
    x: *$X/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        srscl(n, sa, x, 1);
    } else if MT == float64 {
        drscl(n, sa, x, 1);
    } else if MT == C32 {
        csrscl(n, sa, x, 1);
    } else {
        zdrscl(n, sa, x, 1);
    }
}

_gbcon_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> int #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
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

_gbcon :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    ipiv: $IPIV/ArrayType1D,
    anorm: $RT,
    rcond: *RT,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if CHECKS {
        #run {
            assert(AB.data_type   == MT);
            assert(IPIV.data_type == int);
            assert(WORK.data_type == MT);
        }
        #if MT == RT {
            #run assert(AUX.data_type == int);
        } else {
            #run assert(AUX.data_type == RT);
        }
    }

    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == float32 {
        "SGBCON"
    } else if MT == float64 {
        "DGBCON"
    } else if MT == C32 {
        "CGBCON"
    } else {
        "ZGBCON"
    };
    trans_char :: #run if MT == RT then #char "t" else #char "c";

    onenrm := norm == #char "1" || char_eq(norm, #char "o");
    lnoti  := kl > 0;
    kd     := kl + ku + 1;
    normin := #char "n";
    kase1  := ifx onenrm then 1 else 2;

    ainvnm, scale : RT = ---;
    t, sum : MT = ---;
    kase, ix, j, jp, latbs_info, lm : int = ---;
    smlnum : RT = ---;

    info.* = 0;
    if !onenrm && !char_eq(norm, #char "i") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if kl < 0 {
        info.* = -3;
    } else if ku < 0 {
        info.* = -4;
    } else if ldab < 2*kl + ku + 1 {
        info.* = -6;
    } else if anorm < zero_r {
        info.* = -8;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    rcond.* = zero_r;
    if n == 0 {
        rcond.* = one_r;
        return;
    } else if anorm == zero_r {
        return;
    }

    smlnum = _gbcon_lamch_safe_minimum(RT);
    ainvnm = zero_r;
    kase = 0;

    x_view := _gbcon_subvector(work.*, 1, n);
    v_view := _gbcon_subvector(work.*, n + 1, n);
    cnorm_view : [] RT = ---;
    #if MT == RT {
        cnorm_view = _gbcon_subvector(work.*, 2*n + 1, n);
    } else {
        cnorm_view = _gbcon_subvector(aux.*, 1, n);
    }

    isave := array1d(3, 0, 0, 0);

    while true {
        _gbcon_lacn2(n, *v_view, *x_view, aux, *ainvnm, *kase, *isave, RT, MT);
        if kase == 0 then break;

        if kase == kase1 {
            if lnoti {
                for j:1..n-1 {
                    lm = min(kl, n-j);
                    jp = fget(ipiv, j);
                    t = fget(x_view, jp);
                    if jp != j {
                        fset(x_view, jp, fget(x_view, j));
                        fset(x_view, j, t);
                    }
                    for i:1..lm {
                        fset(x_view, j+i, fget(x_view, j+i) - t*fget(ab, kd+i, j));
                    }
                }
            }

            latbs_info = 0;
            _gbcon_latbs(#char "n", normin, n, kl+ku, ab, ldab, *x_view, *scale, *cnorm_view, *latbs_info, RT, MT);
        } else {
            latbs_info = 0;
            _gbcon_latbs(trans_char, normin, n, kl+ku, ab, ldab, *x_view, *scale, *cnorm_view, *latbs_info, RT, MT);

            if lnoti {
                for <j:1..n-1 {
                    lm = min(kl, n-j);
                    sum = zero_t;
                    for i:1..lm {
                        sum += conjugate(fget(ab, kd+i, j))*fget(x_view, j+i);
                    }
                    fset(x_view, j, fget(x_view, j) - sum);
                    jp = fget(ipiv, j);
                    if jp != j {
                        t = fget(x_view, jp);
                        fset(x_view, jp, fget(x_view, j));
                        fset(x_view, j, t);
                    }
                }
            }
        }

        normin = #char "y";
        if scale != one_r {
            ix = _gbcon_iamax_1based(n, x_view, RT, MT);
            if scale < _abs1(fget(x_view, ix))*smlnum || scale == zero_r then break;
            _gbcon_rscl(n, scale, *x_view, RT, MT);
        }
    }

    if ainvnm != zero_r {
        rcond.* = (one_r / ainvnm) / anorm;
    }
}
