// ? LAPACK *gecon()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgecon/cgecon/dgecon/zgecon
// *GECON estimates the reciprocal of the condition number of a
// general matrix using the LU factorization computed by *GETRF.

sgecon :: (
    $$norm: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    anorm: float32,
    rcond: *float32,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gecon(norm, n, a, lda, anorm, rcond, work, iwork, info, RT=float32, MT=float32);
}

dgecon :: (
    $$norm: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    anorm: float64,
    rcond: *float64,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gecon(norm, n, a, lda, anorm, rcond, work, iwork, info, RT=float64, MT=float64);
}

cgecon :: (
    $$norm: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    anorm: float32,
    rcond: *float32,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gecon(norm, n, a, lda, anorm, rcond, work, rwork, info, RT=float32, MT=C32);
}

zgecon :: (
    $$norm: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    anorm: float64,
    rcond: *float64,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gecon(norm, n, a, lda, anorm, rcond, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gecon_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gecon_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gecon_lamch_overflow :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "o");
    } else {
        return dlamch(#char "o");
    }
}

_gecon_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_gecon_lacn2 :: (
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

_gecon_latrs :: (
    $$uplo: u8,
    $$trans: u8,
    $$diag: u8,
    $$normin: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    x: *$X/ArrayType1D,
    scale: *$RT,
    cnorm: *$CNORM/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           X.data_type == MT &&
           CNORM.data_type == RT;
} {
    #if MT == float32 {
        slatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    } else if MT == float64 {
        dlatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    } else if MT == C32 {
        clatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    } else {
        zlatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    }
}

_gecon_rscl :: (
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

_gecon_iamax_1based :: (
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

_gecon :: (
    $$norm: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    anorm: $RT,
    rcond: *RT,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    routine_name :: #run if MT == float32 {
        "SGECON"
    } else if MT == float64 {
        "DGECON"
    } else if MT == C32 {
        "CGECON"
    } else {
        "ZGECON"
    };
    trans_char :: #run if MT == RT then #char "t" else #char "c";

    onenrm := norm == #char "1" || char_eq(norm, #char "o");
    hugeval := _gecon_lamch_overflow(RT);

    ainvnm, scale, sl, smlnum, su : RT = ---;
    kase1, ix, kase, latrs_info : int = ---;
    normin : u8 = ---;

    info.* = 0;
    if !onenrm && !char_eq(norm, #char "i") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if anorm < zero_r {
        info.* = -5;
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
    } else if _gecon_isnan(anorm, RT) {
        rcond.* = anorm;
        info.* = -5;
        return;
    } else if anorm > hugeval {
        info.* = -5;
        return;
    }

    smlnum = _gecon_lamch_safe_minimum(RT);
    ainvnm = zero_r;
    normin = #char "n";
    kase1 = ifx onenrm then 1 else 2;
    kase = 0;

    x_view := _gecon_subvector(work.*, 1, n);
    v_view := _gecon_subvector(work.*, n + 1, n);
    cnorm_l : [] RT = ---;
    cnorm_u : [] RT = ---;

    #if MT == RT {
        cnorm_l = _gecon_subvector(work.*, 2*n + 1, n);
        cnorm_u = _gecon_subvector(work.*, 3*n + 1, n);
    } else {
        cnorm_l = _gecon_subvector(aux.*, 1, n);
        cnorm_u = _gecon_subvector(aux.*, n + 1, n);
    }

    isave := array1d(3, 0, 0, 0);

    while true {
        _gecon_lacn2(n, *v_view, *x_view, aux, *ainvnm, *kase, *isave, RT, MT);
        if kase == 0 then break;

        if kase == kase1 {
            latrs_info = 0;
            _gecon_latrs(#char "l", #char "n", #char "u", normin, n, a, lda, *x_view, *sl, *cnorm_l, *latrs_info, RT, MT);

            latrs_info = 0;
            _gecon_latrs(#char "u", #char "n", #char "n", normin, n, a, lda, *x_view, *su, *cnorm_u, *latrs_info, RT, MT);
        } else {
            latrs_info = 0;
            _gecon_latrs(#char "u", trans_char, #char "n", normin, n, a, lda, *x_view, *su, *cnorm_u, *latrs_info, RT, MT);

            latrs_info = 0;
            _gecon_latrs(#char "l", trans_char, #char "u", normin, n, a, lda, *x_view, *sl, *cnorm_l, *latrs_info, RT, MT);
        }

        scale = sl * su;
        normin = #char "y";
        if scale != one_r {
            ix = _gecon_iamax_1based(n, x_view, RT, MT);
            if scale < _abs1(fget(x_view, ix))*smlnum || scale == zero_r {
                return;
            }
            _gecon_rscl(n, scale, *x_view, RT, MT);
        }
    }

    if ainvnm != zero_r {
        rcond.* = (one_r / ainvnm) / anorm;
    } else {
        info.* = 1;
        return;
    }

    if _gecon_isnan(rcond.*, RT) || rcond.* > hugeval {
        info.* = 1;
    }
}
