// ? LAPACK *gbrfsx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbrfsx/cgbrfsx/dgbrfsx/zgbrfsx
// *GBRFSX improves the computed solution to a system of linear
// equations and provides error bounds and backward error estimates
// for the solution.

sgbrfsx :: (
    $$trans: u8,
    $$equed: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    r: $R/ArrayType1D,
    c: $C/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *float32,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info, RT=float32, MT=float32);
}

dgbrfsx :: (
    $$trans: u8,
    $$equed: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    r: $R/ArrayType1D,
    c: $C/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *float64,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info, RT=float64, MT=float64);
}

cgbrfsx :: (
    $$trans: u8,
    $$equed: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    r: $R/ArrayType1D,
    c: $C/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *float32,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float32, MT=C32);
}

zgbrfsx :: (
    $$trans: u8,
    $$equed: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    r: $R/ArrayType1D,
    c: $C/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *float64,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gbrfsx_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_gbrfsx_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gbrfsx_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
    #if CHECKS {
        assert(ld >= max(1, count));
        assert(col >= 1);
        assert(count >= 0);
    }

    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld;
    view.count = count;
    return view;
}

_gbrfsx_complex_scratch :: inline (x: $X/ArrayType1D, $$count: int, $RT: Type, $MT: Type) -> [] MT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == RT;
} {
    #if CHECKS {
        assert(x.count >= 2*count);
        #run assert(size_of(MT) == 2*size_of(RT));
    }

    view: [] MT;
    view.data = cast(*MT) x.data;
    view.count = count;
    return view;
}

_gbrfsx_ilaprec :: ($RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return ilaprec(#char "d");
    } else {
        return ilaprec(#char "e");
    }
}

_gbrfsx_langb :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    aux: *$AUX/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           ((MT == RT && AUX.data_type == RT) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        return slangb(norm, n, kl, ku, ab, ldab, aux);
    } else if MT == float64 {
        return dlangb(norm, n, kl, ku, ab, ldab, aux);
    } else if MT == C32 {
        return clangb(norm, n, kl, ku, ab, ldab, aux);
    } else {
        return zlangb(norm, n, kl, ku, ab, ldab, aux);
    }
}

_gbrfsx_gbcon :: (
    $$norm: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
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
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        sgbcon(norm, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, aux, info);
    } else if MT == float64 {
        dgbcon(norm, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, aux, info);
    } else if MT == C32 {
        cgbcon(norm, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, aux, info);
    } else {
        zgbcon(norm, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, aux, info);
    }
}

_gbrfsx_extended :: (
    $$prec_type: int,
    $$trans_type: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    $$colequ: bool,
    c: $C/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    y: *$Y/ArrayType2D,
    $$ldy: int,
    berr_out: *$BERR_OUT/ArrayType1D,
    $$n_norms: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    rcond: *$RT,
    $$ithresh: int,
    rthresh: RT,
    unstable_thresh: RT,
    $$ignore_cwise: bool,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           C.data_type == RT &&
           B.data_type == MT &&
           Y.data_type == MT &&
           BERR_OUT.data_type == RT &&
           ERR_BNDS_NORM.data_type == RT &&
           ERR_BNDS_COMP.data_type == RT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        res_view := _gbrfsx_subvector(work.*, n + 1, n);
        ayb_view := _gbrfsx_subvector(work.*, 1, n);
        dy_view  := _gbrfsx_subvector(work.*, 2*n + 1, n);
        y_tail_view := _gbrfsx_subvector(work.*, 1, n);
        sla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, *res_view, *ayb_view, *dy_view, *y_tail_view, *rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info);
    } else if MT == float64 {
        res_view := _gbrfsx_subvector(work.*, n + 1, n);
        ayb_view := _gbrfsx_subvector(work.*, 1, n);
        dy_view  := _gbrfsx_subvector(work.*, 2*n + 1, n);
        y_tail_view := _gbrfsx_subvector(work.*, 1, n);
        dla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, *res_view, *ayb_view, *dy_view, *y_tail_view, *rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info);
    } else if MT == C32 {
        dy_view := _gbrfsx_subvector(work.*, n + 1, n);
        y_tail_view := _gbrfsx_complex_scratch(aux.*, n, RT, MT);
        cla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, work, aux, *dy_view, *y_tail_view, *rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info);
    } else {
        dy_view := _gbrfsx_subvector(work.*, n + 1, n);
        y_tail_view := _gbrfsx_complex_scratch(aux.*, n, RT, MT);
        zla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, work, aux, *dy_view, *y_tail_view, *rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info);
    }
}

_gbrfsx_gbrcond_norm :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    $$cmode: int,
    c: $C/ArrayType1D,
    info: *int,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           C.data_type == RT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        return sla_gbrcond(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, cmode, c, info, work, aux);
    } else if MT == float64 {
        return dla_gbrcond(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, cmode, c, info, work, aux);
    } else if MT == C32 {
        return cla_gbrcond_c(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, cmode == -1, info, work, aux);
    } else {
        return zla_gbrcond_c(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, cmode == -1, info, work, aux);
    }
}

_gbrfsx_gbrcond_x :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    x: $X/ArrayType1D,
    info: *int,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           X.data_type == MT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    #if MT == float32 {
        return sla_gbrcond(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, 1, x, info, work, aux);
    } else if MT == float64 {
        return dla_gbrcond(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, 1, x, info, work, aux);
    } else if MT == C32 {
        return cla_gbrcond_x(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, x, info, work, aux);
    } else {
        return zla_gbrcond_x(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, x, info, work, aux);
    }
}

_gbrfsx :: (
    $$trans: u8,
    $$equed: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    r: $R/ArrayType1D,
    c: $C/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    rcond: *RT,
    berr: *$BERR/ArrayType1D,
    $$n_err_bnds: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    $$nparams: int,
    params: *$PARAMS/ArrayType1D,
    work: *$WORK/ArrayType1D,
    aux: *$AUX/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT &&
           AFB.data_type == MT &&
           IPIV.data_type == int &&
           R.data_type == RT &&
           C.data_type == RT &&
           B.data_type == MT &&
           X.data_type == MT &&
           BERR.data_type == RT &&
           ERR_BNDS_NORM.data_type == RT &&
           ERR_BNDS_COMP.data_type == RT &&
           PARAMS.data_type == RT &&
           WORK.data_type == MT &&
           ((MT == RT && AUX.data_type == int) || (MT == complex_with_inner_type(RT) && AUX.data_type == RT));
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    half_r :: #run ncast(RT, 0.5);
    quarter_r :: #run ncast(RT, 0.25);
    ten_r :: #run ncast(RT, 10);
    routine_name :: #run if MT == float32 {
        "SGBRFSX"
    } else if MT == float64 {
        "DGBRFSX"
    } else if MT == C32 {
        "CGBRFSX"
    } else {
        "ZGBRFSX"
    };

    la_linrx_itref_i :: 1;
    la_linrx_ithresh_i :: 2;
    la_linrx_cwise_i :: 3;
    la_linrx_trust_i :: 1;
    la_linrx_err_i :: 2;
    la_linrx_rcond_i :: 3;

    anorm, cwise_wrong, err_lbnd, illrcond_thresh, rcond_tmp, rthresh, unstable_thresh : RT = ---;
    cwise_param : RT = ---;
    ithresh, j, n_norms, prec_type, ref_type, trans_type : int = ---;
    colequ, ignore_cwise, notrans, rowequ : bool = ---;
    norm_char : u8 = ---;

    info.* = 0;
    trans_type = ilatrans(trans);
    ref_type = 1;
    if nparams >= la_linrx_itref_i {
        if fget(params.*, la_linrx_itref_i) < zero_r {
            fset(params, la_linrx_itref_i, one_r);
        } else {
            ref_type = cast(int) fget(params.*, la_linrx_itref_i);
        }
    }

    illrcond_thresh = ncast(RT, n) * _gbrfsx_lamch_epsilon(RT);
    ithresh = 10;
    rthresh = half_r;
    unstable_thresh = quarter_r;
    ignore_cwise = false;

    if nparams >= la_linrx_ithresh_i {
        if fget(params.*, la_linrx_ithresh_i) < zero_r {
            fset(params, la_linrx_ithresh_i, ncast(RT, ithresh));
        } else {
            ithresh = cast(int) fget(params.*, la_linrx_ithresh_i);
        }
    }
    if nparams >= la_linrx_cwise_i {
        if fget(params.*, la_linrx_cwise_i) < zero_r {
            fset(params, la_linrx_cwise_i, one_r);
        } else {
            ignore_cwise = fget(params.*, la_linrx_cwise_i) == zero_r;
        }
    }

    if nparams >= la_linrx_cwise_i {
        cwise_param = fget(params.*, la_linrx_cwise_i);
    } else {
        cwise_param = ifx ignore_cwise then zero_r else one_r;
    }

    if ref_type == 0 || n_err_bnds == 0 {
        n_norms = 0;
    } else if ignore_cwise {
        n_norms = 1;
    } else {
        n_norms = 2;
    }

    notrans = char_eq(trans, #char "n");
    rowequ = char_eq(equed, #char "r") || char_eq(equed, #char "b");
    colequ = char_eq(equed, #char "c") || char_eq(equed, #char "b");

    if trans_type == -1 {
        info.* = -1;
    } else if !rowequ && !colequ && !char_eq(equed, #char "n") {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if kl < 0 {
        info.* = -4;
    } else if ku < 0 {
        info.* = -5;
    } else if nrhs < 0 {
        info.* = -6;
    } else if ldab < kl + ku + 1 {
        info.* = -8;
    } else if ldafb < 2*kl + ku + 1 {
        info.* = -10;
    } else if ldb < max(1, n) {
        info.* = -13;
    } else if ldx < max(1, n) {
        info.* = -15;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 || nrhs == 0 {
        rcond.* = one_r;
        for j:1..nrhs {
            fset(berr, j, zero_r);
            if n_err_bnds >= 1 {
                fset(err_bnds_norm, j, la_linrx_trust_i, one_r);
                fset(err_bnds_comp, j, la_linrx_trust_i, one_r);
            }
            if n_err_bnds >= 2 {
                fset(err_bnds_norm, j, la_linrx_err_i, zero_r);
                fset(err_bnds_comp, j, la_linrx_err_i, zero_r);
            }
            if n_err_bnds >= 3 {
                fset(err_bnds_norm, j, la_linrx_rcond_i, one_r);
                fset(err_bnds_comp, j, la_linrx_rcond_i, one_r);
            }
        }
        return;
    }

    rcond.* = zero_r;
    for j:1..nrhs {
        fset(berr, j, one_r);
        if n_err_bnds >= 1 {
            fset(err_bnds_norm, j, la_linrx_trust_i, one_r);
            fset(err_bnds_comp, j, la_linrx_trust_i, one_r);
        }
        if n_err_bnds >= 2 {
            fset(err_bnds_norm, j, la_linrx_err_i, one_r);
            fset(err_bnds_comp, j, la_linrx_err_i, one_r);
        }
        if n_err_bnds >= 3 {
            fset(err_bnds_norm, j, la_linrx_rcond_i, zero_r);
            fset(err_bnds_comp, j, la_linrx_rcond_i, zero_r);
        }
    }

    norm_char = ifx notrans then #char "i" else #char "1";
    #if MT == RT {
        anorm = _gbrfsx_langb(norm_char, n, kl, ku, ab, ldab, work, RT, MT);
    } else {
        anorm = _gbrfsx_langb(norm_char, n, kl, ku, ab, ldab, aux, RT, MT);
    }
    _gbrfsx_gbcon(norm_char, n, kl, ku, afb, ldafb, ipiv, anorm, rcond, work, aux, info, RT, MT);

    if ref_type != 0 && info.* == 0 {
        prec_type = _gbrfsx_ilaprec(RT);
        if notrans {
            _gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, x, ldx, berr, n_norms, err_bnds_norm, err_bnds_comp, work, aux, rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info, RT, MT);
        } else {
            _gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, rowequ, r, b, ldb, x, ldx, berr, n_norms, err_bnds_norm, err_bnds_comp, work, aux, rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info, RT, MT);
        }
    }

    err_lbnd = max(ten_r, sqrt(ncast(RT, n))) * _gbrfsx_lamch_epsilon(RT);

    if n_err_bnds >= 1 && n_norms >= 1 {
        if colequ && notrans {
            rcond_tmp = _gbrfsx_gbrcond_norm(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, -1, c, info, work, aux, RT, MT);
        } else if rowequ && !notrans {
            rcond_tmp = _gbrfsx_gbrcond_norm(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, -1, r, info, work, aux, RT, MT);
        } else {
            rcond_tmp = _gbrfsx_gbrcond_norm(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, 0, r, info, work, aux, RT, MT);
        }

        for j:1..nrhs {
            if n_err_bnds >= la_linrx_err_i && fget(err_bnds_norm.*, j, la_linrx_err_i) > one_r {
                fset(err_bnds_norm, j, la_linrx_err_i, one_r);
            }

            if rcond_tmp < illrcond_thresh {
                fset(err_bnds_norm, j, la_linrx_err_i, one_r);
                fset(err_bnds_norm, j, la_linrx_trust_i, zero_r);
                if info.* <= n then info.* = n + j;
            } else if fget(err_bnds_norm.*, j, la_linrx_err_i) < err_lbnd {
                fset(err_bnds_norm, j, la_linrx_err_i, err_lbnd);
                fset(err_bnds_norm, j, la_linrx_trust_i, one_r);
            }

            if n_err_bnds >= la_linrx_rcond_i {
                fset(err_bnds_norm, j, la_linrx_rcond_i, rcond_tmp);
            }
        }
    }

    if n_err_bnds >= 1 && n_norms >= 2 {
        cwise_wrong = sqrt(_gbrfsx_lamch_epsilon(RT));
        for j:1..nrhs {
            if fget(err_bnds_comp.*, j, la_linrx_err_i) < cwise_wrong {
                x_view := _gbrfsx_matrix_column(x.*, ldx, j, n);
                rcond_tmp = _gbrfsx_gbrcond_x(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, x_view, info, work, aux, RT, MT);
            } else {
                rcond_tmp = zero_r;
            }

            if n_err_bnds >= la_linrx_err_i && fget(err_bnds_comp.*, j, la_linrx_err_i) > one_r {
                fset(err_bnds_comp, j, la_linrx_err_i, one_r);
            }

            if rcond_tmp < illrcond_thresh {
                fset(err_bnds_comp, j, la_linrx_err_i, one_r);
                fset(err_bnds_comp, j, la_linrx_trust_i, zero_r);
                if cwise_param == one_r && info.* < n + j then info.* = n + j;
            } else if fget(err_bnds_comp.*, j, la_linrx_err_i) < err_lbnd {
                fset(err_bnds_comp, j, la_linrx_err_i, err_lbnd);
                fset(err_bnds_comp, j, la_linrx_trust_i, one_r);
            }

            if n_err_bnds >= la_linrx_rcond_i {
                fset(err_bnds_comp, j, la_linrx_rcond_i, rcond_tmp);
            }
        }
    }
}
