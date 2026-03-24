// ? LAPACK *la_gbrfsx_extended()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sla_gbrfsx_extended/cla_gbrfsx_extended/dla_gbrfsx_extended/zla_gbrfsx_extended
// *LA_GBRFSX_EXTENDED improves the computed solution to a system of
// linear equations for general banded matrices by performing extra-precise
// iterative refinement and provides error bounds and backward error
// estimates for the solution.

sla_gbrfsx_extended :: #bake_arguments _la_gbrfsx_extended(RT=float32, MT=float32);
dla_gbrfsx_extended :: #bake_arguments _la_gbrfsx_extended(RT=float64, MT=float64);
cla_gbrfsx_extended :: #bake_arguments _la_gbrfsx_extended(RT=float32, MT=C32);
zla_gbrfsx_extended :: #bake_arguments _la_gbrfsx_extended(RT=float64, MT=C64);

#scope_file
_la_gbrfsx_extended_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
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

_la_gbrfsx_extended_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_la_gbrfsx_extended_lamch_overflow :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "o");
    } else {
        return dlamch(#char "o");
    }
}

_la_gbrfsx_extended_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_la_gbrfsx_extended_copy :: (
    $$n: int,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        scopy(n, x, 1, y, 1);
    } else if MT == float64 {
        dcopy(n, x, 1, y, 1);
    } else if MT == C32 {
        ccopy(n, x, 1, y, 1);
    } else {
        zcopy(n, x, 1, y, 1);
    }
}

_la_gbrfsx_extended_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        saxpy(n, alpha, x, 1, y, 1);
    } else if MT == float64 {
        daxpy(n, alpha, x, 1, y, 1);
    } else if MT == C32 {
        caxpy(n, alpha, x, 1, y, 1);
    } else {
        zaxpy(n, alpha, x, 1, y, 1);
    }
}

_la_gbrfsx_extended_gbmv :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    alpha: $MT,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    x: $X/ArrayType1D,
    beta: MT,
    y: *$Y/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AB.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        sgbmv(trans, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1);
    } else if MT == float64 {
        dgbmv(trans, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1);
    } else if MT == C32 {
        cgbmv(trans, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1);
    } else {
        zgbmv(trans, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1);
    }
}

_la_gbrfsx_extended_blas_gbmv_x :: (
    $$trans_type: int,
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    alpha: $MT,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    x: $X/ArrayType1D,
    beta: MT,
    y: *$Y/ArrayType1D,
    $$prec_type: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AB.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        BLAS_SGBMV_X(trans_type, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1, prec_type);
    } else if MT == float64 {
        BLAS_DGBMV_X(trans_type, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1, prec_type);
    } else if MT == C32 {
        BLAS_CGBMV_X(trans_type, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1, prec_type);
    } else {
        BLAS_ZGBMV_X(trans_type, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1, prec_type);
    }
}

_la_gbrfsx_extended_blas_gbmv2_x :: (
    $$trans_type: int,
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    alpha: $MT,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    x_head: $XH/ArrayType1D,
    x_tail: $XT/ArrayType1D,
    beta: MT,
    y: *$Y/ArrayType1D,
    $$prec_type: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AB.data_type == MT && XH.data_type == MT && XT.data_type == MT && Y.data_type == MT;
} {
    #if MT == float32 {
        BLAS_SGBMV2_X(trans_type, m, n, kl, ku, alpha, ab, ldab, x_head, x_tail, 1, beta, y, 1, prec_type);
    } else if MT == float64 {
        BLAS_DGBMV2_X(trans_type, m, n, kl, ku, alpha, ab, ldab, x_head, x_tail, 1, beta, y, 1, prec_type);
    } else if MT == C32 {
        BLAS_CGBMV2_X(trans_type, m, n, kl, ku, alpha, ab, ldab, x_head, x_tail, 1, beta, y, 1, prec_type);
    } else {
        BLAS_ZGBMV2_X(trans_type, m, n, kl, ku, alpha, ab, ldab, x_head, x_tail, 1, beta, y, 1, prec_type);
    }
}

_la_gbrfsx_extended_gbtrs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType1D,
    $$ldb: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           AFB.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    #if MT == float32 {
        sgbtrs(trans, n, kl, ku, 1, afb, ldafb, ipiv, b, ldb, info);
    } else if MT == float64 {
        dgbtrs(trans, n, kl, ku, 1, afb, ldafb, ipiv, b, ldb, info);
    } else if MT == C32 {
        cgbtrs(trans, n, kl, ku, 1, afb, ldafb, ipiv, b, ldb, info);
    } else {
        zgbtrs(trans, n, kl, ku, 1, afb, ldafb, ipiv, b, ldb, info);
    }
}

_la_gbrfsx_extended_wwaddw :: (
    $$n: int,
    x: *$X/ArrayType1D,
    y: *$Y/ArrayType1D,
    w: $W/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT && Y.data_type == MT && W.data_type == MT;
} {
    #if MT == float32 {
        sla_wwaddw(n, x, y, w);
    } else if MT == float64 {
        dla_wwaddw(n, x, y, w);
    } else if MT == C32 {
        cla_wwaddw(n, x, y, w);
    } else {
        zla_wwaddw(n, x, y, w);
    }
}

_la_gbrfsx_extended_gbamv :: (
    $$trans_type: int,
    $$m: int,
    $$n: int,
    $$kl: int,
    $$ku: int,
    alpha: $RT,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    x: $X/ArrayType1D,
    beta: RT,
    y: *$Y/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AB.data_type == MT && X.data_type == MT && Y.data_type == RT;
} {
    #if MT == float32 {
        sla_gbamv(trans_type, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1);
    } else if MT == float64 {
        dla_gbamv(trans_type, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1);
    } else if MT == C32 {
        cla_gbamv(trans_type, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1);
    } else {
        zla_gbamv(trans_type, m, n, kl, ku, alpha, ab, ldab, x, 1, beta, y, 1);
    }
}

_la_gbrfsx_extended_lin_berr :: (
    $$n: int,
    $$nz: int,
    $$nrhs: int,
    res: $RES/ArrayType1D,
    ayb: $AYB/ArrayType1D,
    berr: *$RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           RES.data_type == MT && AYB.data_type == RT;
} {
    berr_buf: [1] RT;
    berr_buf[0] = berr.*;

    #if MT == float32 {
        sla_lin_berr(n, nz, nrhs, res, ayb, berr_buf);
    } else if MT == float64 {
        dla_lin_berr(n, nz, nrhs, res, ayb, berr_buf);
    } else if MT == C32 {
        cla_lin_berr(n, nz, nrhs, res, ayb, berr_buf);
    } else {
        zla_lin_berr(n, nz, nrhs, res, ayb, berr_buf);
    }

    berr.* = berr_buf[0];
}

_la_gbrfsx_extended :: (
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
    res: *$RES/ArrayType1D,
    ayb: *$AYB/ArrayType1D,
    dy: *$DY/ArrayType1D,
    y_tail: *$Y_TAIL/ArrayType1D,
    rcond: RT,
    $$ithresh: int,
    rthresh: RT,
    dz_ub: RT,
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
           RES.data_type == MT &&
           AYB.data_type == RT &&
           DY.data_type == MT &&
           Y_TAIL.data_type == MT;
} {
    unstable_state :: 0;
    working_state  :: 1;
    conv_state     :: 2;
    noprog_state   :: 3;
    base_residual  :: 0;
    extra_residual :: 1;
    extra_y        :: 2;
    la_linrx_err_i :: 2;

    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    zero_t :: #run zero(MT);
    one_t  :: #run one(MT);

    cnt, i, j, m, x_state, z_state, y_prec_state : int = ---;
    dx_x, dxrat, dxratmax, dyk, dz_z, dzrat, dzratmax, eps, final_dx_x, final_dz_z,
    hugeval, incr_thresh, normdx, normx, normy, prev_dz_z, prevnormdx, yk, ymin, berr_j : RT = ---;
    incr_prec : bool = ---;
    trans : u8 = ---;

    if info.* != 0 then return;

    trans = chla_transtype(trans_type);
    eps = _la_gbrfsx_extended_lamch_epsilon(RT);
    hugeval = _la_gbrfsx_extended_lamch_overflow(RT);
    hugeval *= hugeval;
    incr_thresh = ncast(RT, n) * eps;
    m = kl + ku + 1;

    minus_one_t := -one_t;

    for j:1..nrhs {
        y_prec_state = extra_residual;
        if y_prec_state == extra_y {
            for i:1..n {
                fset(y_tail, i, zero_t);
            }
        }

        dxrat = zero_r;
        dxratmax = zero_r;
        dzrat = zero_r;
        dzratmax = zero_r;
        final_dx_x = hugeval;
        final_dz_z = hugeval;
        prevnormdx = hugeval;
        prev_dz_z = hugeval;
        dz_z = hugeval;
        dx_x = hugeval;

        x_state = working_state;
        z_state = unstable_state;
        incr_prec = false;

        b_col := _la_gbrfsx_extended_matrix_column(b, ldb, j, n);
        y_col := _la_gbrfsx_extended_matrix_column(y.*, ldy, j, n);

        for cnt:1..ithresh {
            _la_gbrfsx_extended_copy(n, b_col, res, MT);

            if y_prec_state == base_residual {
                _la_gbrfsx_extended_gbmv(trans, m, n, kl, ku, minus_one_t, ab, ldab, y_col, one_t, res, MT);
            } else if y_prec_state == extra_residual {
                _la_gbrfsx_extended_blas_gbmv_x(trans_type, n, n, kl, ku, minus_one_t, ab, ldab, y_col, one_t, res, prec_type, MT);
            } else {
                _la_gbrfsx_extended_blas_gbmv2_x(trans_type, n, n, kl, ku, minus_one_t, ab, ldab, y_col, y_tail.*, one_t, res, prec_type, MT);
            }

            _la_gbrfsx_extended_copy(n, res.*, dy, MT);
            _la_gbrfsx_extended_gbtrs(trans, n, kl, ku, afb, ldafb, ipiv, dy, n, info, MT);

            normx = zero_r;
            normy = zero_r;
            normdx = zero_r;
            dz_z = zero_r;
            ymin = hugeval;

            for i:1..n {
                yk = _la_gbrfsx_extended_abs1(fget(y_col, i), RT, MT);
                dyk = _la_gbrfsx_extended_abs1(fget(dy.*, i), RT, MT);

                if yk != zero_r {
                    dz_z = max(dz_z, dyk / yk);
                } else if dyk != zero_r {
                    dz_z = hugeval;
                }

                ymin = min(ymin, yk);
                normy = max(normy, yk);

                if colequ {
                    normx = max(normx, yk * fget(c, i));
                    normdx = max(normdx, dyk * fget(c, i));
                } else {
                    normx = normy;
                    normdx = max(normdx, dyk);
                }
            }

            if normx != zero_r {
                dx_x = normdx / normx;
            } else if normdx == zero_r {
                dx_x = zero_r;
            } else {
                dx_x = hugeval;
            }

            dxrat = normdx / prevnormdx;
            dzrat = dz_z / prev_dz_z;

            if !ignore_cwise && ymin*rcond < incr_thresh*normy && y_prec_state < extra_y {
                incr_prec = true;
            }

            if x_state == noprog_state && dxrat <= rthresh {
                x_state = working_state;
            }
            if x_state == working_state {
                if dx_x <= eps {
                    x_state = conv_state;
                } else if dxrat > rthresh {
                    if y_prec_state != extra_y {
                        incr_prec = true;
                    } else {
                        x_state = noprog_state;
                    }
                } else if dxrat > dxratmax {
                    dxratmax = dxrat;
                }
                if x_state > working_state {
                    final_dx_x = dx_x;
                }
            }

            if z_state == unstable_state && dz_z <= dz_ub {
                z_state = working_state;
            }
            if z_state == noprog_state && dzrat <= rthresh {
                z_state = working_state;
            }
            if z_state == working_state {
                if dz_z <= eps {
                    z_state = conv_state;
                } else if dz_z > dz_ub {
                    z_state = unstable_state;
                    dzratmax = zero_r;
                    final_dz_z = hugeval;
                } else if dzrat > rthresh {
                    if y_prec_state != extra_y {
                        incr_prec = true;
                    } else {
                        z_state = noprog_state;
                    }
                } else if dzrat > dzratmax {
                    dzratmax = dzrat;
                }
                if z_state > working_state {
                    final_dz_z = dz_z;
                }
            }

            if x_state != working_state {
                if ignore_cwise then break;
                if z_state == noprog_state || z_state == conv_state then break;
                if z_state == unstable_state && cnt > 1 then break;
            }

            if incr_prec {
                incr_prec = false;
                y_prec_state += 1;
                for i:1..n {
                    fset(y_tail, i, zero_t);
                }
            }

            prevnormdx = normdx;
            prev_dz_z = dz_z;

            if y_prec_state < extra_y {
                _la_gbrfsx_extended_axpy(n, one_t, dy.*, *y_col, MT);
            } else {
                _la_gbrfsx_extended_wwaddw(n, *y_col, y_tail, dy.*, MT);
            }
        }

        if x_state == working_state then final_dx_x = dx_x;
        if z_state == working_state then final_dz_z = dz_z;

        if n_norms >= 1 {
            fset(err_bnds_norm, j, la_linrx_err_i, final_dx_x / (one_r - dxratmax));
        }
        if n_norms >= 2 {
            fset(err_bnds_comp, j, la_linrx_err_i, final_dz_z / (one_r - dzratmax));
        }

        _la_gbrfsx_extended_copy(n, b_col, res, MT);
        _la_gbrfsx_extended_gbmv(trans, n, n, kl, ku, minus_one_t, ab, ldab, y_col, one_t, res, MT);

        for i:1..n {
            fset(ayb, i, _la_gbrfsx_extended_abs1(fget(b_col, i), RT, MT));
        }

        _la_gbrfsx_extended_gbamv(trans_type, n, n, kl, ku, one_r, ab, ldab, y_col, one_r, ayb, RT, MT);
        berr_j = ---;
        _la_gbrfsx_extended_lin_berr(n, n, 1, res.*, ayb.*, *berr_j, RT, MT);
        fset(berr_out, j, berr_j);
    }
}
