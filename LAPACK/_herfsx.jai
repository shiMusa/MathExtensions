// ? LAPACK *herfsx()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cherfsx/zherfsx
// *HERFSX improves the computed solution to a system of linear
// equations when the coefficient matrix is Hermitian indefinite.

cherfsx :: (
    $$uplo: u8,
    $$equed: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    s: $S/ArrayType1D,
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
    _herfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float32, MT=C32);
}

zherfsx :: (
    $$uplo: u8,
    $$equed: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    s: $S/ArrayType1D,
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
    _herfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_herfsx_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_herfsx_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_herfsx_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
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

_herfsx_complex_scratch :: inline (x: $X/ArrayType1D, $$count: int, $RT: Type, $MT: Type) -> [] MT #modify {
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

_herfsx_ilaprec :: ($RT: Type) -> int #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return ilaprec(#char "d");
    } else {
        return ilaprec(#char "e");
    }
}

_herfsx_lanhe :: (
    $$norm: u8,
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    aux: *$AUX/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           AUX.data_type == RT;
} {
    #if MT == C32 {
        return clanhe(norm, uplo, n, a, lda, aux);
    } else {
        return zlanhe(norm, uplo, n, a, lda, aux);
    }
}

_herfsx_hecon :: (
    $$uplo: u8,
    $$n: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    anorm: $RT,
    rcond: *RT,
    work: *$WORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           WORK.data_type == MT;
} {
    #if MT == C32 {
        checon(uplo, n, af, ldaf, ipiv, anorm, rcond, work, info);
    } else {
        zhecon(uplo, n, af, ldaf, ipiv, anorm, rcond, work, info);
    }
}

_herfsx_extended :: (
    $$prec_type: int,
    $$uplo: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    $$rcequ: bool,
    s: $S/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    berr: *$BERR/ArrayType1D,
    $$n_norms: int,
    err_bnds_norm: *$ERR_BNDS_NORM/ArrayType2D,
    err_bnds_comp: *$ERR_BNDS_COMP/ArrayType2D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
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
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           S.data_type == RT &&
           B.data_type == MT &&
           X.data_type == MT &&
           BERR.data_type == RT &&
           ERR_BNDS_NORM.data_type == RT &&
           ERR_BNDS_COMP.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    dy_view := _herfsx_subvector(work.*, n + 1, n);
    y_tail_view := _herfsx_complex_scratch(rwork.*, n, RT, MT);
    #if MT == C32 {
        cla_herfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, rcequ, s, b, ldb, x, ldx, berr, n_norms, err_bnds_norm, err_bnds_comp, work, rwork, *dy_view, *y_tail_view, rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info);
    } else {
        zla_herfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, rcequ, s, b, ldb, x, ldx, berr, n_norms, err_bnds_norm, err_bnds_comp, work, rwork, *dy_view, *y_tail_view, rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info);
    }
}

_herfsx_hercond_c :: (
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    s: $S/ArrayType1D,
    $$apply_scale: bool,
    info: *int,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           S.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    #if MT == C32 {
        return cla_hercond_c(uplo, n, a, lda, af, ldaf, ipiv, s, apply_scale, info, work, rwork);
    } else {
        return zla_hercond_c(uplo, n, a, lda, af, ldaf, ipiv, s, apply_scale, info, work, rwork);
    }
}

_herfsx_hercond_x :: (
    $$uplo: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    x: $X/ArrayType1D,
    info: *int,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    $RT: Type,
    $MT: Type
) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           X.data_type == MT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    #if MT == C32 {
        return cla_hercond_x(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork);
    } else {
        return zla_hercond_x(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork);
    }
}

_herfsx :: (
    $$uplo: u8,
    $$equed: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    s: $S/ArrayType1D,
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
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           S.data_type == RT &&
           B.data_type == MT &&
           X.data_type == MT &&
           BERR.data_type == RT &&
           ERR_BNDS_NORM.data_type == RT &&
           ERR_BNDS_COMP.data_type == RT &&
           PARAMS.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    half_r :: #run ncast(RT, 0.5);
    quarter_r :: #run ncast(RT, 0.25);
    ten_r :: #run ncast(RT, 10);
    routine_name :: #run if MT == C32 { "CHERFSX" } else { "ZHERFSX" };

    la_linrx_itref_i :: 1;
    la_linrx_ithresh_i :: 2;
    la_linrx_cwise_i :: 3;
    la_linrx_trust_i :: 1;
    la_linrx_err_i :: 2;
    la_linrx_rcond_i :: 3;

    anorm, cwise_wrong, err_lbnd, illrcond_thresh, rcond_tmp, rthresh, unstable_thresh : RT = ---;
    ithresh, j, n_norms, prec_type, ref_type : int = ---;
    ignore_cwise, rcequ, upper : bool = ---;
    norm_char : u8 = ---;

    info.* = 0;
    ref_type = 1;
    if nparams >= la_linrx_itref_i {
        if fget(params.*, la_linrx_itref_i) < zero_r {
            fset(params, la_linrx_itref_i, one_r);
        } else {
            ref_type = cast(int) fget(params.*, la_linrx_itref_i);
        }
    }

    illrcond_thresh = ncast(RT, n) * _herfsx_lamch_epsilon(RT);
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

    if ref_type == 0 || n_err_bnds == 0 {
        n_norms = 0;
    } else if ignore_cwise {
        n_norms = 1;
    } else {
        n_norms = 2;
    }

    rcequ = lsame(equed, #char "Y");
    upper = lsame(uplo, #char "U");

    if !upper && !lsame(uplo, #char "L") {
        info.* = -1;
    } else if !rcequ && !lsame(equed, #char "N") {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if nrhs < 0 {
        info.* = -4;
    } else if lda < max(1, n) {
        info.* = -6;
    } else if ldaf < max(1, n) {
        info.* = -8;
    } else if ldb < max(1, n) {
        info.* = -12;
    } else if ldx < max(1, n) {
        info.* = -14;
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

    norm_char = #char "I";
    anorm = _herfsx_lanhe(norm_char, uplo, n, a, lda, rwork, RT, MT);
    _herfsx_hecon(uplo, n, af, ldaf, ipiv, anorm, rcond, work, info, RT, MT);

    if ref_type != 0 {
        prec_type = _herfsx_ilaprec(RT);
        _herfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, rcequ, s, b, ldb, x, ldx, berr, n_norms, err_bnds_norm, err_bnds_comp, work, rwork, rcond, ithresh, rthresh, unstable_thresh, ignore_cwise, info, RT, MT);
    }

    err_lbnd = max(ten_r, sqrt(ncast(RT, n))) * _herfsx_lamch_epsilon(RT);

    if n_err_bnds >= 1 && n_norms >= 1 {
        rcond_tmp = _herfsx_hercond_c(uplo, n, a, lda, af, ldaf, ipiv, s, rcequ, info, work, rwork, RT, MT);
        for j:1..nrhs {
            if n_err_bnds >= la_linrx_err_i && fget(err_bnds_norm.*, j, la_linrx_err_i) > one_r {
                fset(err_bnds_norm, j, la_linrx_err_i, one_r);
            }

            if rcond_tmp < illrcond_thresh {
                fset(err_bnds_norm, j, la_linrx_err_i, one_r);
                fset(err_bnds_norm, j, la_linrx_trust_i, zero_r);
                if info.* <= n {
                    info.* = n + j;
                }
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
        cwise_wrong = sqrt(_herfsx_lamch_epsilon(RT));
        for j:1..nrhs {
            if fget(err_bnds_comp.*, j, la_linrx_err_i) < cwise_wrong {
                x_view := _herfsx_matrix_column(x.*, ldx, j, n);
                rcond_tmp = _herfsx_hercond_x(uplo, n, a, lda, af, ldaf, ipiv, x_view, info, work, rwork, RT, MT);
            } else {
                rcond_tmp = zero_r;
            }

            if n_err_bnds >= la_linrx_err_i && fget(err_bnds_comp.*, j, la_linrx_err_i) > one_r {
                fset(err_bnds_comp, j, la_linrx_err_i, one_r);
            }

            if rcond_tmp < illrcond_thresh {
                fset(err_bnds_comp, j, la_linrx_err_i, one_r);
                fset(err_bnds_comp, j, la_linrx_trust_i, zero_r);
                if !ignore_cwise && info.* < n + j {
                    info.* = n + j;
                }
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
