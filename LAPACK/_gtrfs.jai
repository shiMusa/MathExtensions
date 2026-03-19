// ? LAPACK *gtrfs()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgtrfs/cgtrfs/dgtrfs/zgtrfs
// *GTRFS improves the computed solution to a tridiagonal system and
// provides error bounds and backward error estimates.

sgtrfs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: $DLF/ArrayType1D,
    df: $DF/ArrayType1D,
    duf: $DUF/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gtrfs_real(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, RT=float32);
}

dgtrfs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: $DLF/ArrayType1D,
    df: $DF/ArrayType1D,
    duf: $DUF/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gtrfs_real(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, RT=float64);
}

cgtrfs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: $DLF/ArrayType1D,
    df: $DF/ArrayType1D,
    duf: $DUF/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gtrfs_complex(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float32, MT=C32);
}

zgtrfs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: $DLF/ArrayType1D,
    df: $DF/ArrayType1D,
    duf: $DUF/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int
) {
    _gtrfs_complex(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gtrfs_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gtrfs_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
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

_gtrfs_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_gtrfs_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gtrfs_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_gtrfs_lacn2 :: (
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

_gtrfs_gttrs :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dlf: $DLF/ArrayType1D,
    df: $DF/ArrayType1D,
    duf: $DUF/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DLF.data_type == MT &&
           DF.data_type == MT &&
           DUF.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    #if MT == float32 {
        sgttrs(trans, n, nrhs, dlf, df, duf, du2, ipiv, b, ldb, info);
    } else if MT == float64 {
        dgttrs(trans, n, nrhs, dlf, df, duf, du2, ipiv, b, ldb, info);
    } else if MT == C32 {
        cgttrs(trans, n, nrhs, dlf, df, duf, du2, ipiv, b, ldb, info);
    } else {
        zgttrs(trans, n, nrhs, dlf, df, duf, du2, ipiv, b, ldb, info);
    }
}

_gtrfs_copy :: (
    $$n: int,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT;
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

_gtrfs_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           Y.data_type == MT;
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

_gtrfs_lagtm :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    alpha: $RT,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    x: $X/ArrayType2D,
    $$ldx: int,
    beta: $RT,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           X.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == float32 {
        slagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
    } else if MT == float64 {
        dlagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
    } else if MT == C32 {
        clagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
    } else {
        zlagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
    }
}

_gtrfs_real :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: $DLF/ArrayType1D,
    df: $DF/ArrayType1D,
    duf: $DUF/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    iwork: *$IWORK/ArrayType1D,
    info: *int,
    $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           DL.data_type == RT &&
           D.data_type == RT &&
           DU.data_type == RT &&
           DLF.data_type == RT &&
           DF.data_type == RT &&
           DUF.data_type == RT &&
           DU2.data_type == RT &&
           IPIV.data_type == int &&
           B.data_type == RT &&
           X.data_type == RT &&
           FERR.data_type == RT &&
           BERR.data_type == RT &&
           WORK.data_type == RT &&
           IWORK.data_type == int;
} {
    itmax :: 5;
    zero_r  :: #run zero(RT);
    one_r   :: #run one(RT);
    two_r   :: #run ncast(RT, 2);
    three_r :: #run ncast(RT, 3);
    neg_one_r :: #run -one(RT);
    routine_name :: #run if RT == float32 { "SGTRFS" } else { "DGTRFS" };

    notran : bool = ---;
    count, i, j, kase, nz : int = ---;
    eps, ferr_j, lstres, s, safe1, safe2, safmin, xnorm : RT = ---;
    transn, transt : u8 = ---;

    info.* = 0;
    notran = lsame(trans, #char "N");
    if !notran && !lsame(trans, #char "T") && !lsame(trans, #char "C") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
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
        for j:1..nrhs {
            fset(ferr, j, zero_r);
            fset(berr, j, zero_r);
        }
        return;
    }

    if notran {
        transn = #char "N";
        transt = #char "T";
    } else {
        transn = #char "T";
        transt = #char "N";
    }

    nz = 4;
    eps = _gtrfs_lamch_epsilon(RT);
    safmin = _gtrfs_lamch_safe_minimum(RT);
    safe1 = ncast(RT, nz)*safmin;
    safe2 = safe1 / eps;

    weights_view := _gtrfs_subvector(work.*, 1, n);
    active_view  := _gtrfs_subvector(work.*, n + 1, n);
    v_view       := _gtrfs_subvector(work.*, 2*n + 1, n);

    active_mat := array2d_view(active_view, max(1, n), row_major=false);

    for j:1..nrhs {
        b_view := _gtrfs_matrix_column(b, ldb, j, n);
        x_view := _gtrfs_matrix_column(x.*, ldx, j, n);
        x_mat := array2d_view(x_view, max(1, n), row_major=false);

        count = 1;
        lstres = three_r;

        while true {
            _gtrfs_copy(n, b_view, *active_view, RT, RT);
            _gtrfs_lagtm(trans, n, 1, neg_one_r, dl, d, du, *x_mat, n, one_r, *active_mat, n, RT, RT);

            if notran {
                if n == 1 {
                    fset(weights_view, 1, abs(fget(b_view, 1)) + abs(fget(d, 1)*fget(x_view, 1)));
                } else {
                    fset(weights_view, 1, abs(fget(b_view, 1)) + abs(fget(d, 1)*fget(x_view, 1)) + abs(fget(du, 1)*fget(x_view, 2)));
                    for i:2..n - 1 {
                        fset(weights_view, i, abs(fget(b_view, i)) +
                            abs(fget(dl, i - 1)*fget(x_view, i - 1)) +
                            abs(fget(d, i)*fget(x_view, i)) +
                            abs(fget(du, i)*fget(x_view, i + 1)));
                    }
                    fset(weights_view, n, abs(fget(b_view, n)) +
                        abs(fget(dl, n - 1)*fget(x_view, n - 1)) +
                        abs(fget(d, n)*fget(x_view, n)));
                }
            } else {
                if n == 1 {
                    fset(weights_view, 1, abs(fget(b_view, 1)) + abs(fget(d, 1)*fget(x_view, 1)));
                } else {
                    fset(weights_view, 1, abs(fget(b_view, 1)) + abs(fget(d, 1)*fget(x_view, 1)) + abs(fget(dl, 1)*fget(x_view, 2)));
                    for i:2..n - 1 {
                        fset(weights_view, i, abs(fget(b_view, i)) +
                            abs(fget(du, i - 1)*fget(x_view, i - 1)) +
                            abs(fget(d, i)*fget(x_view, i)) +
                            abs(fget(dl, i)*fget(x_view, i + 1)));
                    }
                    fset(weights_view, n, abs(fget(b_view, n)) +
                        abs(fget(du, n - 1)*fget(x_view, n - 1)) +
                        abs(fget(d, n)*fget(x_view, n)));
                }
            }

            s = zero_r;
            for i:1..n {
                if fget(weights_view, i) > safe2 {
                    s = max(s, abs(fget(active_view, i)) / fget(weights_view, i));
                } else {
                    s = max(s, (abs(fget(active_view, i)) + safe1) / (fget(weights_view, i) + safe1));
                }
            }
            fset(berr, j, s);

            if fget(berr.*, j) > eps && two_r*fget(berr.*, j) <= lstres && count <= itmax {
                _gtrfs_gttrs(trans, n, 1, dlf, df, duf, du2, ipiv, *active_mat, n, info, RT, RT);
                _gtrfs_axpy(n, one_r, active_view, *x_view, RT, RT);
                lstres = fget(berr.*, j);
                count += 1;
                continue;
            }

            for i:1..n {
                if fget(weights_view, i) > safe2 {
                    fset(weights_view, i, abs(fget(active_view, i)) + ncast(RT, nz)*eps*fget(weights_view, i));
                } else {
                    fset(weights_view, i, abs(fget(active_view, i)) + ncast(RT, nz)*eps*fget(weights_view, i) + safe1);
                }
            }

            ferr_j = zero_r;
            kase = 0;
            isave := array1d(3, 0, 0, 0);

            while true {
                _gtrfs_lacn2(n, *v_view, *active_view, iwork, *ferr_j, *kase, *isave, RT, RT);
                if kase == 0 break;

                if kase == 1 {
                    _gtrfs_gttrs(transt, n, 1, dlf, df, duf, du2, ipiv, *active_mat, n, info, RT, RT);
                    for i:1..n {
                        fset(active_view, i, fget(weights_view, i)*fget(active_view, i));
                    }
                } else {
                    for i:1..n {
                        fset(active_view, i, fget(weights_view, i)*fget(active_view, i));
                    }
                    _gtrfs_gttrs(transn, n, 1, dlf, df, duf, du2, ipiv, *active_mat, n, info, RT, RT);
                }
            }

            xnorm = zero_r;
            for i:1..n {
                xnorm = max(xnorm, abs(fget(x_view, i)));
            }
            if xnorm != zero_r {
                ferr_j /= xnorm;
            }
            fset(ferr, j, ferr_j);
            break;
        }
    }
}

_gtrfs_complex :: (
    $$trans: u8,
    $$n: int,
    $$nrhs: int,
    dl: $DL/ArrayType1D,
    d: $D/ArrayType1D,
    du: $DU/ArrayType1D,
    dlf: $DLF/ArrayType1D,
    df: $DF/ArrayType1D,
    duf: $DUF/ArrayType1D,
    du2: $DU2/ArrayType1D,
    ipiv: $IPIV/ArrayType1D,
    b: $B/ArrayType2D,
    $$ldb: int,
    x: *$X/ArrayType2D,
    $$ldx: int,
    ferr: *$FERR/ArrayType1D,
    berr: *$BERR/ArrayType1D,
    work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           DL.data_type == MT &&
           D.data_type == MT &&
           DU.data_type == MT &&
           DLF.data_type == MT &&
           DF.data_type == MT &&
           DUF.data_type == MT &&
           DU2.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT &&
           X.data_type == MT &&
           FERR.data_type == RT &&
           BERR.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    itmax :: 5;
    zero_r  :: #run zero(RT);
    one_r   :: #run one(RT);
    one_t   :: #run cpx(one(RT), zero(RT));
    neg_one_r :: #run -one(RT);
    two_r   :: #run ncast(RT, 2);
    three_r :: #run ncast(RT, 3);
    routine_name :: #run if MT == C32 { "CGTRFS" } else { "ZGTRFS" };

    notran : bool = ---;
    count, i, j, kase, nz : int = ---;
    eps, ferr_j, lstres, s, safe1, safe2, safmin, xnorm : RT = ---;
    transn, transt : u8 = ---;

    info.* = 0;
    notran = lsame(trans, #char "N");
    if !notran && !lsame(trans, #char "T") && !lsame(trans, #char "C") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
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
        for j:1..nrhs {
            fset(ferr, j, zero_r);
            fset(berr, j, zero_r);
        }
        return;
    }

    if notran {
        transn = #char "N";
        transt = #char "C";
    } else {
        transn = #char "C";
        transt = #char "N";
    }

    nz = 4;
    eps = _gtrfs_lamch_epsilon(RT);
    safmin = _gtrfs_lamch_safe_minimum(RT);
    safe1 = ncast(RT, nz)*safmin;
    safe2 = safe1 / eps;

    active_view := _gtrfs_subvector(work.*, 1, n);
    v_view      := _gtrfs_subvector(work.*, n + 1, n);
    active_mat := array2d_view(active_view, max(1, n), row_major=false);

    for j:1..nrhs {
        b_view := _gtrfs_matrix_column(b, ldb, j, n);
        x_view := _gtrfs_matrix_column(x.*, ldx, j, n);
        x_mat := array2d_view(x_view, max(1, n), row_major=false);

        count = 1;
        lstres = three_r;

        while true {
            _gtrfs_copy(n, b_view, *active_view, RT, MT);
            _gtrfs_lagtm(trans, n, 1, neg_one_r, dl, d, du, *x_mat, n, one_r, *active_mat, n, RT, MT);

            if notran {
                if n == 1 {
                    fset(rwork, 1, _gtrfs_abs1(fget(b_view, 1), RT, MT) + _gtrfs_abs1(fget(d, 1), RT, MT)*_gtrfs_abs1(fget(x_view, 1), RT, MT));
                } else {
                    fset(rwork, 1, _gtrfs_abs1(fget(b_view, 1), RT, MT) +
                        _gtrfs_abs1(fget(d, 1), RT, MT)*_gtrfs_abs1(fget(x_view, 1), RT, MT) +
                        _gtrfs_abs1(fget(du, 1), RT, MT)*_gtrfs_abs1(fget(x_view, 2), RT, MT));
                    for i:2..n - 1 {
                        fset(rwork, i,
                            _gtrfs_abs1(fget(b_view, i), RT, MT) +
                            _gtrfs_abs1(fget(dl, i - 1), RT, MT)*_gtrfs_abs1(fget(x_view, i - 1), RT, MT) +
                            _gtrfs_abs1(fget(d, i), RT, MT)*_gtrfs_abs1(fget(x_view, i), RT, MT) +
                            _gtrfs_abs1(fget(du, i), RT, MT)*_gtrfs_abs1(fget(x_view, i + 1), RT, MT));
                    }
                    fset(rwork, n,
                        _gtrfs_abs1(fget(b_view, n), RT, MT) +
                        _gtrfs_abs1(fget(dl, n - 1), RT, MT)*_gtrfs_abs1(fget(x_view, n - 1), RT, MT) +
                        _gtrfs_abs1(fget(d, n), RT, MT)*_gtrfs_abs1(fget(x_view, n), RT, MT));
                }
            } else {
                if n == 1 {
                    fset(rwork, 1, _gtrfs_abs1(fget(b_view, 1), RT, MT) + _gtrfs_abs1(fget(d, 1), RT, MT)*_gtrfs_abs1(fget(x_view, 1), RT, MT));
                } else {
                    fset(rwork, 1, _gtrfs_abs1(fget(b_view, 1), RT, MT) +
                        _gtrfs_abs1(fget(d, 1), RT, MT)*_gtrfs_abs1(fget(x_view, 1), RT, MT) +
                        _gtrfs_abs1(fget(dl, 1), RT, MT)*_gtrfs_abs1(fget(x_view, 2), RT, MT));
                    for i:2..n - 1 {
                        fset(rwork, i,
                            _gtrfs_abs1(fget(b_view, i), RT, MT) +
                            _gtrfs_abs1(fget(du, i - 1), RT, MT)*_gtrfs_abs1(fget(x_view, i - 1), RT, MT) +
                            _gtrfs_abs1(fget(d, i), RT, MT)*_gtrfs_abs1(fget(x_view, i), RT, MT) +
                            _gtrfs_abs1(fget(dl, i), RT, MT)*_gtrfs_abs1(fget(x_view, i + 1), RT, MT));
                    }
                    fset(rwork, n,
                        _gtrfs_abs1(fget(b_view, n), RT, MT) +
                        _gtrfs_abs1(fget(du, n - 1), RT, MT)*_gtrfs_abs1(fget(x_view, n - 1), RT, MT) +
                        _gtrfs_abs1(fget(d, n), RT, MT)*_gtrfs_abs1(fget(x_view, n), RT, MT));
                }
            }

            s = zero_r;
            for i:1..n {
                if fget(rwork.*, i) > safe2 {
                    s = max(s, _gtrfs_abs1(fget(active_view, i), RT, MT) / fget(rwork.*, i));
                } else {
                    s = max(s, (_gtrfs_abs1(fget(active_view, i), RT, MT) + safe1) / (fget(rwork.*, i) + safe1));
                }
            }
            fset(berr, j, s);

            if fget(berr.*, j) > eps && two_r*fget(berr.*, j) <= lstres && count <= itmax {
                _gtrfs_gttrs(trans, n, 1, dlf, df, duf, du2, ipiv, *active_mat, n, info, RT, MT);
                _gtrfs_axpy(n, one_t, active_view, *x_view, RT, MT);
                lstres = fget(berr.*, j);
                count += 1;
                continue;
            }

            for i:1..n {
                if fget(rwork.*, i) > safe2 {
                    fset(rwork, i, _gtrfs_abs1(fget(active_view, i), RT, MT) + ncast(RT, nz)*eps*fget(rwork.*, i));
                } else {
                    fset(rwork, i, _gtrfs_abs1(fget(active_view, i), RT, MT) + ncast(RT, nz)*eps*fget(rwork.*, i) + safe1);
                }
            }

            ferr_j = zero_r;
            kase = 0;
            isave := array1d(3, 0, 0, 0);

            while true {
                _gtrfs_lacn2(n, *v_view, *active_view, rwork, *ferr_j, *kase, *isave, RT, MT);
                if kase == 0 break;

                if kase == 1 {
                    _gtrfs_gttrs(transt, n, 1, dlf, df, duf, du2, ipiv, *active_mat, n, info, RT, MT);
                    for i:1..n {
                        fset(active_view, i, fget(rwork.*, i) * fget(active_view, i));
                    }
                } else {
                    for i:1..n {
                        fset(active_view, i, fget(rwork.*, i) * fget(active_view, i));
                    }
                    _gtrfs_gttrs(transn, n, 1, dlf, df, duf, du2, ipiv, *active_mat, n, info, RT, MT);
                }
            }

            xnorm = zero_r;
            for i:1..n {
                xnorm = max(xnorm, _gtrfs_abs1(fget(x_view, i), RT, MT));
            }
            if xnorm != zero_r {
                ferr_j /= xnorm;
            }
            fset(ferr, j, ferr_j);
            break;
        }
    }
}
