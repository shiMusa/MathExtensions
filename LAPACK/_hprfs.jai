// ? LAPACK *hprfs()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chprfs/zhprfs
// *HPRFS improves the computed solution to a system of linear
// equations when the coefficient matrix is Hermitian indefinite
// and packed.

chprfs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, ap: $AP/ArrayType1D, afp: $AFP/ArrayType1D,
    ipiv: $IPIV/ArrayType1D, b: $B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D, work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D, info: *int
) {
    _hprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float32, MT=C32);
}

zhprfs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, ap: $AP/ArrayType1D, afp: $AFP/ArrayType1D,
    ipiv: $IPIV/ArrayType1D, b: $B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D, work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D, info: *int
) {
    _hprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_hprfs_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_hprfs_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
    view: [] A.data_type;
    view.data = a.data + (col - 1)*ld;
    view.count = count;
    return view;
}

_hprfs_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_hprfs_lamch :: ($$cmach: string, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slamch(cmach.data); } else { return dlamch(cmach.data); }
}

_hprfs_copy :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { ccopy(n, x, incx, y, incy); } else { zcopy(n, x, incx, y, incy); }
}

_hprfs_axpy :: (
    $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { caxpy(n, alpha, x, incx, y, incy); } else { zaxpy(n, alpha, x, incx, y, incy); }
}

_hprfs_hpmv :: (
    $$uplo: u8, $$n: int, alpha: $MT, ap: $AP/ArrayType1D, x: $X/ArrayType1D, $$incx: int,
    beta: MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy); } else { zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy); }
}

_hprfs_hptrs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, afp: $AFP/ArrayType1D, ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AFP.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    #if MT == C32 { chptrs(uplo, n, nrhs, afp, ipiv, b, ldb, info); } else { zhptrs(uplo, n, nrhs, afp, ipiv, b, ldb, info); }
}

_hprfs_lacn2 :: (
    $$n: int, v: *$V/ArrayType1D, x: *$X/ArrayType1D, est: *$RT, kase: *int,
    isave: *$ISAVE/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           V.data_type == MT && X.data_type == MT && ISAVE.data_type == int;
} {
    #if MT == C32 { clacn2(n, v, x, est, kase, isave); } else { zlacn2(n, v, x, est, kase, isave); }
}

_hprfs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, ap: $AP/ArrayType1D, afp: $AFP/ArrayType1D,
    ipiv: $IPIV/ArrayType1D, b: $B/ArrayType2D, $$ldb: int, x: *$X/ArrayType2D, $$ldx: int,
    ferr: *$FERR/ArrayType1D, berr: *$BERR/ArrayType1D, work: *$WORK/ArrayType1D,
    rwork: *$RWORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && AFP.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && X.data_type == MT &&
           FERR.data_type == RT && BERR.data_type == RT &&
           WORK.data_type == MT && RWORK.data_type == RT;
} {
    itmax :: 5;
    zero_r :: #run zero(RT);
    one_t  :: #run one(MT);
    two_r  :: #run ncast(RT, 2);
    three_r :: #run ncast(RT, 3);
    neg_one_t :: #run -one(MT);
    routine_name :: #run if MT == C32 { "CHPRFS" } else { "ZHPRFS" };

    count, i, ik, j, k, kase, kk, nz : int = ---;
    berr_j, eps, ferr_j, lstres, s, safe1, safe2, safmin, xk, xnorm : RT = ---;
    upper : bool = ---;
    isave : [3] int;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !upper && !lsame(uplo, #char "L") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if ldb < max(1, n) {
        info.* = -8;
    } else if ldx < max(1, n) {
        info.* = -10;
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

    nz = n + 1;
    eps = _hprfs_lamch("Epsilon", RT);
    safmin = _hprfs_lamch("Safe minimum", RT);
    safe1 = ncast(RT, nz) * safmin;
    safe2 = safe1 / eps;

    work_res := _hprfs_subvector(work.*, 1, n);
    work_v := _hprfs_subvector(work.*, n + 1, n);
    work_mat := array2d_view(work_res, max(1, n), row_major=false);

    for j:1..nrhs {
        b_view := _hprfs_matrix_column(b, ldb, j, n);
        x_view := _hprfs_matrix_column(x.*, ldx, j, n);

        count = 1;
        lstres = three_r;

        while true {
            _hprfs_copy(n, b_view, 1, *work_res, 1, MT);
            _hprfs_hpmv(uplo, n, neg_one_t, ap, x_view, 1, one_t, *work_res, 1, MT);

            for i:1..n {
                fset(rwork, i, _hprfs_abs1(fget(b_view, i), RT, MT));
            }

            kk = 1;
            if upper {
                for k:1..n {
                    s = zero_r;
                    xk = _hprfs_abs1(fget(x_view, k), RT, MT);
                    ik = kk;
                    for i:1..k - 1 {
                        fset(rwork, i, fget(rwork, i) + _hprfs_abs1(fget(ap, ik), RT, MT) * xk);
                        s += _hprfs_abs1(fget(ap, ik), RT, MT) * _hprfs_abs1(fget(x_view, i), RT, MT);
                        ik += 1;
                    }
                    fset(rwork, k, fget(rwork, k) + abs(real(fget(ap, kk + k - 1))) * xk + s);
                    kk += k;
                }
            } else {
                for k:1..n {
                    s = zero_r;
                    xk = _hprfs_abs1(fget(x_view, k), RT, MT);
                    fset(rwork, k, fget(rwork, k) + abs(real(fget(ap, kk))) * xk);
                    ik = kk + 1;
                    for i:k + 1..n {
                        fset(rwork, i, fget(rwork, i) + _hprfs_abs1(fget(ap, ik), RT, MT) * xk);
                        s += _hprfs_abs1(fget(ap, ik), RT, MT) * _hprfs_abs1(fget(x_view, i), RT, MT);
                        ik += 1;
                    }
                    fset(rwork, k, fget(rwork, k) + s);
                    kk += n - k + 1;
                }
            }

            s = zero_r;
            for i:1..n {
                if fget(rwork, i) > safe2 {
                    s = max(s, _hprfs_abs1(fget(work_res, i), RT, MT) / fget(rwork, i));
                } else {
                    s = max(s, (_hprfs_abs1(fget(work_res, i), RT, MT) + safe1) / (fget(rwork, i) + safe1));
                }
            }
            berr_j = s;
            fset(berr, j, berr_j);

            if berr_j > eps && two_r*berr_j <= lstres && count <= itmax {
                _hprfs_hptrs(uplo, n, 1, afp, ipiv, *work_mat, n, info, MT);
                _hprfs_axpy(n, one_t, work_res, 1, *x_view, 1, MT);
                lstres = berr_j;
                count += 1;
                continue;
            }

            for i:1..n {
                if fget(rwork, i) > safe2 {
                    fset(rwork, i, _hprfs_abs1(fget(work_res, i), RT, MT) + ncast(RT, nz)*eps*fget(rwork, i));
                } else {
                    fset(rwork, i, _hprfs_abs1(fget(work_res, i), RT, MT) + ncast(RT, nz)*eps*fget(rwork, i) + safe1);
                }
            }

            ferr_j = zero_r;
            kase = 0;
            isave = .[0, 0, 0];

            while true {
                _hprfs_lacn2(n, *work_v, *work_res, *ferr_j, *kase, *isave, RT, MT);
                if kase == 0 then break;

                if kase == 1 {
                    _hprfs_hptrs(uplo, n, 1, afp, ipiv, *work_mat, n, info, MT);
                    for i:1..n {
                        fset(work_res, i, fget(rwork, i) * fget(work_res, i));
                    }
                } else if kase == 2 {
                    for i:1..n {
                        fset(work_res, i, fget(rwork, i) * fget(work_res, i));
                    }
                    _hprfs_hptrs(uplo, n, 1, afp, ipiv, *work_mat, n, info, MT);
                }
            }

            xnorm = zero_r;
            for i:1..n {
                xnorm = max(xnorm, _hprfs_abs1(fget(x_view, i), RT, MT));
            }
            if xnorm != zero_r {
                ferr_j /= xnorm;
            }
            fset(ferr, j, ferr_j);
            break;
        }
    }
}
