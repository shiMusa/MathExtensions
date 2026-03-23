// ? LAPACK *herfs()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment cherfs/zherfs
// *HERFS improves the computed solution to a system of linear
// equations when the coefficient matrix is Hermitian indefinite.

cherfs :: (
    $$uplo: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
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
    _herfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float32, MT=C32);
}

zherfs :: (
    $$uplo: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
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
    _herfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_herfs_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_herfs_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
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

_herfs_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_herfs_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch("Epsilon".data);
    } else {
        return dlamch("Epsilon".data);
    }
}

_herfs_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch("Safe minimum".data);
    } else {
        return dlamch("Safe minimum".data);
    }
}

_herfs_lacn2 :: (
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
    #if MT == C32 {
        clacn2(n, v, x, est, kase, isave);
    } else {
        zlacn2(n, v, x, est, kase, isave);
    }
}

_herfs_hetrs :: (
    $$uplo: u8,
    $$n: int,
    $$nrhs: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT;
} {
    #if MT == C32 {
        chetrs(uplo, n, nrhs, af, ldaf, ipiv, b, ldb, info);
    } else {
        zhetrs(uplo, n, nrhs, af, ldaf, ipiv, b, ldb, info);
    }
}

_herfs_copy :: (
    $$n: int,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        ccopy(n, x, 1, y, 1);
    } else {
        zcopy(n, x, 1, y, 1);
    }
}

_herfs_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        caxpy(n, alpha, x, 1, y, 1);
    } else {
        zaxpy(n, alpha, x, 1, y, 1);
    }
}

_herfs_hemv :: (
    $$uplo: u8,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    x: $X/ArrayType1D,
    beta: $MT,
    y: *$Y/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == C32 {
        chemv(uplo, n, alpha, a, lda, x, 1, beta, y, 1);
    } else {
        zhemv(uplo, n, alpha, a, lda, x, 1, beta, y, 1);
    }
}

_herfs :: (
    $$uplo: u8,
    $$n: int,
    $$nrhs: int,
    a: $A/ArrayType2D,
    $$lda: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
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
           A.data_type == MT &&
           AF.data_type == MT &&
           IPIV.data_type == int &&
           B.data_type == MT &&
           X.data_type == MT &&
           FERR.data_type == RT &&
           BERR.data_type == RT &&
           WORK.data_type == MT &&
           RWORK.data_type == RT;
} {
    itmax :: 5;
    zero_r :: #run zero(RT);
    one_t  :: #run one(MT);
    two_r  :: #run ncast(RT, 2);
    three_r :: #run ncast(RT, 3);
    neg_one_t :: #run -one(MT);
    routine_name :: #run if MT == C32 { "CHERFS" } else { "ZHERFS" };

    count, i, j, k, kase, nz : int = ---;
    berr_j, eps, ferr_j, lstres, s, safe1, safe2, safmin, xk, xnorm : RT = ---;
    upper : bool = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !upper && !lsame(uplo, #char "L") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ldaf < max(1, n) {
        info.* = -7;
    } else if ldb < max(1, n) {
        info.* = -10;
    } else if ldx < max(1, n) {
        info.* = -12;
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
    eps = _herfs_lamch_epsilon(RT);
    safmin = _herfs_lamch_safe_minimum(RT);
    safe1 = ncast(RT, nz) * safmin;
    safe2 = safe1 / eps;

    work_res := _herfs_subvector(work.*, 1, n);
    work_v := _herfs_subvector(work.*, n + 1, n);
    work_mat := array2d_view(work_res, max(1, n), row_major=false);

    for j:1..nrhs {
        b_view := _herfs_matrix_column(b, ldb, j, n);
        x_view := _herfs_matrix_column(x.*, ldx, j, n);

        count = 1;
        lstres = three_r;

        while true {
            _herfs_copy(n, b_view, *work_res, MT);
            _herfs_hemv(uplo, n, neg_one_t, a, lda, x_view, one_t, *work_res, MT);

            for i:1..n {
                fset(rwork, i, _herfs_abs1(fget(b_view, i), RT, MT));
            }

            if upper {
                for k:1..n {
                    s = zero_r;
                    xk = _herfs_abs1(fget(x_view, k), RT, MT);
                    for i:1..k - 1 {
                        fset(rwork, i, fget(rwork, i) + _herfs_abs1(fget(a, i, k), RT, MT) * xk);
                        s += _herfs_abs1(fget(a, i, k), RT, MT) * _herfs_abs1(fget(x_view, i), RT, MT);
                    }
                    fset(rwork, k, fget(rwork, k) + abs(real(fget(a, k, k))) * xk + s);
                }
            } else {
                for k:1..n {
                    s = zero_r;
                    xk = _herfs_abs1(fget(x_view, k), RT, MT);
                    fset(rwork, k, fget(rwork, k) + abs(real(fget(a, k, k))) * xk);
                    for i:k + 1..n {
                        fset(rwork, i, fget(rwork, i) + _herfs_abs1(fget(a, i, k), RT, MT) * xk);
                        s += _herfs_abs1(fget(a, i, k), RT, MT) * _herfs_abs1(fget(x_view, i), RT, MT);
                    }
                    fset(rwork, k, fget(rwork, k) + s);
                }
            }

            s = zero_r;
            for i:1..n {
                if fget(rwork, i) > safe2 {
                    s = max(s, _herfs_abs1(fget(work_res, i), RT, MT) / fget(rwork, i));
                } else {
                    s = max(s, (_herfs_abs1(fget(work_res, i), RT, MT) + safe1) / (fget(rwork, i) + safe1));
                }
            }
            berr_j = s;
            fset(berr, j, berr_j);

            if berr_j > eps && two_r*berr_j <= lstres && count <= itmax {
                _herfs_hetrs(uplo, n, 1, af, ldaf, ipiv, *work_mat, n, info, MT);
                _herfs_axpy(n, one_t, work_res, *x_view, MT);
                lstres = berr_j;
                count += 1;
                continue;
            }

            for i:1..n {
                if fget(rwork, i) > safe2 {
                    fset(rwork, i, _herfs_abs1(fget(work_res, i), RT, MT) + ncast(RT, nz)*eps*fget(rwork, i));
                } else {
                    fset(rwork, i, _herfs_abs1(fget(work_res, i), RT, MT) + ncast(RT, nz)*eps*fget(rwork, i) + safe1);
                }
            }

            ferr_j = zero_r;
            kase = 0;
            isave := array1d(3, 0, 0, 0);

            while true {
                _herfs_lacn2(n, *work_v, *work_res, *ferr_j, *kase, *isave, RT, MT);
                if kase == 0 then break;

                if kase == 1 {
                    _herfs_hetrs(uplo, n, 1, af, ldaf, ipiv, *work_mat, n, info, MT);
                    for i:1..n {
                        fset(work_res, i, fget(rwork, i) * fget(work_res, i));
                    }
                } else {
                    for i:1..n {
                        fset(work_res, i, fget(rwork, i) * fget(work_res, i));
                    }
                    _herfs_hetrs(uplo, n, 1, af, ldaf, ipiv, *work_mat, n, info, MT);
                }
            }

            xnorm = zero_r;
            for i:1..n {
                xnorm = max(xnorm, _herfs_abs1(fget(x_view, i), RT, MT));
            }
            if xnorm != zero_r {
                ferr_j /= xnorm;
            }
            fset(ferr, j, ferr_j);
            break;
        }
    }
}
