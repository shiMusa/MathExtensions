// ? LAPACK *gerfs()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgerfs/cgerfs/dgerfs/zgerfs
// *GERFS improves the computed solution to a system of linear
// equations and provides error bounds and backward error estimates
// for the solution.

sgerfs :: (
    $$trans: u8,
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
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, RT=float32, MT=float32);
}

dgerfs :: (
    $$trans: u8,
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
    iwork: *$IWORK/ArrayType1D,
    info: *int
) {
    _gerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, RT=float64, MT=float64);
}

cgerfs :: (
    $$trans: u8,
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
    _gerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float32, MT=C32);
}

zgerfs :: (
    $$trans: u8,
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
    _gerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gerfs_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gerfs_matrix_column :: inline (a: $A/ArrayType2D, $$ld: int, $$col: int, $$count: int) -> [] A.data_type {
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

_gerfs_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_gerfs_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gerfs_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_gerfs_lacn2 :: (
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

_gerfs_getrs :: (
    $$trans: u8,
    $$n: int,
    af: $AF/ArrayType2D,
    $$ldaf: int,
    ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType1D,
    $$ldb: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           AF.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    #if MT == float32 {
        sgetrs(trans, n, 1, af, ldaf, ipiv, b, ldb, info);
    } else if MT == float64 {
        dgetrs(trans, n, 1, af, ldaf, ipiv, b, ldb, info);
    } else if MT == C32 {
        cgetrs(trans, n, 1, af, ldaf, ipiv, b, ldb, info);
    } else {
        zgetrs(trans, n, 1, af, ldaf, ipiv, b, ldb, info);
    }
}

_gerfs_copy :: (
    $$n: int,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
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

_gerfs_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    y: *$Y/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
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

_gerfs_gemv :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    x: $X/ArrayType1D,
    beta: $MT,
    y: *$Y/ArrayType1D,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sgemv(trans, m, n, alpha, a, lda, x, 1, beta, y, 1);
    } else if MT == float64 {
        dgemv(trans, m, n, alpha, a, lda, x, 1, beta, y, 1);
    } else if MT == C32 {
        cgemv(trans, m, n, alpha, a, lda, x, 1, beta, y, 1);
    } else {
        zgemv(trans, m, n, alpha, a, lda, x, 1, beta, y, 1);
    }
}

_gerfs :: (
    $$trans: u8,
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
            assert(A.data_type    == MT);
            assert(AF.data_type   == MT);
            assert(IPIV.data_type == int);
            assert(B.data_type    == MT);
            assert(X.data_type    == MT);
            assert(FERR.data_type == RT);
            assert(BERR.data_type == RT);
            assert(WORK.data_type == MT);
        }
        #if MT == RT {
            #run assert(AUX.data_type == int);
        } else {
            #run assert(AUX.data_type == RT);
        }
    }

    itmax :: 5;
    zero_r :: #run zero(RT);
    one_t  :: #run one(MT);
    two_r  :: #run ncast(RT, 2);
    three_r :: #run ncast(RT, 3);
    neg_one_t :: #run -one(MT);
    routine_name :: #run if MT == float32 {
        "SGERFS"
    } else if MT == float64 {
        "DGERFS"
    } else if MT == C32 {
        "CGERFS"
    } else {
        "ZGERFS"
    };

    count, i, j, k, kase, nz : int = ---;
    berr_j, eps, ferr_j, lstres, s, safe1, safe2, safmin, xk, xnorm : RT = ---;
    notrans : bool = ---;
    transt, transn : u8 = ---;

    info.* = 0;
    notrans = char_eq(trans, #char "n");
    if !notrans && !char_eq(trans, #char "t") && !char_eq(trans, #char "c") {
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

    #if MT == float32 || MT == float64 {
        transn = trans;
        transt = ifx notrans then #char "t" else #char "n";
    } else {
        transn = ifx notrans then #char "n" else #char "c";
        transt = ifx notrans then #char "c" else #char "n";
    }

    nz = n + 1;
    eps = _gerfs_lamch_epsilon(RT);
    safmin = _gerfs_lamch_safe_minimum(RT);
    safe1 = ncast(RT, nz) * safmin;
    safe2 = safe1 / eps;

    weights_view : [] RT = ---;
    active_view  : [] MT = ---;
    v_view       : [] MT = ---;
    #if MT == float32 || MT == float64 {
        weights_view = _gerfs_subvector(work.*, 1, n);
        active_view  = _gerfs_subvector(work.*, n + 1, n);
        v_view       = _gerfs_subvector(work.*, 2*n + 1, n);
    } else {
        weights_view = _gerfs_subvector(aux.*, 1, n);
        active_view  = _gerfs_subvector(work.*, 1, n);
        v_view       = _gerfs_subvector(work.*, n + 1, n);
    }

    for j:1..nrhs {
        b_view := _gerfs_matrix_column(b, ldb, j, n);
        x_view := _gerfs_matrix_column(x.*, ldx, j, n);

        count = 1;
        lstres = three_r;

        while true {
            _gerfs_copy(n, b_view, *active_view, RT, MT);
            _gerfs_gemv(trans, n, n, neg_one_t, a, lda, x_view, one_t, *active_view, RT, MT);

            for i:1..n {
                fset(weights_view, i, _gerfs_abs1(fget(b_view, i), RT, MT));
            }

            if notrans {
                for k:1..n {
                    xk = _gerfs_abs1(fget(x_view, k), RT, MT);
                    for i:1..n {
                        fset(weights_view, i, fget(weights_view, i) + _gerfs_abs1(fget(a, i, k), RT, MT)*xk);
                    }
                }
            } else {
                for k:1..n {
                    s = zero_r;
                    for i:1..n {
                        s += _gerfs_abs1(fget(a, i, k), RT, MT) * _gerfs_abs1(fget(x_view, i), RT, MT);
                    }
                    fset(weights_view, k, fget(weights_view, k) + s);
                }
            }

            s = zero_r;
            for i:1..n {
                if fget(weights_view, i) > safe2 {
                    s = max(s, _gerfs_abs1(fget(active_view, i), RT, MT) / fget(weights_view, i));
                } else {
                    s = max(s, (_gerfs_abs1(fget(active_view, i), RT, MT) + safe1) / (fget(weights_view, i) + safe1));
                }
            }
            berr_j = s;
            fset(berr, j, berr_j);

            if berr_j > eps && two_r*berr_j <= lstres && count <= itmax {
                _gerfs_getrs(trans, n, af, ldaf, ipiv, *active_view, n, info, RT, MT);
                _gerfs_axpy(n, one_t, active_view, *x_view, RT, MT);
                lstres = berr_j;
                count += 1;
                continue;
            }

            for i:1..n {
                if fget(weights_view, i) > safe2 {
                    fset(weights_view, i, _gerfs_abs1(fget(active_view, i), RT, MT) + ncast(RT, nz)*eps*fget(weights_view, i));
                } else {
                    fset(weights_view, i, _gerfs_abs1(fget(active_view, i), RT, MT) + ncast(RT, nz)*eps*fget(weights_view, i) + safe1);
                }
            }

            ferr_j = zero_r;
            kase = 0;
            isave := array1d(3, 0, 0, 0);

            while true {
                _gerfs_lacn2(n, *v_view, *active_view, aux, *ferr_j, *kase, *isave, RT, MT);
                if kase == 0 then break;

                if kase == 1 {
                    _gerfs_getrs(transt, n, af, ldaf, ipiv, *active_view, n, info, RT, MT);
                    for i:1..n {
                        fset(active_view, i, fget(weights_view, i) * fget(active_view, i));
                    }
                } else {
                    for i:1..n {
                        fset(active_view, i, fget(weights_view, i) * fget(active_view, i));
                    }
                    _gerfs_getrs(transn, n, af, ldaf, ipiv, *active_view, n, info, RT, MT);
                }
            }

            xnorm = zero_r;
            for i:1..n {
                xnorm = max(xnorm, _gerfs_abs1(fget(x_view, i), RT, MT));
            }
            if xnorm != zero_r {
                ferr_j /= xnorm;
            }
            fset(ferr, j, ferr_j);
            break;
        }
    }
}
