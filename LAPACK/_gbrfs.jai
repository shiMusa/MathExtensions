// ? LAPACK *gbrfs()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgbrfs/cgbrfs/dgbrfs/zgbrfs
// *GBRFS improves the computed solution to a system of linear
// equations when the coefficient matrix is banded, and provides
// error bounds and backward error estimates for the solution.

sgbrfs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
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
    _gbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, RT=float32, MT=float32);
}

dgbrfs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
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
    _gbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, RT=float64, MT=float64);
}

cgbrfs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
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
    _gbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float32, MT=C32);
}

zgbrfs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
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
    _gbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info, RT=float64, MT=C64);
}

#scope_file
_gbrfs_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    #if CHECKS {
        assert(offset >= 1);
        assert(count >= 0);
    }

    view: [] X.data_type = x;
    view.data += offset - 1;
    view.count = count;
    return view;
}

_gbrfs_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "e");
    } else {
        return dlamch(#char "e");
    }
}

_gbrfs_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_gbrfs_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return abs(z);
    } else {
        return _abs1(z);
    }
}

_gbrfs_lacn2 :: (
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

_gbrfs_gbtrs :: (
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
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
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

_gbrfs :: (
    $$trans: u8,
    $$n: int,
    $$kl: int,
    $$ku: int,
    $$nrhs: int,
    ab: $AB/ArrayType2D,
    $$ldab: int,
    afb: $AFB/ArrayType2D,
    $$ldafb: int,
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
            assert(AB.data_type   == MT);
            assert(AFB.data_type  == MT);
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
    one_r  :: #run one(RT);
    two_r  :: #run ncast(RT, 2);
    three_r :: #run ncast(RT, 3);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == float32 {
        "SGBRFS"
    } else if MT == float64 {
        "DGBRFS"
    } else if MT == C32 {
        "CGBRFS"
    } else {
        "ZGBRFS"
    };

    count, i, j, k, kase, kk, nz : int = ---;
    berr_j, eps, ferr_j, lstres, s, safe1, safe2, safmin, xk, xnorm : RT = ---;
    notrans, trans_no_conj : bool = ---;
    transt, transn : u8 = ---;
    sum_t : MT = ---;

    info.* = 0;
    notrans = char_eq(trans, #char "n");
    if !notrans && !char_eq(trans, #char "t") && !char_eq(trans, #char "c") {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if kl < 0 {
        info.* = -3;
    } else if ku < 0 {
        info.* = -4;
    } else if nrhs < 0 {
        info.* = -5;
    } else if ldab < kl + ku + 1 {
        info.* = -7;
    } else if ldafb < 2*kl + ku + 1 {
        info.* = -9;
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

    nz = min(kl + ku + 2, n + 1);
    eps = _gbrfs_lamch_epsilon(RT);
    safmin = _gbrfs_lamch_safe_minimum(RT);
    safe1 = ncast(RT, nz) * safmin;
    safe2 = safe1 / eps;

    weights_view : [] RT = ---;
    active_view  : [] MT = ---;
    v_view       : [] MT = ---;
    #if MT == float32 || MT == float64 {
        weights_view = _gbrfs_subvector(work.*, 1, n);
        active_view  = _gbrfs_subvector(work.*, n + 1, n);
        v_view       = _gbrfs_subvector(work.*, 2*n + 1, n);
    } else {
        weights_view = _gbrfs_subvector(aux.*, 1, n);
        active_view  = _gbrfs_subvector(work.*, 1, n);
        v_view       = _gbrfs_subvector(work.*, n + 1, n);
    }

    trans_no_conj = char_eq(trans, #char "t");

    for j:1..nrhs {
        count = 1;
        lstres = three_r;

        while true {
            for i:1..n {
                fset(active_view, i, fget(b, i, j));
            }

            if notrans {
                for k:1..n {
                    kk = ku + 1 - k;
                    for i:max(1, k-ku)..min(n, k+kl) {
                        fset(active_view, i, fget(active_view, i) - fget(ab, kk+i, k)*fget(x.*, k, j));
                    }
                }
            } else {
                for k:1..n {
                    kk = ku + 1 - k;
                    sum_t = zero_t;
                    #if MT == float32 || MT == float64 {
                        for i:max(1, k-ku)..min(n, k+kl) {
                            sum_t += fget(ab, kk+i, k) * fget(x.*, i, j);
                        }
                    } else {
                        if trans_no_conj {
                            for i:max(1, k-ku)..min(n, k+kl) {
                                sum_t += fget(ab, kk+i, k) * fget(x.*, i, j);
                            }
                        } else {
                            for i:max(1, k-ku)..min(n, k+kl) {
                                sum_t += conjugate(fget(ab, kk+i, k)) * fget(x.*, i, j);
                            }
                        }
                    }
                    fset(active_view, k, fget(active_view, k) - sum_t);
                }
            }

            for i:1..n {
                fset(weights_view, i, _gbrfs_abs1(fget(b, i, j), RT, MT));
            }

            if notrans {
                for k:1..n {
                    kk = ku + 1 - k;
                    xk = _gbrfs_abs1(fget(x.*, k, j), RT, MT);
                    for i:max(1, k-ku)..min(n, k+kl) {
                        fset(weights_view, i, fget(weights_view, i) + _gbrfs_abs1(fget(ab, kk+i, k), RT, MT)*xk);
                    }
                }
            } else {
                for k:1..n {
                    kk = ku + 1 - k;
                    s = zero_r;
                    for i:max(1, k-ku)..min(n, k+kl) {
                        s += _gbrfs_abs1(fget(ab, kk+i, k), RT, MT) * _gbrfs_abs1(fget(x.*, i, j), RT, MT);
                    }
                    fset(weights_view, k, fget(weights_view, k) + s);
                }
            }

            s = zero_r;
            for i:1..n {
                if fget(weights_view, i) > safe2 {
                    s = max(s, _gbrfs_abs1(fget(active_view, i), RT, MT) / fget(weights_view, i));
                } else {
                    s = max(s, (_gbrfs_abs1(fget(active_view, i), RT, MT) + safe1) / (fget(weights_view, i) + safe1));
                }
            }
            berr_j = s;
            fset(berr, j, berr_j);

            if berr_j > eps && two_r*berr_j <= lstres && count <= itmax {
                _gbrfs_gbtrs(trans, n, kl, ku, afb, ldafb, ipiv, *active_view, n, info, RT, MT);
                for i:1..n {
                    fset(x, i, j, fget(x.*, i, j) + fget(active_view, i));
                }
                lstres = berr_j;
                count += 1;
                continue;
            }

            for i:1..n {
                if fget(weights_view, i) > safe2 {
                    fset(weights_view, i, _gbrfs_abs1(fget(active_view, i), RT, MT) + ncast(RT, nz)*eps*fget(weights_view, i));
                } else {
                    fset(weights_view, i, _gbrfs_abs1(fget(active_view, i), RT, MT) + ncast(RT, nz)*eps*fget(weights_view, i) + safe1);
                }
            }

            ferr_j = zero_r;
            kase = 0;
            isave := array1d(3, 0, 0, 0);

            while true {
                _gbrfs_lacn2(n, *v_view, *active_view, aux, *ferr_j, *kase, *isave, RT, MT);
                if kase == 0 then break;

                if kase == 1 {
                    _gbrfs_gbtrs(transt, n, kl, ku, afb, ldafb, ipiv, *active_view, n, info, RT, MT);
                    for i:1..n {
                        fset(active_view, i, fget(weights_view, i) * fget(active_view, i));
                    }
                } else {
                    for i:1..n {
                        fset(active_view, i, fget(weights_view, i) * fget(active_view, i));
                    }
                    _gbrfs_gbtrs(transn, n, kl, ku, afb, ldafb, ipiv, *active_view, n, info, RT, MT);
                }
            }

            xnorm = zero_r;
            for i:1..n {
                xnorm = max(xnorm, _gbrfs_abs1(fget(x.*, i, j), RT, MT));
            }
            if xnorm != zero_r {
                ferr_j /= xnorm;
            }
            fset(ferr, j, ferr_j);
            break;
        }
    }
}
