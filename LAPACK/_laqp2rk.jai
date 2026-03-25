// ? LAPACK *laqp2rk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqp2rk/dlaqp2rk/claqp2rk/zlaqp2rk
// *LAQP2RK computes a truncated QR factorization with column pivoting
// of a matrix block using Level 2 BLAS.

slaqp2rk :: #bake_arguments _laqp2rk_real(RT=float32);
dlaqp2rk :: #bake_arguments _laqp2rk_real(RT=float64);
claqp2rk :: #bake_arguments _laqp2rk_complex(RT=float32, MT=C32);
zlaqp2rk :: #bake_arguments _laqp2rk_complex(RT=float64, MT=C64);

#scope_file
_laqp2rk_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_laqp2rk_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _laqp2rk_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqp2rk_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _laqp2rk_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqp2rk_iamax_1based :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return isamax(n, x, incx) + 1;
    } else {
        return idamax(n, x, incx) + 1;
    }
}

_laqp2rk_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "E");
    } else {
        return dlamch(#char "E");
    }
}

_laqp2rk_lamch_overflow :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "O");
    } else {
        return dlamch(#char "O");
    }
}

_laqp2rk_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_laqp2rk_nrm2 :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        return snrm2(n, x, incx);
    } else if MT == float64 {
        return dnrm2(n, x, incx);
    } else if MT == C32 {
        return scnrm2(n, x, incx);
    } else {
        return dznrm2(n, x, incx);
    }
}

_laqp2rk_swap :: ($$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $MT: Type) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sswap(n, x, incx, y, incy);
    } else if MT == float64 {
        dswap(n, x, incx, y, incy);
    } else if MT == C32 {
        cswap(n, x, incx, y, incy);
    } else {
        zswap(n, x, incx, y, incy);
    }
}

_laqp2rk_larfg :: ($$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *$TAU, $RT: Type, $MT: Type) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           X.data_type == MT &&
           TAU == MT;
} {
    #if MT == float32 {
        slarfg(n, alpha, x, incx, tau);
    } else if MT == float64 {
        dlarfg(n, alpha, x, incx, tau);
    } else if MT == C32 {
        clarfg(n, alpha, x, incx, tau);
    } else {
        zlarfg(n, alpha, x, incx, tau);
    }
}

_laqp2rk_larf1f :: (
    $$side: u8, $$m: int, $$n: int, v: $V/ArrayType1D, $$incv: int, $$tau: $MT,
    c: *$C/ArrayType2D, $$ldc: int, work: *$WORK/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        slarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == float64 {
        dlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == C32 {
        clarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else {
        zlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    }
}

_laqp2rk_apply_tau :: inline (tau: $MT, $MT: Type) -> MT #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    #if MT == float32 || MT == float64 {
        return tau;
    } else {
        return conjugate(tau);
    }
}

_laqp2rk_tau_nan_real :: (tau: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    if _laqp2rk_isnan(tau, RT) {
        return tau;
    }
    return zero(RT);
}

_laqp2rk_tau_nan_complex :: (tau: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    if _laqp2rk_isnan(real(tau), RT) {
        return real(tau);
    } else if _laqp2rk_isnan(imag(tau), RT) {
        return imag(tau);
    }
    return zero(RT);
}

_laqp2rk_real :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$kmax_: int, $$abstol: $RT, $$reltol: RT,
    $$kp1: int, $$maxc2nrm: RT, a: *$A/ArrayType2D, $$lda: int, k: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, work: *$WORK/ArrayType1D, info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPIV.data_type == int &&
           TAU.data_type == RT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           WORK.data_type == RT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);

    hugeval, taunan, temp, temp2, tol3z : RT = ---;
    i, itemp, j, jmaxc2nrm, kk, kp, minmnfact, minmnupdt, kmax_use : int = ---;

    info.* = 0;
    minmnfact = min(m - ioffset, n);
    minmnupdt = min(m - ioffset, n + nrhs);
    kmax_use = min(kmax_, minmnfact);
    tol3z = sqrt(_laqp2rk_lamch_epsilon(RT));
    hugeval = _laqp2rk_lamch_overflow(RT);

    for kk:1..kmax_use {
        i = ioffset + kk;

        if i == 1 {
            kp = kp1;
        } else {
            kp = (kk - 1) + _laqp2rk_iamax_1based(n - kk + 1, _laqp2rk_raw_slice(vn1.data + kk - 1, n - kk + 1), 1, RT);
            maxc2nrmk.* = fget(vn1.*, kp);

            if _laqp2rk_isnan(maxc2nrmk.*, RT) {
                k.* = kk - 1;
                info.* = k.* + kp;
                relmaxc2nrmk.* = maxc2nrmk.*;
                return;
            }

            if maxc2nrmk.* == zero_r {
                k.* = kk - 1;
                relmaxc2nrmk.* = zero_r;
                for j:kk..minmnfact {
                    fset(tau, j, zero_r);
                }
                return;
            }

            if info.* == 0 && maxc2nrmk.* > hugeval {
                info.* = n + kk - 1 + kp;
            }

            relmaxc2nrmk.* = maxc2nrmk.* / maxc2nrm;
            if maxc2nrmk.* <= abstol || relmaxc2nrmk.* <= reltol {
                k.* = kk - 1;
                for j:kk..minmnfact {
                    fset(tau, j, zero_r);
                }
                return;
            }
        }

        if kp != kk {
            _laqp2rk_swap(m, _laqp2rk_matrix_column(a.*, lda, 1, kp, m), 1, _laqp2rk_matrix_column(a.*, lda, 1, kk, m), 1, RT);
            fset(vn1, kp, fget(vn1.*, kk));
            fset(vn2, kp, fget(vn2.*, kk));
            itemp = fget(jpiv.*, kp);
            fset(jpiv, kp, fget(jpiv.*, kk));
            fset(jpiv, kk, itemp);
        }

        if i < m {
            _laqp2rk_larfg(m - i + 1, a.data + (kk - 1)*lda + (i - 1), *_laqp2rk_matrix_column(a.*, lda, i + 1, kk, m - i), 1, tau.data + (kk - 1), RT, RT);
        } else {
            _laqp2rk_larfg(1, a.data + (kk - 1)*lda + (m - 1), *_laqp2rk_matrix_column(a.*, lda, m, kk, 0), 1, tau.data + (kk - 1), RT, RT);
        }

        taunan = _laqp2rk_tau_nan_real(fget(tau.*, kk), RT);
        if _laqp2rk_isnan(taunan, RT) {
            k.* = kk - 1;
            info.* = kk;
            maxc2nrmk.* = taunan;
            relmaxc2nrmk.* = taunan;
            return;
        }

        if kk < minmnupdt {
            _laqp2rk_larf1f(#char "L", m - i + 1, n + nrhs - kk, _laqp2rk_matrix_column(a.*, lda, i, kk, m - i + 1), 1,
                            fget(tau.*, kk), *_laqp2rk_matrix_storage(a.*, lda, i, kk + 1, m - i + 1, n + nrhs - kk), lda, work, RT, RT);
        }

        if kk < minmnfact {
            for j:kk + 1..n {
                if fget(vn1.*, j) != zero_r {
                    temp = one_r - (abs(fget(a.*, i, j)) / fget(vn1.*, j))*(abs(fget(a.*, i, j)) / fget(vn1.*, j));
                    temp = max(temp, zero_r);
                    temp2 = temp*(fget(vn1.*, j) / fget(vn2.*, j))*(fget(vn1.*, j) / fget(vn2.*, j));
                    if temp2 <= tol3z {
                        fset(vn1, j, _laqp2rk_nrm2(m - i, _laqp2rk_matrix_column(a.*, lda, i + 1, j, m - i), 1, RT, RT));
                        fset(vn2, j, fget(vn1.*, j));
                    } else {
                        fset(vn1, j, fget(vn1.*, j)*sqrt(temp));
                    }
                }
            }
        }
    }

    k.* = kmax_use;
    if k.* < minmnfact {
        jmaxc2nrm = k.* + _laqp2rk_iamax_1based(n - k.*, _laqp2rk_raw_slice(vn1.data + k.*, n - k.*), 1, RT);
        maxc2nrmk.* = fget(vn1.*, jmaxc2nrm);
        if k.* == 0 {
            relmaxc2nrmk.* = one_r;
        } else {
            relmaxc2nrmk.* = maxc2nrmk.* / maxc2nrm;
        }
    } else {
        maxc2nrmk.* = zero_r;
        relmaxc2nrmk.* = zero_r;
    }

    for j:k.* + 1..minmnfact {
        fset(tau, j, zero_r);
    }
}

_laqp2rk_complex :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$kmax_: int, $$abstol: $RT, $$reltol: RT,
    $$kp1: int, $$maxc2nrm: RT, a: *$A/ArrayType2D, $$lda: int, k: *int,
    maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D, tau: *$TAU/ArrayType1D,
    vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D, work: *$WORK/ArrayType1D, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           JPIV.data_type == int &&
           TAU.data_type == MT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           WORK.data_type == MT;
} {
    zero_r :: #run zero(RT);
    zero_t :: #run zero(MT);
    one_r  :: #run one(RT);

    hugeval, taunan, temp, temp2, tol3z : RT = ---;
    i, itemp, j, jmaxc2nrm, kk, kp, minmnfact, minmnupdt, kmax_use : int = ---;

    info.* = 0;
    minmnfact = min(m - ioffset, n);
    minmnupdt = min(m - ioffset, n + nrhs);
    kmax_use = min(kmax_, minmnfact);
    tol3z = sqrt(_laqp2rk_lamch_epsilon(RT));
    hugeval = _laqp2rk_lamch_overflow(RT);

    for kk:1..kmax_use {
        i = ioffset + kk;

        if i == 1 {
            kp = kp1;
        } else {
            kp = (kk - 1) + _laqp2rk_iamax_1based(n - kk + 1, _laqp2rk_raw_slice(vn1.data + kk - 1, n - kk + 1), 1, RT);
            maxc2nrmk.* = fget(vn1.*, kp);

            if _laqp2rk_isnan(maxc2nrmk.*, RT) {
                k.* = kk - 1;
                info.* = k.* + kp;
                relmaxc2nrmk.* = maxc2nrmk.*;
                return;
            }

            if maxc2nrmk.* == zero_r {
                k.* = kk - 1;
                relmaxc2nrmk.* = zero_r;
                for j:kk..minmnfact {
                    fset(tau, j, zero_t);
                }
                return;
            }

            if info.* == 0 && maxc2nrmk.* > hugeval {
                info.* = n + kk - 1 + kp;
            }

            relmaxc2nrmk.* = maxc2nrmk.* / maxc2nrm;
            if maxc2nrmk.* <= abstol || relmaxc2nrmk.* <= reltol {
                k.* = kk - 1;
                for j:kk..minmnfact {
                    fset(tau, j, zero_t);
                }
                return;
            }
        }

        if kp != kk {
            _laqp2rk_swap(m, _laqp2rk_matrix_column(a.*, lda, 1, kp, m), 1, _laqp2rk_matrix_column(a.*, lda, 1, kk, m), 1, MT);
            fset(vn1, kp, fget(vn1.*, kk));
            fset(vn2, kp, fget(vn2.*, kk));
            itemp = fget(jpiv.*, kp);
            fset(jpiv, kp, fget(jpiv.*, kk));
            fset(jpiv, kk, itemp);
        }

        if i < m {
            _laqp2rk_larfg(m - i + 1, a.data + (kk - 1)*lda + (i - 1), *_laqp2rk_matrix_column(a.*, lda, i + 1, kk, m - i), 1, tau.data + (kk - 1), RT, MT);
        } else {
            _laqp2rk_larfg(1, a.data + (kk - 1)*lda + (m - 1), *_laqp2rk_matrix_column(a.*, lda, m, kk, 0), 1, tau.data + (kk - 1), RT, MT);
        }

        taunan = _laqp2rk_tau_nan_complex(fget(tau.*, kk), RT, MT);
        if _laqp2rk_isnan(taunan, RT) {
            k.* = kk - 1;
            info.* = kk;
            maxc2nrmk.* = taunan;
            relmaxc2nrmk.* = taunan;
            return;
        }

        if kk < minmnupdt {
            _laqp2rk_larf1f(#char "L", m - i + 1, n + nrhs - kk, _laqp2rk_matrix_column(a.*, lda, i, kk, m - i + 1), 1,
                            _laqp2rk_apply_tau(fget(tau.*, kk), MT), *_laqp2rk_matrix_storage(a.*, lda, i, kk + 1, m - i + 1, n + nrhs - kk), lda, work, RT, MT);
        }

        if kk < minmnfact {
            for j:kk + 1..n {
                if fget(vn1.*, j) != zero_r {
                    temp = one_r - (abs(fget(a.*, i, j)) / fget(vn1.*, j))*(abs(fget(a.*, i, j)) / fget(vn1.*, j));
                    temp = max(temp, zero_r);
                    temp2 = temp*(fget(vn1.*, j) / fget(vn2.*, j))*(fget(vn1.*, j) / fget(vn2.*, j));
                    if temp2 <= tol3z {
                        fset(vn1, j, _laqp2rk_nrm2(m - i, _laqp2rk_matrix_column(a.*, lda, i + 1, j, m - i), 1, RT, MT));
                        fset(vn2, j, fget(vn1.*, j));
                    } else {
                        fset(vn1, j, fget(vn1.*, j)*sqrt(temp));
                    }
                }
            }
        }
    }

    k.* = kmax_use;
    if k.* < minmnfact {
        jmaxc2nrm = k.* + _laqp2rk_iamax_1based(n - k.*, _laqp2rk_raw_slice(vn1.data + k.*, n - k.*), 1, RT);
        maxc2nrmk.* = fget(vn1.*, jmaxc2nrm);
        if k.* == 0 {
            relmaxc2nrmk.* = one_r;
        } else {
            relmaxc2nrmk.* = maxc2nrmk.* / maxc2nrm;
        }
    } else {
        maxc2nrmk.* = zero_r;
        relmaxc2nrmk.* = zero_r;
    }

    for j:k.* + 1..minmnfact {
        fset(tau, j, zero_t);
    }
}
