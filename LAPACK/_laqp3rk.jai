// ? LAPACK *laqp3rk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slaqp3rk/dlaqp3rk/claqp3rk/zlaqp3rk
// *LAQP3RK computes a blocked step of truncated QR factorization with
// column pivoting.

slaqp3rk :: #bake_arguments _laqp3rk_real(RT=float32);
dlaqp3rk :: #bake_arguments _laqp3rk_real(RT=float64);
claqp3rk :: #bake_arguments _laqp3rk_complex(RT=float32, MT=C32);
zlaqp3rk :: #bake_arguments _laqp3rk_complex(RT=float64, MT=C64);

#scope_file
_laqp3rk_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_laqp3rk_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _laqp3rk_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqp3rk_matrix_row :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _laqp3rk_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_laqp3rk_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _laqp3rk_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_laqp3rk_iamax_1based :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type) -> int #modify {
    return (RT == float32 || RT == float64) && X.data_type == RT;
} {
    #if RT == float32 {
        return isamax(n, x, incx) + 1;
    } else {
        return idamax(n, x, incx) + 1;
    }
}

_laqp3rk_lamch_epsilon :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "E");
    } else {
        return dlamch(#char "E");
    }
}

_laqp3rk_lamch_overflow :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "O");
    } else {
        return dlamch(#char "O");
    }
}

_laqp3rk_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_laqp3rk_nrm2 :: ($$n: int, x: $X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type) -> RT #modify {
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

_laqp3rk_swap :: ($$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $MT: Type) #modify {
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

_laqp3rk_larfg :: ($$n: int, alpha: *$MT, x: *$X/ArrayType1D, $$incx: int, tau: *$TAU, $RT: Type, $MT: Type) #modify {
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

_laqp3rk_gemv :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    x: $X/ArrayType1D,
    $$incx: int,
    beta: $MT,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == float64 {
        dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == C32 {
        cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_laqp3rk_gemm :: (
    $$transa: u8,
    $$transb: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: $MT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
} {
    #if MT == float32 {
        sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == float64 {
        dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_laqp3rk_tau_nan_real :: (tau: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    if _laqp3rk_isnan(tau, RT) {
        return tau;
    }
    return zero(RT);
}

_laqp3rk_tau_nan_complex :: (tau: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT);
} {
    if _laqp3rk_isnan(real(tau), RT) {
        return real(tau);
    } else if _laqp3rk_isnan(imag(tau), RT) {
        return imag(tau);
    }
    return zero(RT);
}

_laqp3rk_real :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$nb_: int, $$abstol: $RT,
    $$reltol: RT, $$kp1: int, $$maxc2nrm: RT, a: *$A/ArrayType2D, $$lda: int,
    done: *bool, kb: *int, maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D,
    tau: *$TAU/ArrayType1D, vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D,
    auxv: *$AUXV/ArrayType1D, f: *$F/ArrayType2D, $$ldf: int, iwork: *$IWORK/ArrayType1D,
    info: *int, $RT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           A.data_type == RT &&
           JPIV.data_type == int &&
           TAU.data_type == RT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           AUXV.data_type == RT &&
           F.data_type == RT &&
           IWORK.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);

    aik, hugeval, maxc2_local, relmax_local, taunan, temp, temp2, tol3z : RT = ---;
    i, if_, itemp, j, k, kp, lsticc, minmnfact, minmnupdt, nb_use : int = ---;

    info.* = 0;
    minmnfact = min(m - ioffset, n);
    minmnupdt = min(m - ioffset, n + nrhs);
    nb_use = min(nb_, minmnfact);
    tol3z = sqrt(_laqp3rk_lamch_epsilon(RT));
    hugeval = _laqp3rk_lamch_overflow(RT);

    k = 0;
    i = ioffset;
    lsticc = 0;
    done.* = false;

    while k < nb_use && lsticc == 0 {
        k += 1;
        i = ioffset + k;

        if i == 1 {
            kp = kp1;
        } else {
            kp = (k - 1) + _laqp3rk_iamax_1based(n - k + 1, _laqp3rk_raw_slice(vn1.data + k - 1, n - k + 1), 1, RT);
            maxc2_local = fget(vn1.*, kp);
            maxc2nrmk.* = maxc2_local;

            if _laqp3rk_isnan(maxc2_local, RT) {
                done.* = true;
                kb.* = k - 1;
                if_ = i - 1;
                info.* = kb.* + kp;
                relmaxc2nrmk.* = maxc2_local;

                if nrhs > 0 && kb.* < (m - ioffset) {
                    a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
                    a_tail := array2d_view(a_tail_raw, lda, row_major=false);
                    f_rhs_raw := _laqp3rk_matrix_storage(f.*, ldf, n + 1, 1, nrhs, kb.*);
                    f_rhs := array2d_view(f_rhs_raw, ldf, row_major=false);
                    b_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, n + 1, m - if_, nrhs);
                    b_tail := array2d_view(b_tail_raw, lda, row_major=false);
                    _laqp3rk_gemm(#char "N", #char "T", m - if_, nrhs, kb.*, -one_r, a_tail.*, lda, f_rhs.*, ldf, one_r, *b_tail, lda, RT);
                }
                return;
            }

            if maxc2_local == zero_r {
                done.* = true;
                kb.* = k - 1;
                if_ = i - 1;
                relmaxc2nrmk.* = zero_r;

                if nrhs > 0 && kb.* < (m - ioffset) {
                    a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
                    a_tail := array2d_view(a_tail_raw, lda, row_major=false);
                    f_rhs_raw := _laqp3rk_matrix_storage(f.*, ldf, n + 1, 1, nrhs, kb.*);
                    f_rhs := array2d_view(f_rhs_raw, ldf, row_major=false);
                    b_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, n + 1, m - if_, nrhs);
                    b_tail := array2d_view(b_tail_raw, lda, row_major=false);
                    _laqp3rk_gemm(#char "N", #char "T", m - if_, nrhs, kb.*, -one_r, a_tail.*, lda, f_rhs.*, ldf, one_r, *b_tail, lda, RT);
                }

                for j:k..minmnfact {
                    fset(tau, j, zero_r);
                }
                return;
            }

            if info.* == 0 && maxc2_local > hugeval {
                info.* = n + k - 1 + kp;
            }

            relmax_local = maxc2_local / maxc2nrm;
            relmaxc2nrmk.* = relmax_local;
            if maxc2_local <= abstol || relmax_local <= reltol {
                done.* = true;
                kb.* = k - 1;
                if_ = i - 1;

                if kb.* < minmnupdt {
                    a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
                    a_tail := array2d_view(a_tail_raw, lda, row_major=false);
                    f_tail_raw := _laqp3rk_matrix_storage(f.*, ldf, kb.* + 1, 1, n + nrhs - kb.*, kb.*);
                    f_tail := array2d_view(f_tail_raw, ldf, row_major=false);
                    residual_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, kb.* + 1, m - if_, n + nrhs - kb.*);
                    residual := array2d_view(residual_raw, lda, row_major=false);
                    _laqp3rk_gemm(#char "N", #char "T", m - if_, n + nrhs - kb.*, kb.*, -one_r, a_tail.*, lda, f_tail.*, ldf, one_r, *residual, lda, RT);
                }

                for j:k..minmnfact {
                    fset(tau, j, zero_r);
                }
                return;
            }
        }

        if kp != k {
            _laqp3rk_swap(m, _laqp3rk_matrix_column(a.*, lda, 1, kp, m), 1, _laqp3rk_matrix_column(a.*, lda, 1, k, m), 1, RT);
            _laqp3rk_swap(k - 1, _laqp3rk_matrix_row(f.*, ldf, kp, 1, k - 1), ldf, _laqp3rk_matrix_row(f.*, ldf, k, 1, k - 1), ldf, RT);
            fset(vn1, kp, fget(vn1.*, k));
            fset(vn2, kp, fget(vn2.*, k));
            itemp = fget(jpiv.*, kp);
            fset(jpiv, kp, fget(jpiv.*, k));
            fset(jpiv, k, itemp);
        }

        if k > 1 {
            a_col_raw := _laqp3rk_matrix_column(a.*, lda, i, 1, m - i + 1 + (k - 2)*lda);
            a_col_prev := array2d_view(a_col_raw, lda, row_major=false);
            f_row := _laqp3rk_matrix_row(f.*, ldf, k, 1, k - 1);
            a_k_col := _laqp3rk_matrix_column(a.*, lda, i, k, m - i + 1);
            _laqp3rk_gemv(#char "N", m - i + 1, k - 1, -one_r, a_col_prev.*, lda, f_row, ldf, one_r, *a_k_col, 1, RT);
        }

        if i < m {
            _laqp3rk_larfg(m - i + 1, a.data + (k - 1)*lda + (i - 1), *_laqp3rk_matrix_column(a.*, lda, i + 1, k, m - i), 1, tau.data + (k - 1), RT, RT);
        } else {
            fset(tau, k, zero_r);
        }

        taunan = _laqp3rk_tau_nan_real(fget(tau.*, k), RT);
        if _laqp3rk_isnan(taunan, RT) {
            done.* = true;
            kb.* = k - 1;
            if_ = i - 1;
            info.* = k;
            maxc2nrmk.* = taunan;
            relmaxc2nrmk.* = taunan;

            if nrhs > 0 && kb.* < (m - ioffset) {
                a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
                a_tail := array2d_view(a_tail_raw, lda, row_major=false);
                f_rhs_raw := _laqp3rk_matrix_storage(f.*, ldf, n + 1, 1, nrhs, kb.*);
                f_rhs := array2d_view(f_rhs_raw, ldf, row_major=false);
                b_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, n + 1, m - if_, nrhs);
                b_tail := array2d_view(b_tail_raw, lda, row_major=false);
                _laqp3rk_gemm(#char "N", #char "T", m - if_, nrhs, kb.*, -one_r, a_tail.*, lda, f_rhs.*, ldf, one_r, *b_tail, lda, RT);
            }
            return;
        }

        aik = fget(a.*, i, k);
        fset(a, i, k, one_r);

        if k < n + nrhs {
            a_trailing_raw := _laqp3rk_matrix_storage(a.*, lda, i, k + 1, m - i + 1, n + nrhs - k);
            a_trailing := array2d_view(a_trailing_raw, lda, row_major=false);
            a_k_col := _laqp3rk_matrix_column(a.*, lda, i, k, m - i + 1);
            f_k_col := _laqp3rk_matrix_column(f.*, ldf, k + 1, k, n + nrhs - k);
            _laqp3rk_gemv(#char "T", m - i + 1, n + nrhs - k, fget(tau.*, k), a_trailing.*, lda, a_k_col, 1, zero_r, *f_k_col, 1, RT);
        }

        for j:1..k {
            fset(f, j, k, zero_r);
        }

        if k > 1 {
            a_prev_raw := _laqp3rk_matrix_storage(a.*, lda, i, 1, m - i + 1, k - 1);
            a_prev := array2d_view(a_prev_raw, lda, row_major=false);
            a_k_col := _laqp3rk_matrix_column(a.*, lda, i, k, m - i + 1);
            auxv_col := _laqp3rk_raw_slice(auxv.data, k - 1);
            _laqp3rk_gemv(#char "T", m - i + 1, k - 1, -fget(tau.*, k), a_prev.*, lda, a_k_col, 1, zero_r, *auxv_col, 1, RT);

            f_prev_raw := _laqp3rk_matrix_storage(f.*, ldf, 1, 1, n + nrhs, k - 1);
            f_prev := array2d_view(f_prev_raw, ldf, row_major=false);
            f_k_col := _laqp3rk_matrix_column(f.*, ldf, 1, k, n + nrhs);
            _laqp3rk_gemv(#char "N", n + nrhs, k - 1, one_r, f_prev.*, ldf, auxv_col, 1, one_r, *f_k_col, 1, RT);
        }

        if k < n + nrhs {
            f_tail_raw := _laqp3rk_matrix_storage(f.*, ldf, k + 1, 1, n + nrhs - k, k);
            f_tail := array2d_view(f_tail_raw, ldf, row_major=false);
            a_i_row := _laqp3rk_matrix_row(a.*, lda, i, 1, k);
            a_i_tail := _laqp3rk_matrix_row(a.*, lda, i, k + 1, n + nrhs - k);
            _laqp3rk_gemv(#char "N", n + nrhs - k, k, -one_r, f_tail.*, ldf, a_i_row, lda, one_r, *a_i_tail, lda, RT);
        }

        fset(a, i, k, aik);

        if k < minmnfact {
            for j:k + 1..n {
                if fget(vn1.*, j) != zero_r {
                    temp = abs(fget(a.*, i, j)) / fget(vn1.*, j);
                    temp = max(zero_r, (one_r + temp)*(one_r - temp));
                    temp2 = temp*(fget(vn1.*, j) / fget(vn2.*, j))*(fget(vn1.*, j) / fget(vn2.*, j));
                    if temp2 <= tol3z {
                        fset(iwork, j - 1, lsticc);
                        lsticc = j;
                    } else {
                        fset(vn1, j, fget(vn1.*, j)*sqrt(temp));
                    }
                }
            }
        }
    }

    kb.* = k;
    if_ = i;

    if kb.* < minmnupdt {
        a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
        a_tail := array2d_view(a_tail_raw, lda, row_major=false);
        f_tail_raw := _laqp3rk_matrix_storage(f.*, ldf, kb.* + 1, 1, n + nrhs - kb.*, kb.*);
        f_tail := array2d_view(f_tail_raw, ldf, row_major=false);
        residual_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, kb.* + 1, m - if_, n + nrhs - kb.*);
        residual := array2d_view(residual_raw, lda, row_major=false);
        _laqp3rk_gemm(#char "N", #char "T", m - if_, n + nrhs - kb.*, kb.*, -one_r, a_tail.*, lda, f_tail.*, ldf, one_r, *residual, lda, RT);
    }

    while lsticc > 0 {
        itemp = fget(iwork.*, lsticc - 1);
        temp = _laqp3rk_nrm2(m - if_, _laqp3rk_matrix_column(a.*, lda, if_ + 1, lsticc, m - if_), 1, RT, RT);
        fset(vn1, lsticc, temp);
        fset(vn2, lsticc, temp);
        lsticc = itemp;
    }
}

_laqp3rk_complex :: (
    $$m: int, $$n: int, $$nrhs: int, $$ioffset: int, $$nb_: int, $$abstol: $RT,
    $$reltol: RT, $$kp1: int, $$maxc2nrm: RT, a: *$A/ArrayType2D, $$lda: int,
    done: *bool, kb: *int, maxc2nrmk: *RT, relmaxc2nrmk: *RT, jpiv: *$JPIV/ArrayType1D,
    tau: *$TAU/ArrayType1D, vn1: *$VN1/ArrayType1D, vn2: *$VN2/ArrayType1D,
    auxv: *$AUXV/ArrayType1D, f: *$F/ArrayType2D, $$ldf: int, iwork: *$IWORK/ArrayType1D,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           JPIV.data_type == int &&
           TAU.data_type == MT &&
           VN1.data_type == RT &&
           VN2.data_type == RT &&
           AUXV.data_type == MT &&
           F.data_type == MT &&
           IWORK.data_type == int;
} {
    zero_r :: #run zero(RT);
    zero_t :: #run zero(MT);
    one_r  :: #run one(RT);
    one_t  :: #run one(MT);

    aik : MT = ---;
    hugeval, maxc2_local, relmax_local, taunan, temp, temp2, tol3z : RT = ---;
    i, if_, itemp, j, k, kp, lsticc, minmnfact, minmnupdt, nb_use : int = ---;

    info.* = 0;
    minmnfact = min(m - ioffset, n);
    minmnupdt = min(m - ioffset, n + nrhs);
    nb_use = min(nb_, minmnfact);
    tol3z = sqrt(_laqp3rk_lamch_epsilon(RT));
    hugeval = _laqp3rk_lamch_overflow(RT);

    k = 0;
    i = ioffset;
    lsticc = 0;
    done.* = false;

    while k < nb_use && lsticc == 0 {
        k += 1;
        i = ioffset + k;

        if i == 1 {
            kp = kp1;
        } else {
            kp = (k - 1) + _laqp3rk_iamax_1based(n - k + 1, _laqp3rk_raw_slice(vn1.data + k - 1, n - k + 1), 1, RT);
            maxc2_local = fget(vn1.*, kp);
            maxc2nrmk.* = maxc2_local;

            if _laqp3rk_isnan(maxc2_local, RT) {
                done.* = true;
                kb.* = k - 1;
                if_ = i - 1;
                info.* = kb.* + kp;
                relmaxc2nrmk.* = maxc2_local;

                if nrhs > 0 && kb.* < (m - ioffset) {
                    a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
                    a_tail := array2d_view(a_tail_raw, lda, row_major=false);
                    f_rhs_raw := _laqp3rk_matrix_storage(f.*, ldf, n + 1, 1, nrhs, kb.*);
                    f_rhs := array2d_view(f_rhs_raw, ldf, row_major=false);
                    b_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, n + 1, m - if_, nrhs);
                    b_tail := array2d_view(b_tail_raw, lda, row_major=false);
                    _laqp3rk_gemm(#char "N", #char "C", m - if_, nrhs, kb.*, -one_t, a_tail.*, lda, f_rhs.*, ldf, one_t, *b_tail, lda, MT);
                }
                return;
            }

            if maxc2_local == zero_r {
                done.* = true;
                kb.* = k - 1;
                if_ = i - 1;
                relmaxc2nrmk.* = zero_r;

                if nrhs > 0 && kb.* < (m - ioffset) {
                    a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
                    a_tail := array2d_view(a_tail_raw, lda, row_major=false);
                    f_rhs_raw := _laqp3rk_matrix_storage(f.*, ldf, n + 1, 1, nrhs, kb.*);
                    f_rhs := array2d_view(f_rhs_raw, ldf, row_major=false);
                    b_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, n + 1, m - if_, nrhs);
                    b_tail := array2d_view(b_tail_raw, lda, row_major=false);
                    _laqp3rk_gemm(#char "N", #char "C", m - if_, nrhs, kb.*, -one_t, a_tail.*, lda, f_rhs.*, ldf, one_t, *b_tail, lda, MT);
                }

                for j:k..minmnfact {
                    fset(tau, j, zero_t);
                }
                return;
            }

            if info.* == 0 && maxc2_local > hugeval {
                info.* = n + k - 1 + kp;
            }

            relmax_local = maxc2_local / maxc2nrm;
            relmaxc2nrmk.* = relmax_local;
            if maxc2_local <= abstol || relmax_local <= reltol {
                done.* = true;
                kb.* = k - 1;
                if_ = i - 1;

                if kb.* < minmnupdt {
                    a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
                    a_tail := array2d_view(a_tail_raw, lda, row_major=false);
                    f_tail_raw := _laqp3rk_matrix_storage(f.*, ldf, kb.* + 1, 1, n + nrhs - kb.*, kb.*);
                    f_tail := array2d_view(f_tail_raw, ldf, row_major=false);
                    residual_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, kb.* + 1, m - if_, n + nrhs - kb.*);
                    residual := array2d_view(residual_raw, lda, row_major=false);
                    _laqp3rk_gemm(#char "N", #char "C", m - if_, n + nrhs - kb.*, kb.*, -one_t, a_tail.*, lda, f_tail.*, ldf, one_t, *residual, lda, MT);
                }

                for j:k..minmnfact {
                    fset(tau, j, zero_t);
                }
                return;
            }
        }

        if kp != k {
            _laqp3rk_swap(m, _laqp3rk_matrix_column(a.*, lda, 1, kp, m), 1, _laqp3rk_matrix_column(a.*, lda, 1, k, m), 1, MT);
            _laqp3rk_swap(k - 1, _laqp3rk_matrix_row(f.*, ldf, kp, 1, k - 1), ldf, _laqp3rk_matrix_row(f.*, ldf, k, 1, k - 1), ldf, MT);
            fset(vn1, kp, fget(vn1.*, k));
            fset(vn2, kp, fget(vn2.*, k));
            itemp = fget(jpiv.*, kp);
            fset(jpiv, kp, fget(jpiv.*, k));
            fset(jpiv, k, itemp);
        }

        if k > 1 {
            for j:1..k - 1 {
                fset(f, k, j, conjugate(fget(f.*, k, j)));
            }
            a_prev_raw := _laqp3rk_matrix_storage(a.*, lda, i, 1, m - i + 1, k - 1);
            a_prev := array2d_view(a_prev_raw, lda, row_major=false);
            f_k_row := _laqp3rk_matrix_row(f.*, ldf, k, 1, k - 1);
            a_k_col := _laqp3rk_matrix_column(a.*, lda, i, k, m - i + 1);
            _laqp3rk_gemv(#char "N", m - i + 1, k - 1, -one_t, a_prev.*, lda, f_k_row, ldf, one_t, *a_k_col, 1, MT);
            for j:1..k - 1 {
                fset(f, k, j, conjugate(fget(f.*, k, j)));
            }
        }

        if i < m {
            _laqp3rk_larfg(m - i + 1, a.data + (k - 1)*lda + (i - 1), *_laqp3rk_matrix_column(a.*, lda, i + 1, k, m - i), 1, tau.data + (k - 1), RT, MT);
        } else {
            fset(tau, k, zero_t);
        }

        taunan = _laqp3rk_tau_nan_complex(fget(tau.*, k), RT, MT);
        if _laqp3rk_isnan(taunan, RT) {
            done.* = true;
            kb.* = k - 1;
            if_ = i - 1;
            info.* = k;
            maxc2nrmk.* = taunan;
            relmaxc2nrmk.* = taunan;

            if nrhs > 0 && kb.* < (m - ioffset) {
                a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
                a_tail := array2d_view(a_tail_raw, lda, row_major=false);
                f_rhs_raw := _laqp3rk_matrix_storage(f.*, ldf, n + 1, 1, nrhs, kb.*);
                f_rhs := array2d_view(f_rhs_raw, ldf, row_major=false);
                b_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, n + 1, m - if_, nrhs);
                b_tail := array2d_view(b_tail_raw, lda, row_major=false);
                _laqp3rk_gemm(#char "N", #char "C", m - if_, nrhs, kb.*, -one_t, a_tail.*, lda, f_rhs.*, ldf, one_t, *b_tail, lda, MT);
            }
            return;
        }

        aik = fget(a.*, i, k);
        fset(a, i, k, one_t);

        if k < n + nrhs {
            a_trailing_raw := _laqp3rk_matrix_storage(a.*, lda, i, k + 1, m - i + 1, n + nrhs - k);
            a_trailing := array2d_view(a_trailing_raw, lda, row_major=false);
            a_k_col := _laqp3rk_matrix_column(a.*, lda, i, k, m - i + 1);
            f_k_col := _laqp3rk_matrix_column(f.*, ldf, k + 1, k, n + nrhs - k);
            _laqp3rk_gemv(#char "C", m - i + 1, n + nrhs - k, fget(tau.*, k), a_trailing.*, lda, a_k_col, 1, zero_t, *f_k_col, 1, MT);
        }

        for j:1..k {
            fset(f, j, k, zero_t);
        }

        if k > 1 {
            a_prev_raw := _laqp3rk_matrix_storage(a.*, lda, i, 1, m - i + 1, k - 1);
            a_prev := array2d_view(a_prev_raw, lda, row_major=false);
            a_k_col := _laqp3rk_matrix_column(a.*, lda, i, k, m - i + 1);
            auxv_col := _laqp3rk_raw_slice(auxv.data, k - 1);
            _laqp3rk_gemv(#char "C", m - i + 1, k - 1, -fget(tau.*, k), a_prev.*, lda, a_k_col, 1, zero_t, *auxv_col, 1, MT);

            f_prev_raw := _laqp3rk_matrix_storage(f.*, ldf, 1, 1, n + nrhs, k - 1);
            f_prev := array2d_view(f_prev_raw, ldf, row_major=false);
            f_k_col := _laqp3rk_matrix_column(f.*, ldf, 1, k, n + nrhs);
            _laqp3rk_gemv(#char "N", n + nrhs, k - 1, one_t, f_prev.*, ldf, auxv_col, 1, one_t, *f_k_col, 1, MT);
        }

        if k < n + nrhs {
            a_i_block_raw := _laqp3rk_matrix_storage(a.*, lda, i, 1, 1, k);
            a_i_block := array2d_view(a_i_block_raw, lda, row_major=false);
            f_tail_raw := _laqp3rk_matrix_storage(f.*, ldf, k + 1, 1, n + nrhs - k, k);
            f_tail := array2d_view(f_tail_raw, ldf, row_major=false);
            a_i_tail_raw := _laqp3rk_matrix_storage(a.*, lda, i, k + 1, 1, n + nrhs - k);
            a_i_tail := array2d_view(a_i_tail_raw, lda, row_major=false);
            _laqp3rk_gemm(#char "N", #char "C", 1, n + nrhs - k, k, -one_t, a_i_block.*, lda, f_tail.*, ldf, one_t, *a_i_tail, lda, MT);
        }

        fset(a, i, k, aik);

        if k < minmnfact {
            for j:k + 1..n {
                if fget(vn1.*, j) != zero_r {
                    temp = abs(fget(a.*, i, j)) / fget(vn1.*, j);
                    temp = max(zero_r, (one_r + temp)*(one_r - temp));
                    temp2 = temp*(fget(vn1.*, j) / fget(vn2.*, j))*(fget(vn1.*, j) / fget(vn2.*, j));
                    if temp2 <= tol3z {
                        fset(iwork, j - 1, lsticc);
                        lsticc = j;
                    } else {
                        fset(vn1, j, fget(vn1.*, j)*sqrt(temp));
                    }
                }
            }
        }
    }

    kb.* = k;
    if_ = i;

    if kb.* < minmnupdt {
        a_tail_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, 1, m - if_, kb.*);
        a_tail := array2d_view(a_tail_raw, lda, row_major=false);
        f_tail_raw := _laqp3rk_matrix_storage(f.*, ldf, kb.* + 1, 1, n + nrhs - kb.*, kb.*);
        f_tail := array2d_view(f_tail_raw, ldf, row_major=false);
        residual_raw := _laqp3rk_matrix_storage(a.*, lda, if_ + 1, kb.* + 1, m - if_, n + nrhs - kb.*);
        residual := array2d_view(residual_raw, lda, row_major=false);
        _laqp3rk_gemm(#char "N", #char "C", m - if_, n + nrhs - kb.*, kb.*, -one_t, a_tail.*, lda, f_tail.*, ldf, one_t, *residual, lda, MT);
    }

    while lsticc > 0 {
        itemp = fget(iwork.*, lsticc - 1);
        temp = _laqp3rk_nrm2(m - if_, _laqp3rk_matrix_column(a.*, lda, if_ + 1, lsticc, m - if_), 1, RT, MT);
        fset(vn1, lsticc, temp);
        fset(vn2, lsticc, temp);
        lsticc = itemp;
    }
}
