// ? LAPACK *lahef_rk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clahef_rk/zlahef_rk
// *LAHEF_RK computes a partial factorization of a complex Hermitian
// indefinite matrix using the bounded Bunch-Kaufman rook method.

clahef_rk :: #bake_arguments _lahef_rk(MT=C32);
zlahef_rk :: #bake_arguments _lahef_rk(MT=C64);

#scope_file
_lahef_rk_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_lahef_rk_matrix_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else lda*(cols - 1) + rows;
    return _lahef_rk_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lahef_rk_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _lahef_rk_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lahef_rk_matrix_row :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _lahef_rk_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_lahef_rk_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_lahef_rk_one :: inline ($MT: Type) -> MT {
    RT :: #run inner_type_of_complex(MT);
    return _lahef_rk_complex_from_real(one(RT), RT, MT);
}

_lahef_rk_cabs1 :: inline (z: $MT, $MT: Type) -> inner_type_of_complex(MT) #modify {
    return MT == C32 || MT == C64;
} {
    RT :: #run inner_type_of_complex(MT);
    return abs(real(z)) + abs(imag(z));
}

_lahef_rk_copy :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_lahef_rk_gemv :: (
    $$trans: u8, $$m: int, $$n: int, alpha: $MT, a: *$A/ArrayType2D, $$lda: int,
    x: *$X/ArrayType1D, $$incx: int, beta: MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_lahef_rk_lacgv :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        clacgv(n, x, incx);
    } else {
        zlacgv(n, x, incx);
    }
}

_lahef_rk_scal_real :: (
    $$n: int, $$alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_lahef_rk_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        cswap(n, x, incx, y, incy);
    } else {
        zswap(n, x, incx, y, incy);
    }
}

_lahef_rk_iamax :: ($$n: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type) -> int #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        return icamax(n, x, incx);
    } else {
        return izamax(n, x, incx);
    }
}

_lahef_rk_gemmtr :: (
    $$uplo: u8, $$transa: u8, $$transb: u8, $$n: int, $$k: int,
    alpha: $MT, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    beta: MT, c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 {
        cgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_lahef_rk :: (
    $$uplo: u8,
    $$n: int,
    $$nb: int,
    kb: *int,
    a: *$A/ArrayType2D,
    $$lda: int,
    e: *$E/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    w: *$W/ArrayType2D,
    $$ldw: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           E.data_type == MT &&
           IPIV.data_type == int &&
           W.data_type == MT;
} {
    RT :: #run inner_type_of_complex(MT);
    zero_r :: #run zero(RT);
    zero_t :: #run zero(MT);
    one_r  :: #run one(RT);
    cone   :: #run _lahef_rk_one(MT);
    eight  :: #run ncast(RT, 8.0);
    sevten :: #run ncast(RT, 17.0);

    alpha : RT = (one_r + sqrt(sevten)) / eight;
    sfmin : RT = ifx RT == float32 then slamch("S") else dlamch("S");
    absakk, colmax, d11, d22, r1, rowmax, stemp, t : RT = ---;
    d21 : MT = ---;
    done : bool = ---;
    i, ii, imax, itemp, j, jmax, k, kk, kkw, kp, kstep, kw, p : int = ---;

    info.* = 0;

    if lsame(uplo, #char "U") {
        fset(e, 1, zero_t);
        k = n;
        while true {
            kw = nb + k - n;
            if (k <= n - nb + 1 && nb < n) || k < 1 {
                break;
            }

            kstep = 1;
            p = k;

            if k > 1 {
                a_col := _lahef_rk_matrix_column(a.*, lda, 1, k, k - 1);
                w_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw, k - 1);
                _lahef_rk_copy(k - 1, *a_col, 1, *w_col, 1, MT);
            }
            fset(w, k, kw, _lahef_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
            if k < n {
                a_trail_raw := _lahef_rk_matrix_view(a.*, lda, 1, k + 1, k, n - k);
                a_trail := array2d_view(a_trail_raw, lda, row_major=false);
                w_trail := _lahef_rk_matrix_column(w.*, ldw, k, kw + 1, n - k);
                w_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw, k);
                _lahef_rk_gemv(#char "N", k, n - k, -cone, *a_trail, lda, *w_trail, ldw, cone, *w_col, 1, MT);
                fset(w, k, kw, _lahef_rk_complex_from_real(real(fget(w.*, k, kw)), RT, MT));
            }

            absakk = abs(real(fget(w.*, k, kw)));
            if k > 1 {
                w_search := _lahef_rk_matrix_column(w.*, ldw, 1, kw, k - 1);
                imax = _lahef_rk_iamax(k - 1, *w_search, 1, MT);
                colmax = _lahef_rk_cabs1(fget(w.*, imax, kw), MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r {
                if info.* == 0 {
                    info.* = k;
                }
                kp = k;
                fset(a, k, k, _lahef_rk_complex_from_real(real(fget(w.*, k, kw)), RT, MT));
                if k > 1 {
                    w_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw, k - 1);
                    a_col := _lahef_rk_matrix_column(a.*, lda, 1, k, k - 1);
                    _lahef_rk_copy(k - 1, *w_col, 1, *a_col, 1, MT);
                    fset(e, k, zero_t);
                }
            } else {
                if !(absakk < alpha*colmax) {
                    kp = k;
                } else {
                    done = false;
                    while !done {
                        if imax > 1 {
                            a_imax_head := _lahef_rk_matrix_column(a.*, lda, 1, imax, imax - 1);
                            w_prev_head := _lahef_rk_matrix_column(w.*, ldw, 1, kw - 1, imax - 1);
                            _lahef_rk_copy(imax - 1, *a_imax_head, 1, *w_prev_head, 1, MT);
                        }
                        fset(w, imax, kw - 1, _lahef_rk_complex_from_real(real(fget(a.*, imax, imax)), RT, MT));
                        a_imax_row := _lahef_rk_matrix_row(a.*, lda, imax, imax + 1, k - imax);
                        w_prev_tail := _lahef_rk_matrix_column(w.*, ldw, imax + 1, kw - 1, k - imax);
                        _lahef_rk_copy(k - imax, *a_imax_row, lda, *w_prev_tail, 1, MT);
                        _lahef_rk_lacgv(k - imax, *w_prev_tail, 1, MT);

                        if k < n {
                            a_trail_raw := _lahef_rk_matrix_view(a.*, lda, 1, k + 1, k, n - k);
                            a_trail := array2d_view(a_trail_raw, lda, row_major=false);
                            w_imax_trail := _lahef_rk_matrix_column(w.*, ldw, imax, kw + 1, n - k);
                            w_prev_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw - 1, k);
                            _lahef_rk_gemv(#char "N", k, n - k, -cone, *a_trail, lda, *w_imax_trail, ldw, cone, *w_prev_col, 1, MT);
                            fset(w, imax, kw - 1, _lahef_rk_complex_from_real(real(fget(w.*, imax, kw - 1)), RT, MT));
                        }

                        if imax != k {
                            w_prev_imax_tail := _lahef_rk_matrix_column(w.*, ldw, imax + 1, kw - 1, k - imax);
                            jmax = imax + _lahef_rk_iamax(k - imax, *w_prev_imax_tail, 1, MT);
                            rowmax = _lahef_rk_cabs1(fget(w.*, jmax, kw - 1), MT);
                        } else {
                            rowmax = zero_r;
                        }
                        if imax > 1 {
                            w_prev_imax_head := _lahef_rk_matrix_column(w.*, ldw, 1, kw - 1, imax - 1);
                            itemp = _lahef_rk_iamax(imax - 1, *w_prev_imax_head, 1, MT);
                            stemp = _lahef_rk_cabs1(fget(w.*, itemp, kw - 1), MT);
                            if stemp > rowmax {
                                rowmax = stemp;
                                jmax = itemp;
                            }
                        }

                        if !(abs(real(fget(w.*, imax, kw - 1))) < alpha*rowmax) {
                            kp = imax;
                            w_prev_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw - 1, k);
                            w_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw, k);
                            _lahef_rk_copy(k, *w_prev_col, 1, *w_col, 1, MT);
                            done = true;
                        } else if (p == jmax) || (rowmax <= colmax) {
                            kp = imax;
                            kstep = 2;
                            done = true;
                        } else {
                            p = imax;
                            colmax = rowmax;
                            imax = jmax;
                            w_prev_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw - 1, k);
                            w_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw, k);
                            _lahef_rk_copy(k, *w_prev_col, 1, *w_col, 1, MT);
                        }
                    }
                }

                kk = k - kstep + 1;
                kkw = nb + kk - n;

                if (kstep == 2) && (p != k) {
                    fset(a, p, p, _lahef_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                    a_p_src := _lahef_rk_matrix_column(a.*, lda, p + 1, k, k - 1 - p);
                    a_p_dst := _lahef_rk_matrix_row(a.*, lda, p, p + 1, k - 1 - p);
                    _lahef_rk_copy(k - 1 - p, *a_p_src, 1, *a_p_dst, lda, MT);
                    _lahef_rk_lacgv(k - 1 - p, *a_p_dst, lda, MT);
                    if p > 1 {
                        a_head_src := _lahef_rk_matrix_column(a.*, lda, 1, k, p - 1);
                        a_head_dst := _lahef_rk_matrix_column(a.*, lda, 1, p, p - 1);
                        _lahef_rk_copy(p - 1, *a_head_src, 1, *a_head_dst, 1, MT);
                    }
                    if k < n {
                        a_rowk := _lahef_rk_matrix_row(a.*, lda, k, k + 1, n - k);
                        a_rowp := _lahef_rk_matrix_row(a.*, lda, p, k + 1, n - k);
                        _lahef_rk_swap(n - k, *a_rowk, lda, *a_rowp, lda, MT);
                    }
                    w_rowk := _lahef_rk_matrix_row(w.*, ldw, k, kkw, n - kk + 1);
                    w_rowp := _lahef_rk_matrix_row(w.*, ldw, p, kkw, n - kk + 1);
                    _lahef_rk_swap(n - kk + 1, *w_rowk, ldw, *w_rowp, ldw, MT);
                }

                if kp != kk {
                    fset(a, kp, kp, _lahef_rk_complex_from_real(real(fget(a.*, kk, kk)), RT, MT));
                    a_kp_src := _lahef_rk_matrix_column(a.*, lda, kp + 1, kk, kk - 1 - kp);
                    a_kp_dst := _lahef_rk_matrix_row(a.*, lda, kp, kp + 1, kk - 1 - kp);
                    _lahef_rk_copy(kk - 1 - kp, *a_kp_src, 1, *a_kp_dst, lda, MT);
                    _lahef_rk_lacgv(kk - 1 - kp, *a_kp_dst, lda, MT);
                    if kp > 1 {
                        a_head_src := _lahef_rk_matrix_column(a.*, lda, 1, kk, kp - 1);
                        a_head_dst := _lahef_rk_matrix_column(a.*, lda, 1, kp, kp - 1);
                        _lahef_rk_copy(kp - 1, *a_head_src, 1, *a_head_dst, 1, MT);
                    }
                    if k < n {
                        a_rowkk := _lahef_rk_matrix_row(a.*, lda, kk, k + 1, n - k);
                        a_rowkp := _lahef_rk_matrix_row(a.*, lda, kp, k + 1, n - k);
                        _lahef_rk_swap(n - k, *a_rowkk, lda, *a_rowkp, lda, MT);
                    }
                    w_rowkk := _lahef_rk_matrix_row(w.*, ldw, kk, kkw, n - kk + 1);
                    w_rowkp := _lahef_rk_matrix_row(w.*, ldw, kp, kkw, n - kk + 1);
                    _lahef_rk_swap(n - kk + 1, *w_rowkk, ldw, *w_rowkp, ldw, MT);
                }

                if kstep == 1 {
                    w_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw, k);
                    a_col := _lahef_rk_matrix_column(a.*, lda, 1, k, k);
                    _lahef_rk_copy(k, *w_col, 1, *a_col, 1, MT);
                    if k > 1 {
                        t = real(fget(a.*, k, k));
                        if abs(t) >= sfmin {
                            r1 = one_r / t;
                            a_scale := _lahef_rk_matrix_column(a.*, lda, 1, k, k - 1);
                            _lahef_rk_scal_real(k - 1, r1, *a_scale, 1, RT, MT);
                        } else {
                            for ii:1..k - 1 {
                                fset(a, ii, k, fget(a.*, ii, k) / _lahef_rk_complex_from_real(t, RT, MT));
                            }
                        }
                        _lahef_rk_lacgv(k - 1, *w_col, 1, MT);
                        fset(e, k, zero_t);
                    }
                } else {
                    if k > 2 {
                        d21 = fget(w.*, k - 1, kw);
                        d11 = real(fget(w.*, k, kw) / conjugate(d21));
                        d22 = real(fget(w.*, k - 1, kw - 1) / d21);
                        t = one_r / (real(d11*d22) - one_r);
                        for j:1..k - 2 {
                            fset(a, j, k - 1, _lahef_rk_complex_from_real(t, RT, MT)*(( _lahef_rk_complex_from_real(d11, RT, MT)*fget(w.*, j, kw - 1) - fget(w.*, j, kw) ) / d21));
                            fset(a, j, k, _lahef_rk_complex_from_real(t, RT, MT)*(( _lahef_rk_complex_from_real(d22, RT, MT)*fget(w.*, j, kw) - fget(w.*, j, kw - 1) ) / conjugate(d21)));
                        }
                    }
                    fset(a, k - 1, k - 1, fget(w.*, k - 1, kw - 1));
                    fset(a, k - 1, k, zero_t);
                    fset(a, k, k, fget(w.*, k, kw));
                    fset(e, k, fget(w.*, k - 1, kw));
                    fset(e, k - 1, zero_t);
                    w_col := _lahef_rk_matrix_column(w.*, ldw, 1, kw, k - 1);
                    _lahef_rk_lacgv(k - 1, *w_col, 1, MT);
                    w_prev := _lahef_rk_matrix_column(w.*, ldw, 1, kw - 1, k - 2);
                    _lahef_rk_lacgv(k - 2, *w_prev, 1, MT);
                }
            }

            if kstep == 1 {
                fset(ipiv, k, kp);
            } else {
                fset(ipiv, k, -p);
                fset(ipiv, k - 1, -kp);
            }

            k -= kstep;
        }

        a_upd_raw := _lahef_rk_matrix_view(a.*, lda, 1, 1, k, k);
        a_upd := array2d_view(a_upd_raw, lda, row_major=false);
        a_panel_raw := _lahef_rk_matrix_view(a.*, lda, 1, k + 1, k, n - k);
        a_panel := array2d_view(a_panel_raw, lda, row_major=false);
        w_panel_raw := _lahef_rk_matrix_view(w.*, ldw, 1, kw + 1, k, n - k);
        w_panel := array2d_view(w_panel_raw, ldw, row_major=false);
        _lahef_rk_gemmtr(#char "U", #char "N", #char "T", k, n - k, -cone, *a_panel, lda, *w_panel, ldw, cone, *a_upd, lda, MT);

        kb.* = n - k;
    } else {
        fset(e, n, zero_t);
        k = 1;
        while true {
            if (k >= nb && nb < n) || k > n {
                break;
            }

            kstep = 1;
            p = k;

            fset(w, k, k, _lahef_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
            if k < n {
                a_col := _lahef_rk_matrix_column(a.*, lda, k + 1, k, n - k);
                w_col := _lahef_rk_matrix_column(w.*, ldw, k + 1, k, n - k);
                _lahef_rk_copy(n - k, *a_col, 1, *w_col, 1, MT);
            }
            if k > 1 {
                a_prev_raw := _lahef_rk_matrix_view(a.*, lda, k, 1, n - k + 1, k - 1);
                a_prev := array2d_view(a_prev_raw, lda, row_major=false);
                w_prev := _lahef_rk_matrix_row(w.*, ldw, k, 1, k - 1);
                w_col := _lahef_rk_matrix_column(w.*, ldw, k, k, n - k + 1);
                _lahef_rk_gemv(#char "N", n - k + 1, k - 1, -cone, *a_prev, lda, *w_prev, ldw, cone, *w_col, 1, MT);
                fset(w, k, k, _lahef_rk_complex_from_real(real(fget(w.*, k, k)), RT, MT));
            }

            absakk = abs(real(fget(w.*, k, k)));
            if k < n {
                w_tail := _lahef_rk_matrix_column(w.*, ldw, k + 1, k, n - k);
                imax = k + _lahef_rk_iamax(n - k, *w_tail, 1, MT);
                colmax = _lahef_rk_cabs1(fget(w.*, imax, k), MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r {
                if info.* == 0 {
                    info.* = k;
                }
                kp = k;
                fset(a, k, k, _lahef_rk_complex_from_real(real(fget(w.*, k, k)), RT, MT));
                if k < n {
                    w_col := _lahef_rk_matrix_column(w.*, ldw, k + 1, k, n - k);
                    a_col := _lahef_rk_matrix_column(a.*, lda, k + 1, k, n - k);
                    _lahef_rk_copy(n - k, *w_col, 1, *a_col, 1, MT);
                    fset(e, k, zero_t);
                }
            } else {
                if !(absakk < alpha*colmax) {
                    kp = k;
                } else {
                    done = false;
                    while !done {
                        a_imax_row := _lahef_rk_matrix_row(a.*, lda, imax, k, imax - k);
                        w_next_head := _lahef_rk_matrix_column(w.*, ldw, k, k + 1, imax - k);
                        _lahef_rk_copy(imax - k, *a_imax_row, lda, *w_next_head, 1, MT);
                        _lahef_rk_lacgv(imax - k, *w_next_head, 1, MT);
                        fset(w, imax, k + 1, _lahef_rk_complex_from_real(real(fget(a.*, imax, imax)), RT, MT));
                        if imax < n {
                            a_imax_tail := _lahef_rk_matrix_column(a.*, lda, imax + 1, imax, n - imax);
                            w_next_tail := _lahef_rk_matrix_column(w.*, ldw, imax + 1, k + 1, n - imax);
                            _lahef_rk_copy(n - imax, *a_imax_tail, 1, *w_next_tail, 1, MT);
                        }
                        if k > 1 {
                            a_prev_raw := _lahef_rk_matrix_view(a.*, lda, k, 1, n - k + 1, k - 1);
                            a_prev := array2d_view(a_prev_raw, lda, row_major=false);
                            w_imax_prev := _lahef_rk_matrix_row(w.*, ldw, imax, 1, k - 1);
                            w_next_col := _lahef_rk_matrix_column(w.*, ldw, k, k + 1, n - k + 1);
                            _lahef_rk_gemv(#char "N", n - k + 1, k - 1, -cone, *a_prev, lda, *w_imax_prev, ldw, cone, *w_next_col, 1, MT);
                            fset(w, imax, k + 1, _lahef_rk_complex_from_real(real(fget(w.*, imax, k + 1)), RT, MT));
                        }

                        if imax != k {
                            w_next_head := _lahef_rk_matrix_column(w.*, ldw, k, k + 1, imax - k);
                            jmax = k - 1 + _lahef_rk_iamax(imax - k, *w_next_head, 1, MT);
                            rowmax = _lahef_rk_cabs1(fget(w.*, jmax, k + 1), MT);
                        } else {
                            rowmax = zero_r;
                        }
                        if imax < n {
                            w_next_tail := _lahef_rk_matrix_column(w.*, ldw, imax + 1, k + 1, n - imax);
                            itemp = imax + _lahef_rk_iamax(n - imax, *w_next_tail, 1, MT);
                            stemp = _lahef_rk_cabs1(fget(w.*, itemp, k + 1), MT);
                            if stemp > rowmax {
                                rowmax = stemp;
                                jmax = itemp;
                            }
                        }

                        if !(abs(real(fget(w.*, imax, k + 1))) < alpha*rowmax) {
                            kp = imax;
                            w_next_col := _lahef_rk_matrix_column(w.*, ldw, k, k + 1, n - k + 1);
                            w_col := _lahef_rk_matrix_column(w.*, ldw, k, k, n - k + 1);
                            _lahef_rk_copy(n - k + 1, *w_next_col, 1, *w_col, 1, MT);
                            done = true;
                        } else if (p == jmax) || (rowmax <= colmax) {
                            kp = imax;
                            kstep = 2;
                            done = true;
                        } else {
                            p = imax;
                            colmax = rowmax;
                            imax = jmax;
                            w_next_col := _lahef_rk_matrix_column(w.*, ldw, k, k + 1, n - k + 1);
                            w_col := _lahef_rk_matrix_column(w.*, ldw, k, k, n - k + 1);
                            _lahef_rk_copy(n - k + 1, *w_next_col, 1, *w_col, 1, MT);
                        }
                    }
                }

                kk = k + kstep - 1;

                if (kstep == 2) && (p != k) {
                    fset(a, p, p, _lahef_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                    a_p_src := _lahef_rk_matrix_column(a.*, lda, k + 1, k, p - k - 1);
                    a_p_dst := _lahef_rk_matrix_row(a.*, lda, p, k + 1, p - k - 1);
                    _lahef_rk_copy(p - k - 1, *a_p_src, 1, *a_p_dst, lda, MT);
                    _lahef_rk_lacgv(p - k - 1, *a_p_dst, lda, MT);
                    if p < n {
                        a_tail_src := _lahef_rk_matrix_column(a.*, lda, p + 1, k, n - p);
                        a_tail_dst := _lahef_rk_matrix_column(a.*, lda, p + 1, p, n - p);
                        _lahef_rk_copy(n - p, *a_tail_src, 1, *a_tail_dst, 1, MT);
                    }
                    if k > 1 {
                        a_rowk := _lahef_rk_matrix_row(a.*, lda, k, 1, k - 1);
                        a_rowp := _lahef_rk_matrix_row(a.*, lda, p, 1, k - 1);
                        _lahef_rk_swap(k - 1, *a_rowk, lda, *a_rowp, lda, MT);
                    }
                    w_rowk := _lahef_rk_matrix_row(w.*, ldw, k, 1, kk);
                    w_rowp := _lahef_rk_matrix_row(w.*, ldw, p, 1, kk);
                    _lahef_rk_swap(kk, *w_rowk, ldw, *w_rowp, ldw, MT);
                }

                if kp != kk {
                    fset(a, kp, kp, _lahef_rk_complex_from_real(real(fget(a.*, kk, kk)), RT, MT));
                    a_kp_src := _lahef_rk_matrix_column(a.*, lda, kk + 1, kk, kp - kk - 1);
                    a_kp_dst := _lahef_rk_matrix_row(a.*, lda, kp, kk + 1, kp - kk - 1);
                    _lahef_rk_copy(kp - kk - 1, *a_kp_src, 1, *a_kp_dst, lda, MT);
                    _lahef_rk_lacgv(kp - kk - 1, *a_kp_dst, lda, MT);
                    if kp < n {
                        a_tail_src := _lahef_rk_matrix_column(a.*, lda, kp + 1, kk, n - kp);
                        a_tail_dst := _lahef_rk_matrix_column(a.*, lda, kp + 1, kp, n - kp);
                        _lahef_rk_copy(n - kp, *a_tail_src, 1, *a_tail_dst, 1, MT);
                    }
                    if k > 1 {
                        a_rowkk := _lahef_rk_matrix_row(a.*, lda, kk, 1, k - 1);
                        a_rowkp := _lahef_rk_matrix_row(a.*, lda, kp, 1, k - 1);
                        _lahef_rk_swap(k - 1, *a_rowkk, lda, *a_rowkp, lda, MT);
                    }
                    w_rowkk := _lahef_rk_matrix_row(w.*, ldw, kk, 1, kk);
                    w_rowkp := _lahef_rk_matrix_row(w.*, ldw, kp, 1, kk);
                    _lahef_rk_swap(kk, *w_rowkk, ldw, *w_rowkp, ldw, MT);
                }

                if kstep == 1 {
                    w_col := _lahef_rk_matrix_column(w.*, ldw, k, k, n - k + 1);
                    a_col := _lahef_rk_matrix_column(a.*, lda, k, k, n - k + 1);
                    _lahef_rk_copy(n - k + 1, *w_col, 1, *a_col, 1, MT);
                    if k < n {
                        t = real(fget(a.*, k, k));
                        if abs(t) >= sfmin {
                            r1 = one_r / t;
                            a_scale := _lahef_rk_matrix_column(a.*, lda, k + 1, k, n - k);
                            _lahef_rk_scal_real(n - k, r1, *a_scale, 1, RT, MT);
                        } else {
                            for ii:k + 1..n {
                                fset(a, ii, k, fget(a.*, ii, k) / _lahef_rk_complex_from_real(t, RT, MT));
                            }
                        }
                        _lahef_rk_lacgv(n - k, *w_col, 1, MT);
                        fset(e, k, zero_t);
                    }
                } else {
                    if k < n - 1 {
                        d21 = fget(w.*, k + 1, k);
                        d11 = real(fget(w.*, k + 1, k + 1) / d21);
                        d22 = real(fget(w.*, k, k) / conjugate(d21));
                        t = one_r / (real(d11*d22) - one_r);
                        for j:k + 2..n {
                            fset(a, j, k, _lahef_rk_complex_from_real(t, RT, MT)*(( _lahef_rk_complex_from_real(d11, RT, MT)*fget(w.*, j, k) - fget(w.*, j, k + 1) ) / conjugate(d21)));
                            fset(a, j, k + 1, _lahef_rk_complex_from_real(t, RT, MT)*(( _lahef_rk_complex_from_real(d22, RT, MT)*fget(w.*, j, k + 1) - fget(w.*, j, k) ) / d21));
                        }
                    }
                    fset(a, k, k, fget(w.*, k, k));
                    fset(a, k + 1, k, zero_t);
                    fset(a, k + 1, k + 1, fget(w.*, k + 1, k + 1));
                    fset(e, k, fget(w.*, k + 1, k));
                    fset(e, k + 1, zero_t);
                    w_col := _lahef_rk_matrix_column(w.*, ldw, k + 1, k, n - k);
                    _lahef_rk_lacgv(n - k, *w_col, 1, MT);
                    w_next := _lahef_rk_matrix_column(w.*, ldw, k + 2, k + 1, n - k - 1);
                    _lahef_rk_lacgv(n - k - 1, *w_next, 1, MT);
                }
            }

            if kstep == 1 {
                fset(ipiv, k, kp);
            } else {
                fset(ipiv, k, -p);
                fset(ipiv, k + 1, -kp);
            }

            k += kstep;
        }

        a_upd_raw := _lahef_rk_matrix_view(a.*, lda, k, k, n - k + 1, n - k + 1);
        a_upd := array2d_view(a_upd_raw, lda, row_major=false);
        a_panel_raw := _lahef_rk_matrix_view(a.*, lda, k, 1, n - k + 1, k - 1);
        a_panel := array2d_view(a_panel_raw, lda, row_major=false);
        w_panel_raw := _lahef_rk_matrix_view(w.*, ldw, k, 1, n - k + 1, k - 1);
        w_panel := array2d_view(w_panel_raw, ldw, row_major=false);
        _lahef_rk_gemmtr(#char "L", #char "N", #char "T", n - k + 1, k - 1, -cone, *a_panel, lda, *w_panel, ldw, cone, *a_upd, lda, MT);

        kb.* = k - 1;
    }
}
