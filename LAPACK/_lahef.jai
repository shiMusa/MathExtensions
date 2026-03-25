// ? LAPACK *lahef()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment clahef/zlahef
// *LAHEF computes a partial factorization of a complex Hermitian
// indefinite matrix using the blocked Bunch-Kaufman method.

clahef :: #bake_arguments _lahef(MT=C32);
zlahef :: #bake_arguments _lahef(MT=C64);

#scope_file
_lahef_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_lahef_matrix_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    count := ifx rows <= 0 || cols <= 0 then 0 else lda*(cols - 1) + rows;
    return _lahef_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lahef_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _lahef_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lahef_matrix_row :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _lahef_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_lahef_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_lahef_real_type :: inline ($MT: Type) -> Type {
    return inner_type_of_complex(MT);
}

_lahef_one :: inline ($MT: Type) -> MT {
    RT :: #run inner_type_of_complex(MT);
    return _lahef_complex_from_real(one(RT), RT, MT);
}

_lahef_cabs1 :: inline (z: $MT, $MT: Type) -> inner_type_of_complex(MT) #modify {
    return MT == C32 || MT == C64;
} {
    RT :: #run inner_type_of_complex(MT);
    return abs(real(z)) + abs(imag(z));
}

_lahef_copy :: (
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

_lahef_gemv :: (
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

_lahef_lacgv :: (
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

_lahef_scal_real :: (
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

_lahef_swap :: (
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

_lahef_iamax :: ($$n: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type) -> int #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        return icamax(n, x, incx);
    } else {
        return izamax(n, x, incx);
    }
}

_lahef_gemmtr :: (
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

_lahef :: (
    $$uplo: u8,
    $$n: int,
    $$nb: int,
    kb: *int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    w: *$W/ArrayType2D,
    $$ldw: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           IPIV.data_type == int &&
           W.data_type == MT;
} {
    RT :: #run inner_type_of_complex(MT);
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    cone   :: #run _lahef_one(MT);
    eight  :: #run ncast(RT, 8.0);
    sevten :: #run ncast(RT, 17.0);
    alpha  : RT = (one_r + sqrt(sevten)) / eight;
    absakk, colmax, d11, d22, r1, rowmax, t : RT = ---;
    d21, z : MT = ---;
    i, imax, j, jj, jmax, jp, k, kk, kkw, kp, kstep, kw : int = ---;

    info.* = 0;
    if lsame(uplo, #char "U") {
        k = n;
        while true {
            kw = nb + k - n;
            if (k <= n - nb + 1 && nb < n) || k < 1 {
                break;
            }

            kstep = 1;
            a_col := _lahef_matrix_column(a.*, lda, 1, k, k - 1);
            w_col := _lahef_matrix_column(w.*, ldw, 1, kw, k - 1);
            _lahef_copy(k - 1, *a_col, 1, *w_col, 1, MT);
            fset(w, k, kw, _lahef_complex_from_real(real(fget(a.*, k, k)), RT, MT));
            if k < n {
                a_trail_raw := _lahef_matrix_view(a.*, lda, 1, k + 1, k, n - k);
                a_trail := array2d_view(a_trail_raw, lda, row_major=false);
                w_trail := _lahef_matrix_column(w.*, ldw, k, kw + 1, n - k);
                _lahef_gemv(#char "N", k, n - k, -cone, *a_trail, lda, *w_trail, ldw, cone, *w_col, 1, MT);
                fset(w, k, kw, _lahef_complex_from_real(real(fget(w.*, k, kw)), RT, MT));
            }

            absakk = abs(real(fget(w.*, k, kw)));
            if k > 1 {
                w_search := _lahef_matrix_column(w.*, ldw, 1, kw, k - 1);
                imax = _lahef_iamax(k - 1, *w_search, 1, MT);
                colmax = _lahef_cabs1(fget(w.*, imax, kw), MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r {
                if info.* == 0 {
                    info.* = k;
                }
                kp = k;
                fset(a, k, k, _lahef_complex_from_real(real(fget(a.*, k, k)), RT, MT));
            } else {
                if absakk >= alpha*colmax {
                    kp = k;
                } else {
                    w_imax_col := _lahef_matrix_column(w.*, ldw, 1, kw - 1, imax - 1);
                    a_imax_col1 := _lahef_matrix_column(a.*, lda, 1, imax, imax - 1);
                    _lahef_copy(imax - 1, *a_imax_col1, 1, *w_imax_col, 1, MT);
                    fset(w, imax, kw - 1, _lahef_complex_from_real(real(fget(a.*, imax, imax)), RT, MT));
                    a_imax_row := _lahef_matrix_row(a.*, lda, imax, imax + 1, k - imax);
                    w_imax_tail := _lahef_matrix_column(w.*, ldw, imax + 1, kw - 1, k - imax);
                    _lahef_copy(k - imax, *a_imax_row, lda, *w_imax_tail, 1, MT);
                    _lahef_lacgv(k - imax, *w_imax_tail, 1, MT);
                    if k < n {
                        a_trail_raw := _lahef_matrix_view(a.*, lda, 1, k + 1, k, n - k);
                        a_trail := array2d_view(a_trail_raw, lda, row_major=false);
                        w_imax_trail := _lahef_matrix_column(w.*, ldw, imax, kw + 1, n - k);
                        _lahef_gemv(#char "N", k, n - k, -cone, *a_trail, lda, *w_imax_trail, ldw, cone, *w_imax_col, 1, MT);
                        fset(w, imax, kw - 1, _lahef_complex_from_real(real(fget(w.*, imax, kw - 1)), RT, MT));
                    }

                    jmax = imax + _lahef_iamax(k - imax, *w_imax_tail, 1, MT);
                    rowmax = _lahef_cabs1(fget(w.*, jmax, kw - 1), MT);
                    if imax > 1 {
                        w_imax_head := _lahef_matrix_column(w.*, ldw, 1, kw - 1, imax - 1);
                        jmax = _lahef_iamax(imax - 1, *w_imax_head, 1, MT);
                        rowmax = max(rowmax, _lahef_cabs1(fget(w.*, jmax, kw - 1), MT));
                    }

                    if absakk >= alpha*colmax*(colmax / rowmax) {
                        kp = k;
                    } else if abs(real(fget(w.*, imax, kw - 1))) >= alpha*rowmax {
                        kp = imax;
                        _lahef_copy(k, *w_imax_col, 1, *w_col, 1, MT);
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                }

                kk = k - kstep + 1;
                kkw = nb + kk - n;

                if kp != kk {
                    fset(a, kp, kp, _lahef_complex_from_real(real(fget(a.*, kk, kk)), RT, MT));
                    a_src := _lahef_matrix_column(a.*, lda, kp + 1, kk, kk - 1 - kp);
                    a_dst := _lahef_matrix_row(a.*, lda, kp, kp + 1, kk - 1 - kp);
                    _lahef_copy(kk - 1 - kp, *a_src, 1, *a_dst, lda, MT);
                    _lahef_lacgv(kk - 1 - kp, *a_dst, lda, MT);
                    if kp > 1 {
                        a_head_src := _lahef_matrix_column(a.*, lda, 1, kk, kp - 1);
                        a_head_dst := _lahef_matrix_column(a.*, lda, 1, kp, kp - 1);
                        _lahef_copy(kp - 1, *a_head_src, 1, *a_head_dst, 1, MT);
                    }
                    if k < n {
                        a_row1 := _lahef_matrix_row(a.*, lda, kk, k + 1, n - k);
                        a_row2 := _lahef_matrix_row(a.*, lda, kp, k + 1, n - k);
                        _lahef_swap(n - k, *a_row1, lda, *a_row2, lda, MT);
                    }
                    w_row1 := _lahef_matrix_row(w.*, ldw, kk, kkw, n - kk + 1);
                    w_row2 := _lahef_matrix_row(w.*, ldw, kp, kkw, n - kk + 1);
                    _lahef_swap(n - kk + 1, *w_row1, ldw, *w_row2, ldw, MT);
                }

                if kstep == 1 {
                    a_store := _lahef_matrix_column(a.*, lda, 1, k, k);
                    _lahef_copy(k, *w_col, 1, *a_store, 1, MT);
                    if k > 1 {
                        r1 = one_r / real(fget(a.*, k, k));
                        a_scale := _lahef_matrix_column(a.*, lda, 1, k, k - 1);
                        _lahef_scal_real(k - 1, r1, *a_scale, 1, RT, MT);
                        _lahef_lacgv(k - 1, *w_col, 1, MT);
                    }
                } else {
                    if k > 2 {
                        d21 = fget(w.*, k - 1, kw);
                        d11 = real(fget(w.*, k, kw) / conjugate(d21));
                        d22 = real(fget(w.*, k - 1, kw - 1) / d21);
                        t = one_r / (d11*d22 - one_r);
                        z = _lahef_complex_from_real(t, RT, MT) / d21;
                        for j:1..k - 2 {
                            fset(a, j, k - 1, z*( _lahef_complex_from_real(d11, RT, MT)*fget(w.*, j, kw - 1) - fget(w.*, j, kw) ));
                            fset(a, j, k, conjugate(z)*( _lahef_complex_from_real(d22, RT, MT)*fget(w.*, j, kw) - fget(w.*, j, kw - 1) ));
                        }
                    }
                    fset(a, k - 1, k - 1, fget(w.*, k - 1, kw - 1));
                    fset(a, k - 1, k, fget(w.*, k - 1, kw));
                    fset(a, k, k, fget(w.*, k, kw));
                    _lahef_lacgv(k - 1, *w_col, 1, MT);
                    w_prev := _lahef_matrix_column(w.*, ldw, 1, kw - 1, k - 2);
                    _lahef_lacgv(k - 2, *w_prev, 1, MT);
                }
            }

            if kstep == 1 {
                fset(ipiv, k, kp);
            } else {
                fset(ipiv, k, -kp);
                fset(ipiv, k - 1, -kp);
            }

            k -= kstep;
        }

        kw = nb + k - n;
        a11_raw := _lahef_matrix_view(a.*, lda, 1, k + 1, k, n - k);
        a11 := array2d_view(a11_raw, lda, row_major=false);
        w11_raw := _lahef_matrix_view(w.*, ldw, 1, kw + 1, k, n - k);
        w11 := array2d_view(w11_raw, ldw, row_major=false);
        a_upd_raw := _lahef_matrix_view(a.*, lda, 1, 1, k, k);
        a_upd := array2d_view(a_upd_raw, lda, row_major=false);
        _lahef_gemmtr(#char "U", #char "N", #char "T", k, n - k, -cone, *a11, lda, *w11, ldw, cone, *a_upd, lda, MT);

        j = k + 1;
        while j <= n {
            jj = j;
            jp = fget(ipiv, j);
            if jp < 0 {
                jp = -jp;
                j += 1;
            }
            j += 1;
            if jp != jj && j <= n {
                row1 := _lahef_matrix_row(a.*, lda, jp, j, n - j + 1);
                row2 := _lahef_matrix_row(a.*, lda, jj, j, n - j + 1);
                _lahef_swap(n - j + 1, *row1, lda, *row2, lda, MT);
            }
        }

        kb.* = n - k;
    } else {
        k = 1;
        while true {
            if (k >= nb && nb < n) || k > n {
                break;
            }

            kstep = 1;
            fset(w, k, k, _lahef_complex_from_real(real(fget(a.*, k, k)), RT, MT));
            if k < n {
                a_col := _lahef_matrix_column(a.*, lda, k + 1, k, n - k);
                w_col := _lahef_matrix_column(w.*, ldw, k + 1, k, n - k);
                _lahef_copy(n - k, *a_col, 1, *w_col, 1, MT);
            }
            a_prev_raw := _lahef_matrix_view(a.*, lda, k, 1, n - k + 1, k - 1);
            a_prev := array2d_view(a_prev_raw, lda, row_major=false);
            w_prev := _lahef_matrix_row(w.*, ldw, k, 1, k - 1);
            w_colk := _lahef_matrix_column(w.*, ldw, k, k, n - k + 1);
            _lahef_gemv(#char "N", n - k + 1, k - 1, -cone, *a_prev, lda, *w_prev, ldw, cone, *w_colk, 1, MT);
            fset(w, k, k, _lahef_complex_from_real(real(fget(w.*, k, k)), RT, MT));

            absakk = abs(real(fget(w.*, k, k)));
            if k < n {
                w_tail := _lahef_matrix_column(w.*, ldw, k + 1, k, n - k);
                imax = k + _lahef_iamax(n - k, *w_tail, 1, MT);
                colmax = _lahef_cabs1(fget(w.*, imax, k), MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r {
                if info.* == 0 {
                    info.* = k;
                }
                kp = k;
                fset(a, k, k, _lahef_complex_from_real(real(fget(a.*, k, k)), RT, MT));
            } else {
                if absakk >= alpha*colmax {
                    kp = k;
                } else {
                    a_row1 := _lahef_matrix_row(a.*, lda, imax, k, imax - k);
                    w_col2 := _lahef_matrix_column(w.*, ldw, k, k + 1, imax - k);
                    _lahef_copy(imax - k, *a_row1, lda, *w_col2, 1, MT);
                    _lahef_lacgv(imax - k, *w_col2, 1, MT);
                    fset(w, imax, k + 1, _lahef_complex_from_real(real(fget(a.*, imax, imax)), RT, MT));
                    if imax < n {
                        a_col := _lahef_matrix_column(a.*, lda, imax + 1, imax, n - imax);
                        w_tail := _lahef_matrix_column(w.*, ldw, imax + 1, k + 1, n - imax);
                        _lahef_copy(n - imax, *a_col, 1, *w_tail, 1, MT);
                    }
                    a_prev2_raw := _lahef_matrix_view(a.*, lda, k, 1, n - k + 1, k - 1);
                    a_prev2 := array2d_view(a_prev2_raw, lda, row_major=false);
                    w_prev2 := _lahef_matrix_row(w.*, ldw, imax, 1, k - 1);
                    _lahef_gemv(#char "N", n - k + 1, k - 1, -cone, *a_prev2, lda, *w_prev2, ldw, cone, *w_col2, 1, MT);
                    fset(w, imax, k + 1, _lahef_complex_from_real(real(fget(w.*, imax, k + 1)), RT, MT));

                    jmax = k - 1 + _lahef_iamax(imax - k, *w_col2, 1, MT);
                    rowmax = _lahef_cabs1(fget(w.*, jmax, k + 1), MT);
                    if imax < n {
                        w_tail2 := _lahef_matrix_column(w.*, ldw, imax + 1, k + 1, n - imax);
                        jmax = imax + _lahef_iamax(n - imax, *w_tail2, 1, MT);
                        rowmax = max(rowmax, _lahef_cabs1(fget(w.*, jmax, k + 1), MT));
                    }

                    if absakk >= alpha*colmax*(colmax / rowmax) {
                        kp = k;
                    } else if abs(real(fget(w.*, imax, k + 1))) >= alpha*rowmax {
                        kp = imax;
                        _lahef_copy(n - k + 1, *w_col2, 1, *w_colk, 1, MT);
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                }

                kk = k + kstep - 1;
                if kp != kk {
                    fset(a, kp, kp, _lahef_complex_from_real(real(fget(a.*, kk, kk)), RT, MT));
                    a_src := _lahef_matrix_column(a.*, lda, kk + 1, kk, kp - kk - 1);
                    a_dst := _lahef_matrix_row(a.*, lda, kp, kk + 1, kp - kk - 1);
                    _lahef_copy(kp - kk - 1, *a_src, 1, *a_dst, lda, MT);
                    _lahef_lacgv(kp - kk - 1, *a_dst, lda, MT);
                    if kp < n {
                        a_tail_src := _lahef_matrix_column(a.*, lda, kp + 1, kk, n - kp);
                        a_tail_dst := _lahef_matrix_column(a.*, lda, kp + 1, kp, n - kp);
                        _lahef_copy(n - kp, *a_tail_src, 1, *a_tail_dst, 1, MT);
                    }
                    if k > 1 {
                        row1 := _lahef_matrix_row(a.*, lda, kk, 1, k - 1);
                        row2 := _lahef_matrix_row(a.*, lda, kp, 1, k - 1);
                        _lahef_swap(k - 1, *row1, lda, *row2, lda, MT);
                    }
                    w_row1 := _lahef_matrix_row(w.*, ldw, kk, 1, kk);
                    w_row2 := _lahef_matrix_row(w.*, ldw, kp, 1, kk);
                    _lahef_swap(kk, *w_row1, ldw, *w_row2, ldw, MT);
                }

                if kstep == 1 {
                    a_store := _lahef_matrix_column(a.*, lda, k, k, n - k + 1);
                    _lahef_copy(n - k + 1, *w_colk, 1, *a_store, 1, MT);
                    if k < n {
                        r1 = one_r / real(fget(a.*, k, k));
                        a_scale := _lahef_matrix_column(a.*, lda, k + 1, k, n - k);
                        _lahef_scal_real(n - k, r1, *a_scale, 1, RT, MT);
                        w_scale := _lahef_matrix_column(w.*, ldw, k + 1, k, n - k);
                        _lahef_lacgv(n - k, *w_scale, 1, MT);
                    }
                } else {
                    if k < n - 1 {
                        d21 = fget(w.*, k + 1, k);
                        d11 = real(fget(w.*, k + 1, k + 1) / d21);
                        d22 = real(fget(w.*, k, k) / conjugate(d21));
                        t = one_r / (d11*d22 - one_r);
                        z = _lahef_complex_from_real(t, RT, MT) / d21;
                        for j:k + 2..n {
                            fset(a, j, k, conjugate(z)*( _lahef_complex_from_real(d11, RT, MT)*fget(w.*, j, k) - fget(w.*, j, k + 1) ));
                            fset(a, j, k + 1, z*( _lahef_complex_from_real(d22, RT, MT)*fget(w.*, j, k + 1) - fget(w.*, j, k) ));
                        }
                    }
                    fset(a, k, k, fget(w.*, k, k));
                    fset(a, k + 1, k, fget(w.*, k + 1, k));
                    fset(a, k + 1, k + 1, fget(w.*, k + 1, k + 1));
                    w_k := _lahef_matrix_column(w.*, ldw, k + 1, k, n - k);
                    _lahef_lacgv(n - k, *w_k, 1, MT);
                    w_kp1 := _lahef_matrix_column(w.*, ldw, k + 2, k + 1, n - k - 1);
                    _lahef_lacgv(n - k - 1, *w_kp1, 1, MT);
                }
            }

            if kstep == 1 {
                fset(ipiv, k, kp);
            } else {
                fset(ipiv, k, -kp);
                fset(ipiv, k + 1, -kp);
            }

            k += kstep;
        }

        a22_raw := _lahef_matrix_view(a.*, lda, k, 1, n - k + 1, k - 1);
        a22 := array2d_view(a22_raw, lda, row_major=false);
        w22_raw := _lahef_matrix_view(w.*, ldw, k, 1, n - k + 1, k - 1);
        w22 := array2d_view(w22_raw, ldw, row_major=false);
        a_upd_raw := _lahef_matrix_view(a.*, lda, k, k, n - k + 1, n - k + 1);
        a_upd := array2d_view(a_upd_raw, lda, row_major=false);
        _lahef_gemmtr(#char "L", #char "N", #char "T", n - k + 1, k - 1, -cone, *a22, lda, *w22, ldw, cone, *a_upd, lda, MT);

        j = k - 1;
        while j >= 1 {
            jj = j;
            jp = fget(ipiv, j);
            if jp < 0 {
                jp = -jp;
                j -= 1;
            }
            j -= 1;
            if jp != jj && j >= 1 {
                row1 := _lahef_matrix_row(a.*, lda, jp, 1, j);
                row2 := _lahef_matrix_row(a.*, lda, jj, 1, j);
                _lahef_swap(j, *row1, lda, *row2, lda, MT);
            }
        }

        kb.* = k - 1;
    }
}
