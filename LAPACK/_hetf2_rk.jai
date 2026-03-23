// ? LAPACK *hetf2_rk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetf2_rk/zhetf2_rk
// *HETF2_RK computes the factorization of a complex Hermitian matrix
// using the bounded Bunch-Kaufman rook pivoting method with a separate E.

chetf2_rk :: #bake_arguments _hetf2_rk(RT=float32, MT=C32);
zhetf2_rk :: #bake_arguments _hetf2_rk(RT=float64, MT=C64);

#scope_file
_hetf2_rk_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetf2_rk_matrix_column :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    return _hetf2_rk_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetf2_rk_matrix_row :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _hetf2_rk_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_hetf2_rk_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetf2_rk_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_hetf2_rk_iamax_1based :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    $MT: Type
) -> int #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        return icamax(n, x, incx) + 1;
    } else {
        return izamax(n, x, incx) + 1;
    }
}

_hetf2_rk_lamch_safe_minimum :: ($RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slamch(#char "s");
    } else {
        return dlamch(#char "s");
    }
}

_hetf2_rk_lapy2 :: (x: $RT, y: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slapy2(x, y);
    } else {
        return dlapy2(x, y);
    }
}

_hetf2_rk_swap :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        cswap(n, x, incx, y, incy);
    } else {
        zswap(n, x, incx, y, incy);
    }
}

_hetf2_rk_scal_real :: (
    $$n: int,
    alpha: $RT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        csscal(n, alpha, x, incx);
    } else {
        zdscal(n, alpha, x, incx);
    }
}

_hetf2_rk_her :: (
    $$uplo: u8,
    $$n: int,
    alpha: $RT,
    x: $X/ArrayType1D,
    $$incx: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT &&
           A.data_type == MT;
} {
    #if MT == C32 {
        cher(uplo, n, alpha, x, incx, a, lda);
    } else {
        zher(uplo, n, alpha, x, incx, a, lda);
    }
}

_hetf2_rk :: (
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    e: *$E/ArrayType1D,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           E.data_type == MT &&
           IPIV.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    eight_r :: #run ncast(RT) 8.0;
    sevten_r :: #run ncast(RT) 17.0;
    zero_t :: #run zero(MT);
    sfmin : RT = ---;
    routine_name :: #run if MT == C32 { "CHETF2_RK" } else { "ZHETF2_RK" };

    done, upper : bool = ---;
    i, ii, imax, itemp, j, jmax, k, kk, kp, kstep, p : int = ---;
    absakk, alpha, colmax, d, d11, d22, r1, stemp, rowmax, tt : RT = ---;
    d12, d21, t, wk, wkm1, wkp1 : MT = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    alpha = (one_r + sqrt(sevten_r)) / eight_r;
    sfmin = _hetf2_rk_lamch_safe_minimum(RT);

    if upper {
        fset(e, 1, zero_t);
        k = n;
        while true {
            if k < 1 then break;
            kstep = 1;
            p = k;

            absakk = abs(real(fget(a.*, k, k)));
            if k > 1 {
                colk := _hetf2_rk_matrix_column(a.*, lda, 1, k, k - 1);
                imax = _hetf2_rk_iamax_1based(k - 1, colk, 1, MT);
                colmax = _hetf2_rk_abs1(fget(a.*, imax, k), RT, MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r {
                if info.* == 0 {
                    info.* = k;
                }
                kp = k;
                fset(a, k, k, _hetf2_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                if k > 1 {
                    fset(e, k, zero_t);
                }
            } else {
                if !(absakk < alpha*colmax) {
                    kp = k;
                } else {
                    done = false;
                    while !done {
                        if imax != k {
                            row_imax := _hetf2_rk_matrix_row(a.*, lda, imax, imax + 1, k - imax);
                            jmax = imax + _hetf2_rk_iamax_1based(k - imax, row_imax, lda, MT);
                            rowmax = _hetf2_rk_abs1(fget(a.*, imax, jmax), RT, MT);
                        } else {
                            rowmax = zero_r;
                        }

                        if imax > 1 {
                            col_imax := _hetf2_rk_matrix_column(a.*, lda, 1, imax, imax - 1);
                            itemp = _hetf2_rk_iamax_1based(imax - 1, col_imax, 1, MT);
                            stemp = _hetf2_rk_abs1(fget(a.*, itemp, imax), RT, MT);
                            if stemp > rowmax {
                                rowmax = stemp;
                                jmax = itemp;
                            }
                        }

                        if !(abs(real(fget(a.*, imax, imax))) < alpha*rowmax) {
                            kp = imax;
                            done = true;
                        } else if p == jmax || rowmax <= colmax {
                            kp = imax;
                            kstep = 2;
                            done = true;
                        } else {
                            p = imax;
                            colmax = rowmax;
                            imax = jmax;
                        }
                    }
                }

                kk = k - kstep + 1;
                if kstep == 2 && p != k {
                    if p > 1 {
                        vk := _hetf2_rk_matrix_column(a.*, lda, 1, k, p - 1);
                        vp := _hetf2_rk_matrix_column(a.*, lda, 1, p, p - 1);
                        _hetf2_rk_swap(p - 1, *vk, 1, *vp, 1, MT);
                    }
                    for j:p + 1..k - 1 {
                        t = conjugate(fget(a.*, j, k));
                        fset(a, j, k, conjugate(fget(a.*, p, j)));
                        fset(a, p, j, t);
                    }
                    fset(a, p, k, conjugate(fget(a.*, p, k)));
                    r1 = real(fget(a.*, k, k));
                    fset(a, k, k, _hetf2_rk_complex_from_real(real(fget(a.*, p, p)), RT, MT));
                    fset(a, p, p, _hetf2_rk_complex_from_real(r1, RT, MT));
                    if k < n {
                        rowk := _hetf2_rk_matrix_row(a.*, lda, k, k + 1, n - k);
                        rowp := _hetf2_rk_matrix_row(a.*, lda, p, k + 1, n - k);
                        _hetf2_rk_swap(n - k, *rowk, lda, *rowp, lda, MT);
                    }
                }

                if kp != kk {
                    if kp > 1 {
                        vkk := _hetf2_rk_matrix_column(a.*, lda, 1, kk, kp - 1);
                        vkp := _hetf2_rk_matrix_column(a.*, lda, 1, kp, kp - 1);
                        _hetf2_rk_swap(kp - 1, *vkk, 1, *vkp, 1, MT);
                    }
                    for j:kp + 1..kk - 1 {
                        t = conjugate(fget(a.*, j, kk));
                        fset(a, j, kk, conjugate(fget(a.*, kp, j)));
                        fset(a, kp, j, t);
                    }
                    fset(a, kp, kk, conjugate(fget(a.*, kp, kk)));
                    r1 = real(fget(a.*, kk, kk));
                    fset(a, kk, kk, _hetf2_rk_complex_from_real(real(fget(a.*, kp, kp)), RT, MT));
                    fset(a, kp, kp, _hetf2_rk_complex_from_real(r1, RT, MT));
                    if kstep == 2 {
                        fset(a, k, k, _hetf2_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                        t = fget(a.*, k - 1, k);
                        fset(a, k - 1, k, fget(a.*, kp, k));
                        fset(a, kp, k, t);
                    }
                    if k < n {
                        rowkk := _hetf2_rk_matrix_row(a.*, lda, kk, k + 1, n - k);
                        rowkp := _hetf2_rk_matrix_row(a.*, lda, kp, k + 1, n - k);
                        _hetf2_rk_swap(n - k, *rowkk, lda, *rowkp, lda, MT);
                    }
                } else {
                    fset(a, k, k, _hetf2_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                    if kstep == 2 {
                        fset(a, k - 1, k - 1, _hetf2_rk_complex_from_real(real(fget(a.*, k - 1, k - 1)), RT, MT));
                    }
                }

                if kstep == 1 {
                    if k > 1 {
                        vk := _hetf2_rk_matrix_column(a.*, lda, 1, k, k - 1);
                        a11_raw := _hetf2_rk_raw_slice(a.data, (k - 2)*lda + (k - 1));
                        a11 := array2d_view(a11_raw, lda, row_major=false);
                        if abs(real(fget(a.*, k, k))) >= sfmin {
                            d11 = one_r / real(fget(a.*, k, k));
                            _hetf2_rk_her(uplo, k - 1, -d11, vk, 1, *a11, lda, RT, MT);
                            _hetf2_rk_scal_real(k - 1, d11, *vk, 1, RT, MT);
                        } else {
                            d11 = real(fget(a.*, k, k));
                            for ii:1..k - 1 {
                                fset(a, ii, k, fget(a.*, ii, k) / d11);
                            }
                            _hetf2_rk_her(uplo, k - 1, -d11, vk, 1, *a11, lda, RT, MT);
                        }
                        fset(e, k, zero_t);
                    }
                } else if k > 2 {
                    d = _hetf2_rk_lapy2(real(fget(a.*, k - 1, k)), imag(fget(a.*, k - 1, k)), RT);
                    d11 = real(fget(a.*, k, k) / d);
                    d22 = real(fget(a.*, k - 1, k - 1) / d);
                    d12 = fget(a.*, k - 1, k) / d;
                    tt = one_r / (d11*d22 - one_r);

                    for <j:1..k - 2 {
                        wkm1 = tt*(d11*fget(a.*, j, k - 1) - conjugate(d12)*fget(a.*, j, k));
                        wk = tt*(d22*fget(a.*, j, k) - d12*fget(a.*, j, k - 1));
                        for <i:1..j {
                            fset(a, i, j, fget(a.*, i, j) - (fget(a.*, i, k) / d)*conjugate(wk) - (fget(a.*, i, k - 1) / d)*conjugate(wkm1));
                        }
                        fset(a, j, k, wk / d);
                        fset(a, j, k - 1, wkm1 / d);
                        fset(a, j, j, _hetf2_rk_complex_from_real(real(fget(a.*, j, j)), RT, MT));
                    }
                    fset(e, k, fget(a.*, k - 1, k));
                    fset(e, k - 1, zero_t);
                    fset(a, k - 1, k, zero_t);
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
    } else {
        fset(e, n, zero_t);
        k = 1;
        while true {
            if k > n then break;
            kstep = 1;
            p = k;

            absakk = abs(real(fget(a.*, k, k)));
            if k < n {
                colk := _hetf2_rk_matrix_column(a.*, lda, k + 1, k, n - k);
                imax = k + _hetf2_rk_iamax_1based(n - k, colk, 1, MT);
                colmax = _hetf2_rk_abs1(fget(a.*, imax, k), RT, MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r {
                if info.* == 0 {
                    info.* = k;
                }
                kp = k;
                fset(a, k, k, _hetf2_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                if k < n {
                    fset(e, k, zero_t);
                }
            } else {
                if !(absakk < alpha*colmax) {
                    kp = k;
                } else {
                    done = false;
                    while !done {
                        if imax != k {
                            row_imax := _hetf2_rk_matrix_row(a.*, lda, imax, k, imax - k);
                            jmax = k - 1 + _hetf2_rk_iamax_1based(imax - k, row_imax, lda, MT);
                            rowmax = _hetf2_rk_abs1(fget(a.*, imax, jmax), RT, MT);
                        } else {
                            rowmax = zero_r;
                        }

                        if imax < n {
                            col_imax := _hetf2_rk_matrix_column(a.*, lda, imax + 1, imax, n - imax);
                            itemp = imax + _hetf2_rk_iamax_1based(n - imax, col_imax, 1, MT);
                            stemp = _hetf2_rk_abs1(fget(a.*, itemp, imax), RT, MT);
                            if stemp > rowmax {
                                rowmax = stemp;
                                jmax = itemp;
                            }
                        }

                        if !(abs(real(fget(a.*, imax, imax))) < alpha*rowmax) {
                            kp = imax;
                            done = true;
                        } else if p == jmax || rowmax <= colmax {
                            kp = imax;
                            kstep = 2;
                            done = true;
                        } else {
                            p = imax;
                            colmax = rowmax;
                            imax = jmax;
                        }
                    }
                }

                kk = k + kstep - 1;
                if kstep == 2 && p != k {
                    if p < n {
                        vk := _hetf2_rk_matrix_column(a.*, lda, p + 1, k, n - p);
                        vp := _hetf2_rk_matrix_column(a.*, lda, p + 1, p, n - p);
                        _hetf2_rk_swap(n - p, *vk, 1, *vp, 1, MT);
                    }
                    for j:k + 1..p - 1 {
                        t = conjugate(fget(a.*, j, k));
                        fset(a, j, k, conjugate(fget(a.*, p, j)));
                        fset(a, p, j, t);
                    }
                    fset(a, p, k, conjugate(fget(a.*, p, k)));
                    r1 = real(fget(a.*, k, k));
                    fset(a, k, k, _hetf2_rk_complex_from_real(real(fget(a.*, p, p)), RT, MT));
                    fset(a, p, p, _hetf2_rk_complex_from_real(r1, RT, MT));
                    if k > 1 {
                        rowk := _hetf2_rk_matrix_row(a.*, lda, k, 1, k - 1);
                        rowp := _hetf2_rk_matrix_row(a.*, lda, p, 1, k - 1);
                        _hetf2_rk_swap(k - 1, *rowk, lda, *rowp, lda, MT);
                    }
                }

                if kp != kk {
                    if kp < n {
                        vkk := _hetf2_rk_matrix_column(a.*, lda, kp + 1, kk, n - kp);
                        vkp := _hetf2_rk_matrix_column(a.*, lda, kp + 1, kp, n - kp);
                        _hetf2_rk_swap(n - kp, *vkk, 1, *vkp, 1, MT);
                    }
                    for j:kk + 1..kp - 1 {
                        t = conjugate(fget(a.*, j, kk));
                        fset(a, j, kk, conjugate(fget(a.*, kp, j)));
                        fset(a, kp, j, t);
                    }
                    fset(a, kp, kk, conjugate(fget(a.*, kp, kk)));
                    r1 = real(fget(a.*, kk, kk));
                    fset(a, kk, kk, _hetf2_rk_complex_from_real(real(fget(a.*, kp, kp)), RT, MT));
                    fset(a, kp, kp, _hetf2_rk_complex_from_real(r1, RT, MT));
                    if kstep == 2 {
                        fset(a, k, k, _hetf2_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                        t = fget(a.*, k + 1, k);
                        fset(a, k + 1, k, fget(a.*, kp, k));
                        fset(a, kp, k, t);
                    }
                    if k > 1 {
                        rowkk := _hetf2_rk_matrix_row(a.*, lda, kk, 1, k - 1);
                        rowkp := _hetf2_rk_matrix_row(a.*, lda, kp, 1, k - 1);
                        _hetf2_rk_swap(k - 1, *rowkk, lda, *rowkp, lda, MT);
                    }
                } else {
                    fset(a, k, k, _hetf2_rk_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                    if kstep == 2 {
                        fset(a, k + 1, k + 1, _hetf2_rk_complex_from_real(real(fget(a.*, k + 1, k + 1)), RT, MT));
                    }
                }

                if kstep == 1 {
                    if k < n {
                        vk := _hetf2_rk_matrix_column(a.*, lda, k + 1, k, n - k);
                        a22_raw := _hetf2_rk_raw_slice(a.data + k*lda + k, (n - k - 1)*lda + (n - k));
                        a22 := array2d_view(a22_raw, lda, row_major=false);
                        if abs(real(fget(a.*, k, k))) >= sfmin {
                            d11 = one_r / real(fget(a.*, k, k));
                            _hetf2_rk_her(uplo, n - k, -d11, vk, 1, *a22, lda, RT, MT);
                            _hetf2_rk_scal_real(n - k, d11, *vk, 1, RT, MT);
                        } else {
                            d11 = real(fget(a.*, k, k));
                            for ii:k + 1..n {
                                fset(a, ii, k, fget(a.*, ii, k) / d11);
                            }
                            _hetf2_rk_her(uplo, n - k, -d11, vk, 1, *a22, lda, RT, MT);
                        }
                        fset(e, k, zero_t);
                    }
                } else if k < n - 1 {
                    d = _hetf2_rk_lapy2(real(fget(a.*, k + 1, k)), imag(fget(a.*, k + 1, k)), RT);
                    d11 = real(fget(a.*, k + 1, k + 1)) / d;
                    d22 = real(fget(a.*, k, k)) / d;
                    d21 = fget(a.*, k + 1, k) / d;
                    tt = one_r / (d11*d22 - one_r);

                    for j:k + 2..n {
                        wk = tt*(d11*fget(a.*, j, k) - d21*fget(a.*, j, k + 1));
                        wkp1 = tt*(d22*fget(a.*, j, k + 1) - conjugate(d21)*fget(a.*, j, k));
                        for i:j..n {
                            fset(a, i, j, fget(a.*, i, j) - (fget(a.*, i, k) / d)*conjugate(wk) - (fget(a.*, i, k + 1) / d)*conjugate(wkp1));
                        }
                        fset(a, j, k, wk / d);
                        fset(a, j, k + 1, wkp1 / d);
                        fset(a, j, j, _hetf2_rk_complex_from_real(real(fget(a.*, j, j)), RT, MT));
                    }
                    fset(e, k, fget(a.*, k + 1, k));
                    fset(e, k + 1, zero_t);
                    fset(a, k + 1, k, zero_t);
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
    }
}
