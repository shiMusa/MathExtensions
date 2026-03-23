// ? LAPACK *hetf2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetf2/zhetf2
// *HETF2 computes the factorization of a complex Hermitian matrix
// using the diagonal pivoting method (unblocked algorithm).

chetf2 :: #bake_arguments _hetf2(RT=float32, MT=C32);
zhetf2 :: #bake_arguments _hetf2(RT=float64, MT=C64);

#scope_file
_hetf2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetf2_matrix_column :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    return _hetf2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetf2_matrix_row :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _hetf2_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_hetf2_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetf2_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_hetf2_iamax_1based :: (
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

_hetf2_lapy2 :: (x: $RT, y: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return slapy2(x, y);
    } else {
        return dlapy2(x, y);
    }
}

_hetf2_isnan :: (x: $RT, $RT: Type) -> bool #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 {
        return sisnan(x);
    } else {
        return disnan(x);
    }
}

_hetf2_swap :: (
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

_hetf2_scal_real :: (
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

_hetf2_her :: (
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

_hetf2 :: (
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    ipiv: *$IPIV/ArrayType1D,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           IPIV.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    eight_r :: #run ncast(RT) 8.0;
    sevten_r :: #run ncast(RT) 17.0;
    routine_name :: #run if MT == C32 { "CHETF2" } else { "ZHETF2" };

    upper : bool = ---;
    i, imax, j, jmax, k, kk, kp, kstep : int = ---;
    absakk, alpha, colmax, d, d11, d22, r1, rowmax, tt : RT = ---;
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

    if upper {
        k = n;
        while true {
            if k < 1 then break;
            kstep = 1;

            absakk = abs(real(fget(a.*, k, k)));
            if k > 1 {
                colk := _hetf2_matrix_column(a.*, lda, 1, k, k - 1);
                imax = _hetf2_iamax_1based(k - 1, colk, 1, MT);
                colmax = _hetf2_abs1(fget(a.*, imax, k), RT, MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r || _hetf2_isnan(absakk, RT) {
                if info.* == 0 {
                    info.* = k;
                }
                kp = k;
                fset(a, k, k, _hetf2_complex_from_real(real(fget(a.*, k, k)), RT, MT));
            } else {
                if absakk >= alpha*colmax {
                    kp = k;
                } else {
                    row_imax := _hetf2_matrix_row(a.*, lda, imax, imax + 1, k - imax);
                    jmax = imax + _hetf2_iamax_1based(k - imax, row_imax, lda, MT);
                    rowmax = _hetf2_abs1(fget(a.*, imax, jmax), RT, MT);
                    if imax > 1 {
                        col_imax := _hetf2_matrix_column(a.*, lda, 1, imax, imax - 1);
                        jmax = _hetf2_iamax_1based(imax - 1, col_imax, 1, MT);
                        rowmax = max(rowmax, _hetf2_abs1(fget(a.*, jmax, imax), RT, MT));
                    }

                    if absakk >= alpha*colmax*(colmax / rowmax) {
                        kp = k;
                    } else if abs(real(fget(a.*, imax, imax))) >= alpha*rowmax {
                        kp = imax;
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                }

                kk = k - kstep + 1;
                if kp != kk {
                    if kp > 1 {
                        vkk := _hetf2_matrix_column(a.*, lda, 1, kk, kp - 1);
                        vkp := _hetf2_matrix_column(a.*, lda, 1, kp, kp - 1);
                        _hetf2_swap(kp - 1, *vkk, 1, *vkp, 1, MT);
                    }
                    for j:kp + 1..kk - 1 {
                        t = conjugate(fget(a.*, j, kk));
                        fset(a, j, kk, conjugate(fget(a.*, kp, j)));
                        fset(a, kp, j, t);
                    }
                    fset(a, kp, kk, conjugate(fget(a.*, kp, kk)));
                    r1 = real(fget(a.*, kk, kk));
                    fset(a, kk, kk, _hetf2_complex_from_real(real(fget(a.*, kp, kp)), RT, MT));
                    fset(a, kp, kp, _hetf2_complex_from_real(r1, RT, MT));
                    if kstep == 2 {
                        fset(a, k, k, _hetf2_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                        t = fget(a.*, k - 1, k);
                        fset(a, k - 1, k, fget(a.*, kp, k));
                        fset(a, kp, k, t);
                    }
                } else {
                    fset(a, k, k, _hetf2_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                    if kstep == 2 {
                        fset(a, k - 1, k - 1, _hetf2_complex_from_real(real(fget(a.*, k - 1, k - 1)), RT, MT));
                    }
                }

                if kstep == 1 {
                    if k > 1 {
                        r1 = one_r / real(fget(a.*, k, k));
                        vk := _hetf2_matrix_column(a.*, lda, 1, k, k - 1);
                        a11_raw := _hetf2_raw_slice(a.data, (k - 2)*lda + (k - 1));
                        a11 := array2d_view(a11_raw, lda, row_major=false);
                        _hetf2_her(uplo, k - 1, -r1, vk, 1, *a11, lda, RT, MT);
                        _hetf2_scal_real(k - 1, r1, *vk, 1, RT, MT);
                    }
                } else if k > 2 {
                    d = _hetf2_lapy2(real(fget(a.*, k - 1, k)), imag(fget(a.*, k - 1, k)), RT);
                    d22 = real(fget(a.*, k - 1, k - 1)) / d;
                    d11 = real(fget(a.*, k, k)) / d;
                    tt = one_r / (d11*d22 - one_r);
                    d12 = fget(a.*, k - 1, k) / d;
                    d = tt / d;

                    for <j:1..k - 2 {
                        wkm1 = d*(d11*fget(a.*, j, k - 1) - conjugate(d12)*fget(a.*, j, k));
                        wk = d*(d22*fget(a.*, j, k) - d12*fget(a.*, j, k - 1));
                        for <i:1..j {
                            fset(a, i, j, fget(a.*, i, j) - fget(a.*, i, k)*conjugate(wk) - fget(a.*, i, k - 1)*conjugate(wkm1));
                        }
                        fset(a, j, k, wk);
                        fset(a, j, k - 1, wkm1);
                        fset(a, j, j, _hetf2_complex_from_real(real(fget(a.*, j, j)), RT, MT));
                    }
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
    } else {
        k = 1;
        while true {
            if k > n then break;
            kstep = 1;

            absakk = abs(real(fget(a.*, k, k)));
            if k < n {
                colk := _hetf2_matrix_column(a.*, lda, k + 1, k, n - k);
                imax = k + _hetf2_iamax_1based(n - k, colk, 1, MT);
                colmax = _hetf2_abs1(fget(a.*, imax, k), RT, MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r || _hetf2_isnan(absakk, RT) {
                if info.* == 0 {
                    info.* = k;
                }
                kp = k;
                fset(a, k, k, _hetf2_complex_from_real(real(fget(a.*, k, k)), RT, MT));
            } else {
                if absakk >= alpha*colmax {
                    kp = k;
                } else {
                    row_imax := _hetf2_matrix_row(a.*, lda, imax, k, imax - k);
                    jmax = k - 1 + _hetf2_iamax_1based(imax - k, row_imax, lda, MT);
                    rowmax = _hetf2_abs1(fget(a.*, imax, jmax), RT, MT);
                    if imax < n {
                        col_imax := _hetf2_matrix_column(a.*, lda, imax + 1, imax, n - imax);
                        jmax = imax + _hetf2_iamax_1based(n - imax, col_imax, 1, MT);
                        rowmax = max(rowmax, _hetf2_abs1(fget(a.*, jmax, imax), RT, MT));
                    }

                    if absakk >= alpha*colmax*(colmax / rowmax) {
                        kp = k;
                    } else if abs(real(fget(a.*, imax, imax))) >= alpha*rowmax {
                        kp = imax;
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                }

                kk = k + kstep - 1;
                if kp != kk {
                    if kp < n {
                        vkk := _hetf2_matrix_column(a.*, lda, kp + 1, kk, n - kp);
                        vkp := _hetf2_matrix_column(a.*, lda, kp + 1, kp, n - kp);
                        _hetf2_swap(n - kp, *vkk, 1, *vkp, 1, MT);
                    }
                    for j:kk + 1..kp - 1 {
                        t = conjugate(fget(a.*, j, kk));
                        fset(a, j, kk, conjugate(fget(a.*, kp, j)));
                        fset(a, kp, j, t);
                    }
                    fset(a, kp, kk, conjugate(fget(a.*, kp, kk)));
                    r1 = real(fget(a.*, kk, kk));
                    fset(a, kk, kk, _hetf2_complex_from_real(real(fget(a.*, kp, kp)), RT, MT));
                    fset(a, kp, kp, _hetf2_complex_from_real(r1, RT, MT));
                    if kstep == 2 {
                        fset(a, k, k, _hetf2_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                        t = fget(a.*, k + 1, k);
                        fset(a, k + 1, k, fget(a.*, kp, k));
                        fset(a, kp, k, t);
                    }
                } else {
                    fset(a, k, k, _hetf2_complex_from_real(real(fget(a.*, k, k)), RT, MT));
                    if kstep == 2 {
                        fset(a, k + 1, k + 1, _hetf2_complex_from_real(real(fget(a.*, k + 1, k + 1)), RT, MT));
                    }
                }

                if kstep == 1 {
                    if k < n {
                        r1 = one_r / real(fget(a.*, k, k));
                        vk := _hetf2_matrix_column(a.*, lda, k + 1, k, n - k);
                        a22_raw := _hetf2_raw_slice(a.data + k*lda + k, (n - k - 1)*lda + (n - k));
                        a22 := array2d_view(a22_raw, lda, row_major=false);
                        _hetf2_her(uplo, n - k, -r1, vk, 1, *a22, lda, RT, MT);
                        _hetf2_scal_real(n - k, r1, *vk, 1, RT, MT);
                    }
                } else if k < n - 1 {
                    d = _hetf2_lapy2(real(fget(a.*, k + 1, k)), imag(fget(a.*, k + 1, k)), RT);
                    d11 = real(fget(a.*, k + 1, k + 1)) / d;
                    d22 = real(fget(a.*, k, k)) / d;
                    tt = one_r / (d11*d22 - one_r);
                    d21 = fget(a.*, k + 1, k) / d;
                    d = tt / d;

                    for j:k + 2..n {
                        wk = d*(d11*fget(a.*, j, k) - d21*fget(a.*, j, k + 1));
                        wkp1 = d*(d22*fget(a.*, j, k + 1) - conjugate(d21)*fget(a.*, j, k));
                        for i:j..n {
                            fset(a, i, j, fget(a.*, i, j) - fget(a.*, i, k)*conjugate(wk) - fget(a.*, i, k + 1)*conjugate(wkp1));
                        }
                        fset(a, j, k, wk);
                        fset(a, j, k + 1, wkp1);
                        fset(a, j, j, _hetf2_complex_from_real(real(fget(a.*, j, j)), RT, MT));
                    }
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
    }
}
