// ? LAPACK *hptrf()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chptrf/zhptrf
// *HPTRF computes the factorization of a complex Hermitian packed
// matrix using the Bunch-Kaufman diagonal pivoting method.

chptrf :: #bake_arguments _hptrf(RT=float32, MT=C32);
zhptrf :: #bake_arguments _hptrf(RT=float64, MT=C64);

#scope_file
_hptrf_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hptrf_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hptrf_raw_slice(x.data + offset - 1, count);
}

_hptrf_upper_col_start :: inline ($$j: int) -> int {
    return 1 + (j - 1)*j/2;
}

_hptrf_lower_col_start :: inline ($$n: int, $$j: int) -> int {
    return 1 + (j - 1)*(2*n - j + 2)/2;
}

_hptrf_upper_index :: inline ($$i: int, $$j: int) -> int {
    return i + (j - 1)*j/2;
}

_hptrf_lower_index :: inline ($$n: int, $$i: int, $$j: int) -> int {
    return i + (j - 1)*(2*n - j)/2;
}

_hptrf_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hptrf_abs1 :: inline (z: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return _abs1(z);
}

_hptrf_iamax_1based :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, $MT: Type
) -> int #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { return icamax(n, x, incx) + 1; } else { return izamax(n, x, incx) + 1; }
}

_hptrf_lapy2 :: (x: $RT, y: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return slapy2(x, y); } else { return dlapy2(x, y); }
}

_hptrf_hpr :: (
    $$uplo: u8, $$n: int, alpha: $RT, x: $X/ArrayType1D, $$incx: int, ap: *$AP/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT && AP.data_type == MT;
} {
    #if MT == C32 { chpr(uplo, n, alpha, x, incx, ap); } else { zhpr(uplo, n, alpha, x, incx, ap); }
}

_hptrf_scal_real :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_hptrf_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hptrf :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, ipiv: *$IPIV/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && IPIV.data_type == int;
} {
    zero_r :: #run zero(RT);
    one_r :: #run one(RT);
    eight_r :: #run ncast(RT) 8.0;
    sevten_r :: #run ncast(RT) 17.0;
    routine_name :: #run if MT == C32 { "CHPTRF" } else { "ZHPTRF" };

    upper : bool = ---;
    i, imax, j, jmax, k, kc, kk, knc, kp, kpc, kstep, kx, npp : int = ---;
    absakk, alpha, colmax, d, d11, d22, r1, rowmax, tt : RT = ---;
    d12, d21, t, wk, wkm1, wkp1 : MT = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    alpha = (one_r + sqrt(sevten_r)) / eight_r;

    if upper {
        k = n;
        kc = (n - 1)*n/2 + 1;
        while true {
            knc = kc;
            if k < 1 then break;
            kstep = 1;

            absakk = abs(real(fget(ap.*, kc + k - 1)));
            if k > 1 {
                col_k := _hptrf_subvector(ap.*, kc, k - 1);
                imax = _hptrf_iamax_1based(k - 1, col_k, 1, MT);
                colmax = _hptrf_abs1(fget(ap.*, kc + imax - 1), RT, MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r {
                if info.* == 0 { info.* = k; }
                kp = k;
                fset(ap, kc + k - 1, _hptrf_complex_from_real(real(fget(ap.*, kc + k - 1)), RT, MT));
            } else {
                if absakk >= alpha*colmax {
                    kp = k;
                } else {
                    rowmax = zero_r;
                    jmax = imax;
                    kx = imax*(imax + 1)/2 + imax;
                    for j:imax + 1..k {
                        if _hptrf_abs1(fget(ap.*, kx), RT, MT) > rowmax {
                            rowmax = _hptrf_abs1(fget(ap.*, kx), RT, MT);
                            jmax = j;
                        }
                        kx += j;
                    }
                    kpc = _hptrf_upper_col_start(imax);
                    if imax > 1 {
                        jmax = _hptrf_iamax_1based(imax - 1, _hptrf_subvector(ap.*, kpc, imax - 1), 1, MT);
                        rowmax = max(rowmax, _hptrf_abs1(fget(ap.*, kpc + jmax - 1), RT, MT));
                    }

                    if absakk >= alpha*colmax*(colmax / rowmax) {
                        kp = k;
                    } else if abs(real(fget(ap.*, kpc + imax - 1))) >= alpha*rowmax {
                        kp = imax;
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                }

                kk = k - kstep + 1;
                if kstep == 2 {
                    knc = knc - k + 1;
                }
                if kp != kk {
                    kpc = _hptrf_upper_col_start(kp);
                    if kp > 1 {
                        col_kk := _hptrf_subvector(ap.*, knc, kp - 1);
                        col_kp := _hptrf_subvector(ap.*, kpc, kp - 1);
                        _hptrf_swap(kp - 1, *col_kk, 1, *col_kp, 1, MT);
                    }
                    kx = kpc + kp - 1;
                    for j:kp + 1..kk - 1 {
                        kx += j - 1;
                        t = conjugate(fget(ap.*, knc + j - 1));
                        fset(ap, knc + j - 1, conjugate(fget(ap.*, kx)));
                        fset(ap, kx, t);
                    }
                    fset(ap, kx + kk - 1, conjugate(fget(ap.*, kx + kk - 1)));
                    r1 = real(fget(ap.*, knc + kk - 1));
                    fset(ap, knc + kk - 1, _hptrf_complex_from_real(real(fget(ap.*, kpc + kp - 1)), RT, MT));
                    fset(ap, kpc + kp - 1, _hptrf_complex_from_real(r1, RT, MT));
                    if kstep == 2 {
                        fset(ap, kc + k - 1, _hptrf_complex_from_real(real(fget(ap.*, kc + k - 1)), RT, MT));
                        t = fget(ap.*, kc + k - 2);
                        fset(ap, kc + k - 2, fget(ap.*, kc + kp - 1));
                        fset(ap, kc + kp - 1, t);
                    }
                } else {
                    fset(ap, kc + k - 1, _hptrf_complex_from_real(real(fget(ap.*, kc + k - 1)), RT, MT));
                    if kstep == 2 {
                        fset(ap, kc - 1, _hptrf_complex_from_real(real(fget(ap.*, kc - 1)), RT, MT));
                    }
                }

                if kstep == 1 {
                    r1 = one_r / real(fget(ap.*, kc + k - 1));
                    if k > 1 {
                        col_k := _hptrf_subvector(ap.*, kc, k - 1);
                        ap_head := _hptrf_subvector(ap.*, 1, (k - 1)*k/2);
                        _hptrf_hpr(uplo, k - 1, -r1, col_k, 1, *ap_head, RT, MT);
                        _hptrf_scal_real(k - 1, r1, *col_k, 1, RT, MT);
                    }
                } else if k > 2 {
                    d = _hptrf_lapy2(real(fget(ap.*, kc + k - 2)), imag(fget(ap.*, kc + k - 2)), RT);
                    d22 = real(fget(ap.*, kc - 1)) / d;
                    d11 = real(fget(ap.*, kc + k - 1)) / d;
                    tt = one_r / (d11*d22 - one_r);
                    d12 = fget(ap.*, kc + k - 2) / d;
                    d = tt / d;

                    for <j:1..k - 2 {
                        idx_j_km1 :: _hptrf_upper_index(j, k - 1);
                        idx_j_k :: _hptrf_upper_index(j, k);
                        wkm1 = d*(d11*fget(ap.*, idx_j_km1) - conjugate(d12)*fget(ap.*, idx_j_k));
                        wk = d*(d22*fget(ap.*, idx_j_k) - d12*fget(ap.*, idx_j_km1));
                        for <i:1..j {
                            idx_ij :: _hptrf_upper_index(i, j);
                            fset(ap, idx_ij, fget(ap.*, idx_ij) - fget(ap.*, _hptrf_upper_index(i, k))*conjugate(wk) - fget(ap.*, _hptrf_upper_index(i, k - 1))*conjugate(wkm1));
                        }
                        fset(ap, idx_j_k, wk);
                        fset(ap, idx_j_km1, wkm1);
                        diag_j :: _hptrf_upper_index(j, j);
                        fset(ap, diag_j, _hptrf_complex_from_real(real(fget(ap.*, diag_j)), RT, MT));
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
            kc = knc - k;
        }
    } else {
        k = 1;
        kc = 1;
        npp = n*(n + 1)/2;
        while true {
            knc = kc;
            if k > n then break;
            kstep = 1;

            absakk = abs(real(fget(ap.*, kc)));
            if k < n {
                col_k := _hptrf_subvector(ap.*, kc + 1, n - k);
                imax = k + _hptrf_iamax_1based(n - k, col_k, 1, MT);
                colmax = _hptrf_abs1(fget(ap.*, kc + imax - k), RT, MT);
            } else {
                colmax = zero_r;
            }

            if max(absakk, colmax) == zero_r {
                if info.* == 0 { info.* = k; }
                kp = k;
                fset(ap, kc, _hptrf_complex_from_real(real(fget(ap.*, kc)), RT, MT));
            } else {
                if absakk >= alpha*colmax {
                    kp = k;
                } else {
                    rowmax = zero_r;
                    kx = kc + imax - k;
                    for j:k..imax - 1 {
                        if _hptrf_abs1(fget(ap.*, kx), RT, MT) > rowmax {
                            rowmax = _hptrf_abs1(fget(ap.*, kx), RT, MT);
                            jmax = j;
                        }
                        kx += n - j;
                    }
                    kpc = npp - (n - imax + 1)*(n - imax + 2)/2 + 1;
                    if imax < n {
                        jmax = imax + _hptrf_iamax_1based(n - imax, _hptrf_subvector(ap.*, kpc + 1, n - imax), 1, MT);
                        rowmax = max(rowmax, _hptrf_abs1(fget(ap.*, kpc + jmax - imax), RT, MT));
                    }

                    if absakk >= alpha*colmax*(colmax / rowmax) {
                        kp = k;
                    } else if abs(real(fget(ap.*, kpc))) >= alpha*rowmax {
                        kp = imax;
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                }

                kk = k + kstep - 1;
                if kstep == 2 {
                    knc += n - k + 1;
                }
                if kp != kk {
                    kpc = npp - (n - kp + 1)*(n - kp + 2)/2 + 1;
                    if kp < n {
                        col_kk_tail := _hptrf_subvector(ap.*, knc + kp - kk + 1, n - kp);
                        col_kp_tail := _hptrf_subvector(ap.*, kpc + 1, n - kp);
                        _hptrf_swap(n - kp, *col_kk_tail, 1, *col_kp_tail, 1, MT);
                    }
                    kx = knc + kp - kk;
                    for j:kk + 1..kp - 1 {
                        kx += n - j + 1;
                        t = conjugate(fget(ap.*, knc + j - kk));
                        fset(ap, knc + j - kk, conjugate(fget(ap.*, kx)));
                        fset(ap, kx, t);
                    }
                    fset(ap, knc + kp - kk, conjugate(fget(ap.*, knc + kp - kk)));
                    r1 = real(fget(ap.*, knc));
                    fset(ap, knc, _hptrf_complex_from_real(real(fget(ap.*, kpc)), RT, MT));
                    fset(ap, kpc, _hptrf_complex_from_real(r1, RT, MT));
                    if kstep == 2 {
                        fset(ap, kc, _hptrf_complex_from_real(real(fget(ap.*, kc)), RT, MT));
                        t = fget(ap.*, kc + 1);
                        fset(ap, kc + 1, fget(ap.*, kc + kp - k));
                        fset(ap, kc + kp - k, t);
                    }
                } else {
                    fset(ap, kc, _hptrf_complex_from_real(real(fget(ap.*, kc)), RT, MT));
                    if kstep == 2 {
                        fset(ap, knc, _hptrf_complex_from_real(real(fget(ap.*, knc)), RT, MT));
                    }
                }

                if kstep == 1 {
                    if k < n {
                        r1 = one_r / real(fget(ap.*, kc));
                        col_k := _hptrf_subvector(ap.*, kc + 1, n - k);
                        ap_tail := _hptrf_subvector(ap.*, kc + (n - k + 1), (n - k)*(n - k + 1)/2);
                        _hptrf_hpr(uplo, n - k, -r1, col_k, 1, *ap_tail, RT, MT);
                        _hptrf_scal_real(n - k, r1, *col_k, 1, RT, MT);
                    }
                } else if k < n - 1 {
                    d = _hptrf_lapy2(real(fget(ap.*, kc + 1)), imag(fget(ap.*, kc + 1)), RT);
                    d11 = real(fget(ap.*, knc)) / d;
                    d22 = real(fget(ap.*, kc)) / d;
                    tt = one_r / (d11*d22 - one_r);
                    d21 = fget(ap.*, kc + 1) / d;
                    d = tt / d;

                    for j:k + 2..n {
                        idx_j_k :: _hptrf_lower_index(n, j, k);
                        idx_j_kp1 :: _hptrf_lower_index(n, j, k + 1);
                        wk = d*(d11*fget(ap.*, idx_j_k) - d21*fget(ap.*, idx_j_kp1));
                        wkp1 = d*(d22*fget(ap.*, idx_j_kp1) - conjugate(d21)*fget(ap.*, idx_j_k));
                        for i:j..n {
                            idx_ij :: _hptrf_lower_index(n, i, j);
                            fset(ap, idx_ij, fget(ap.*, idx_ij) - fget(ap.*, _hptrf_lower_index(n, i, k))*conjugate(wk) - fget(ap.*, _hptrf_lower_index(n, i, k + 1))*conjugate(wkp1));
                        }
                        fset(ap, idx_j_k, wk);
                        fset(ap, idx_j_kp1, wkp1);
                        diag_j :: _hptrf_lower_index(n, j, j);
                        fset(ap, diag_j, _hptrf_complex_from_real(real(fget(ap.*, diag_j)), RT, MT));
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
            kc = knc + n - k + 2;
        }
    }
}
