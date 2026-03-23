// ? LAPACK *hptri()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chptri/zhptri
// *HPTRI computes the inverse of a complex Hermitian indefinite matrix
// in packed storage using the factorization computed by *HPTRF.

chptri :: #bake_arguments _hptri(RT=float32, MT=C32);
zhptri :: #bake_arguments _hptri(RT=float64, MT=C64);

#scope_file
_hptri_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hptri_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hptri_raw_slice(x.data + offset - 1, count);
}

_hptri_upper_col_start :: inline ($$j: int) -> int {
    return 1 + (j - 1)*j/2;
}

_hptri_lower_col_start :: inline ($$n: int, $$j: int) -> int {
    return 1 + (j - 1)*(2*n - j + 2)/2;
}

_hptri_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hptri_copy :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { ccopy(n, x, incx, y, incy); } else { zcopy(n, x, incx, y, incy); }
}

_hptri_hpmv :: (
    $$uplo: u8, $$n: int, alpha: $MT, ap: $AP/ArrayType1D, x: $X/ArrayType1D, $$incx: int,
    beta: MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AP.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy); } else { zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy); }
}

_hptri_dotc :: (
    $$n: int, x: $X/ArrayType1D, $$incx: int, y: $Y/ArrayType1D, $$incy: int, $MT: Type
) -> MT #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { return cdotc(n, x, incx, y, incy); } else { return zdotc(n, x, incx, y, incy); }
}

_hptri_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hptri :: (
    $$uplo: u8, $$n: int, ap: *$AP/ArrayType1D, ipiv: $IPIV/ArrayType1D,
    work: *$WORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    one_r  :: #run one(RT);
    one_t  :: #run _hptri_complex_from_real(one_r, RT, MT);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == C32 { "CHPTRI" } else { "ZHPTRI" };

    upper : bool = ---;
    j, k, kc, kcnext, kp, kpc, kstep, kx, npp : int = ---;
    ak, akp1, d, t, r1 : RT = ---;
    akkp1, temp : MT = ---;

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

    if n == 0 then return;

    if upper {
        kp = n*(n + 1)/2;
        for <j:1..n {
            if fget(ipiv, j) > 0 && fget(ap.*, kp) == zero_t {
                info.* = j;
                return;
            }
            kp -= j;
        }
    } else {
        kp = 1;
        for j:1..n {
            if fget(ipiv, j) > 0 && fget(ap.*, kp) == zero_t {
                info.* = j;
                return;
            }
            kp += n - j + 1;
        }
    }
    info.* = 0;

    if upper {
        k = 1;
        kc = 1;
        while true {
            if k > n then break;

            kcnext = kc + k;
            if fget(ipiv, k) > 0 {
                fset(ap, kc + k - 1, _hptri_complex_from_real(one_r / real(fget(ap.*, kc + k - 1)), RT, MT));
                if k > 1 {
                    work_head := _hptri_subvector(work.*, 1, k - 1);
                    col_k := _hptri_subvector(ap.*, kc, k - 1);
                    ap_head := _hptri_subvector(ap.*, 1, (k - 1)*k/2);
                    _hptri_copy(k - 1, col_k, 1, *work_head, 1, MT);
                    _hptri_hpmv(uplo, k - 1, -one_t, ap_head, work_head, 1, zero_t, *col_k, 1, MT);
                    fset(ap, kc + k - 1, fget(ap.*, kc + k - 1) - _hptri_complex_from_real(real(_hptri_dotc(k - 1, work_head, 1, col_k, 1, MT)), RT, MT));
                }
                kstep = 1;
            } else {
                t = abs(fget(ap.*, kcnext + k - 1));
                ak = real(fget(ap.*, kc + k - 1)) / t;
                akp1 = real(fget(ap.*, kcnext + k)) / t;
                akkp1 = fget(ap.*, kcnext + k - 1) / t;
                d = t*(ak*akp1 - one_r);
                fset(ap, kc + k - 1, _hptri_complex_from_real(akp1 / d, RT, MT));
                fset(ap, kcnext + k, _hptri_complex_from_real(ak / d, RT, MT));
                fset(ap, kcnext + k - 1, -akkp1 / d);

                if k > 1 {
                    work_head := _hptri_subvector(work.*, 1, k - 1);
                    col_k := _hptri_subvector(ap.*, kc, k - 1);
                    col_kp1 := _hptri_subvector(ap.*, kcnext, k - 1);
                    ap_head := _hptri_subvector(ap.*, 1, (k - 1)*k/2);

                    _hptri_copy(k - 1, col_k, 1, *work_head, 1, MT);
                    _hptri_hpmv(uplo, k - 1, -one_t, ap_head, work_head, 1, zero_t, *col_k, 1, MT);
                    fset(ap, kc + k - 1, fget(ap.*, kc + k - 1) - _hptri_complex_from_real(real(_hptri_dotc(k - 1, work_head, 1, col_k, 1, MT)), RT, MT));
                    fset(ap, kcnext + k - 1, fget(ap.*, kcnext + k - 1) - _hptri_dotc(k - 1, col_k, 1, col_kp1, 1, MT));

                    _hptri_copy(k - 1, col_kp1, 1, *work_head, 1, MT);
                    _hptri_hpmv(uplo, k - 1, -one_t, ap_head, work_head, 1, zero_t, *col_kp1, 1, MT);
                    fset(ap, kcnext + k, fget(ap.*, kcnext + k) - _hptri_complex_from_real(real(_hptri_dotc(k - 1, work_head, 1, col_kp1, 1, MT)), RT, MT));
                }

                kstep = 2;
                kcnext += k + 1;
            }

            kp = abs(fget(ipiv, k));
            if kp != k {
                kpc = _hptri_upper_col_start(kp);
                if kp > 1 {
                    col_k_pre := _hptri_subvector(ap.*, kc, kp - 1);
                    col_kp_pre := _hptri_subvector(ap.*, kpc, kp - 1);
                    _hptri_swap(kp - 1, *col_k_pre, 1, *col_kp_pre, 1, MT);
                }
                kx = kpc + kp - 1;
                for j:kp + 1..k - 1 {
                    kx += j - 1;
                    temp = conjugate(fget(ap.*, kc + j - 1));
                    fset(ap, kc + j - 1, conjugate(fget(ap.*, kx)));
                    fset(ap, kx, temp);
                }
                fset(ap, kc + kp - 1, conjugate(fget(ap.*, kc + kp - 1)));
                temp = fget(ap.*, kc + k - 1);
                fset(ap, kc + k - 1, fget(ap.*, kpc + kp - 1));
                fset(ap, kpc + kp - 1, temp);
                if kstep == 2 {
                    temp = fget(ap.*, kc + 2*k - 1);
                    fset(ap, kc + 2*k - 1, fget(ap.*, kc + k + kp - 1));
                    fset(ap, kc + k + kp - 1, temp);
                }
            }

            k += kstep;
            kc = kcnext;
        }
    } else {
        npp = n*(n + 1)/2;
        k = n;
        kc = npp;
        while true {
            if k < 1 then break;

            kcnext = kc - (n - k + 2);
            if fget(ipiv, k) > 0 {
                fset(ap, kc, _hptri_complex_from_real(one_r / real(fget(ap.*, kc)), RT, MT));
                if k < n {
                    work_tail := _hptri_subvector(work.*, 1, n - k);
                    col_k_tail := _hptri_subvector(ap.*, kc + 1, n - k);
                    ap_tail := _hptri_subvector(ap.*, kc + (n - k + 1), (n - k)*(n - k + 1)/2);
                    _hptri_copy(n - k, col_k_tail, 1, *work_tail, 1, MT);
                    _hptri_hpmv(uplo, n - k, -one_t, ap_tail, work_tail, 1, zero_t, *col_k_tail, 1, MT);
                    fset(ap, kc, fget(ap.*, kc) - _hptri_complex_from_real(real(_hptri_dotc(n - k, work_tail, 1, col_k_tail, 1, MT)), RT, MT));
                }
                kstep = 1;
            } else {
                t = abs(fget(ap.*, kcnext + 1));
                ak = real(fget(ap.*, kcnext)) / t;
                akp1 = real(fget(ap.*, kc)) / t;
                akkp1 = fget(ap.*, kcnext + 1) / t;
                d = t*(ak*akp1 - one_r);
                fset(ap, kcnext, _hptri_complex_from_real(akp1 / d, RT, MT));
                fset(ap, kc, _hptri_complex_from_real(ak / d, RT, MT));
                fset(ap, kcnext + 1, -akkp1 / d);

                if k < n {
                    work_tail := _hptri_subvector(work.*, 1, n - k);
                    col_k_tail := _hptri_subvector(ap.*, kc + 1, n - k);
                    col_km1_tail := _hptri_subvector(ap.*, kcnext + 2, n - k);
                    ap_tail := _hptri_subvector(ap.*, kc + (n - k + 1), (n - k)*(n - k + 1)/2);

                    _hptri_copy(n - k, col_k_tail, 1, *work_tail, 1, MT);
                    _hptri_hpmv(uplo, n - k, -one_t, ap_tail, work_tail, 1, zero_t, *col_k_tail, 1, MT);
                    fset(ap, kc, fget(ap.*, kc) - _hptri_complex_from_real(real(_hptri_dotc(n - k, work_tail, 1, col_k_tail, 1, MT)), RT, MT));
                    fset(ap, kcnext + 1, fget(ap.*, kcnext + 1) - _hptri_dotc(n - k, col_k_tail, 1, col_km1_tail, 1, MT));

                    _hptri_copy(n - k, col_km1_tail, 1, *work_tail, 1, MT);
                    _hptri_hpmv(uplo, n - k, -one_t, ap_tail, work_tail, 1, zero_t, *col_km1_tail, 1, MT);
                    fset(ap, kcnext, fget(ap.*, kcnext) - _hptri_complex_from_real(real(_hptri_dotc(n - k, work_tail, 1, col_km1_tail, 1, MT)), RT, MT));
                }

                kstep = 2;
                kcnext -= n - k + 3;
            }

            kp = abs(fget(ipiv, k));
            if kp != k {
                kpc = _hptri_lower_col_start(n, kp);
                if kp < n {
                    col_k_tail := _hptri_subvector(ap.*, kc + kp - k + 1, n - kp);
                    col_kp_tail := _hptri_subvector(ap.*, kpc + 1, n - kp);
                    _hptri_swap(n - kp, *col_k_tail, 1, *col_kp_tail, 1, MT);
                }
                kx = kc + kp - k;
                for j:k + 1..kp - 1 {
                    kx += n - j + 1;
                    temp = conjugate(fget(ap.*, kc + j - k));
                    fset(ap, kc + j - k, conjugate(fget(ap.*, kx)));
                    fset(ap, kx, temp);
                }
                fset(ap, kc + kp - k, conjugate(fget(ap.*, kc + kp - k)));
                temp = fget(ap.*, kc);
                fset(ap, kc, fget(ap.*, kpc));
                fset(ap, kpc, temp);
                if kstep == 2 {
                    temp = fget(ap.*, kc - n + k - 1);
                    fset(ap, kc - n + k - 1, fget(ap.*, kc - n + kp - 1));
                    fset(ap, kc - n + kp - 1, temp);
                }
            }

            k -= kstep;
            kc = kcnext;
        }
    }
}
