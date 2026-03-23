// ? LAPACK *hptrs()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chptrs/zhptrs
// *HPTRS solves a system of linear equations A*X = B with a complex
// Hermitian packed matrix using the factorization from *HPTRF.

chptrs :: #bake_arguments _hptrs(RT=float32, MT=C32);
zhptrs :: #bake_arguments _hptrs(RT=float64, MT=C64);

#scope_file
_hptrs_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hptrs_subvector :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hptrs_raw_slice(x.data + offset - 1, count);
}

_hptrs_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hptrs_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hptrs_matrix_row :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    n := ifx count <= 0 then 0 else (count - 1)*lda + 1;
    return _hptrs_raw_slice(a.data + (col - 1)*lda + (row - 1), n);
}

_hptrs_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hptrs_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hptrs_scal :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_hptrs_geru :: (
    $$m: int, $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int,
    y: $Y/ArrayType1D, $$incy: int, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT && A.data_type == MT;
} {
    #if MT == C32 { cgeru(m, n, alpha, x, incx, y, incy, a, lda); } else { zgeru(m, n, alpha, x, incx, y, incy, a, lda); }
}

_hptrs_gemv :: (
    $$trans: u8, $$m: int, $$n: int, alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    x: $X/ArrayType1D, $$incx: int, beta: $MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); } else { zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }
}

_hptrs_lacgv :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { clacgv(n, x, incx); } else { zlacgv(n, x, incx); }
}

_hptrs :: (
    $$uplo: u8, $$n: int, $$nrhs: int, ap: $AP/ArrayType1D, ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           AP.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    one_t :: #run _hptrs_complex_from_real(one(RT), RT, MT);
    neg_one_t :: #run -one(MT);
    routine_name :: #run if MT == C32 { "CHPTRS" } else { "ZHPTRS" };

    upper : bool = ---;
    j, k, kc, kp : int = ---;
    s : RT = ---;
    ak, akm1, akm1k, bk, bkm1, denom : MT = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if ldb < max(1, n) {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 || nrhs == 0 then return;

    if upper {
        k = n;
        kc = n*(n + 1)/2 + 1;
        while true {
            if k < 1 then break;

            kc -= k;
            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hptrs_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hptrs_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }

                if k > 1 {
                    ucol := _hptrs_subvector(ap, kc, k - 1);
                    brow := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    btop_raw := _hptrs_matrix_storage(b.*, ldb, 1, 1, k - 1, nrhs);
                    btop := array2d_view(btop_raw, ldb, row_major=false);
                    _hptrs_geru(k - 1, nrhs, neg_one_t, ucol, 1, brow, ldb, *btop, ldb, MT);
                }

                s = one(RT) / real(fget(ap, kc + k - 1));
                row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                _hptrs_scal(nrhs, s, *row_k, ldb, RT, MT);
                k -= 1;
            } else {
                kp = -fget(ipiv, k);
                if kp != k - 1 {
                    row_km1 := _hptrs_matrix_row(b.*, ldb, k - 1, 1, nrhs);
                    row_kp := _hptrs_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hptrs_swap(nrhs, *row_km1, ldb, *row_kp, ldb, MT);
                }

                if k > 2 {
                    col_k := _hptrs_subvector(ap, kc, k - 2);
                    col_km1 := _hptrs_subvector(ap, kc - (k - 1), k - 2);
                    brow_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    brow_km1 := _hptrs_matrix_row(b.*, ldb, k - 1, 1, nrhs);
                    btop_raw := _hptrs_matrix_storage(b.*, ldb, 1, 1, k - 2, nrhs);
                    btop := array2d_view(btop_raw, ldb, row_major=false);
                    _hptrs_geru(k - 2, nrhs, neg_one_t, col_k, 1, brow_k, ldb, *btop, ldb, MT);
                    _hptrs_geru(k - 2, nrhs, neg_one_t, col_km1, 1, brow_km1, ldb, *btop, ldb, MT);
                }

                akm1k = fget(ap, kc + k - 2);
                akm1 = fget(ap, kc - 1) / akm1k;
                ak = fget(ap, kc + k - 1) / conjugate(akm1k);
                denom = akm1*ak - one_t;
                for j:1..nrhs {
                    bkm1 = fget(b.*, k - 1, j) / akm1k;
                    bk = fget(b.*, k, j) / conjugate(akm1k);
                    fset(b, k - 1, j, (ak*bkm1 - bk) / denom);
                    fset(b, k, j, (akm1*bk - bkm1) / denom);
                }
                kc = kc - k + 1;
                k -= 2;
            }
        }

        k = 1;
        kc = 1;
        while true {
            if k > n then break;

            if fget(ipiv, k) > 0 {
                if k > 1 {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    btop_raw := _hptrs_matrix_storage(b.*, ldb, 1, 1, k - 1, nrhs);
                    btop := array2d_view(btop_raw, ldb, row_major=false);
                    ucol := _hptrs_subvector(ap, kc, k - 1);
                    _hptrs_lacgv(nrhs, *row_k, ldb, MT);
                    _hptrs_gemv(#char "C", k - 1, nrhs, neg_one_t, *btop, ldb, ucol, 1, one_t, *row_k, ldb, MT);
                    _hptrs_lacgv(nrhs, *row_k, ldb, MT);
                }

                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hptrs_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hptrs_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                kc += k;
                k += 1;
            } else {
                if k > 1 {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp1 := _hptrs_matrix_row(b.*, ldb, k + 1, 1, nrhs);
                    btop_raw := _hptrs_matrix_storage(b.*, ldb, 1, 1, k - 1, nrhs);
                    btop := array2d_view(btop_raw, ldb, row_major=false);
                    col_k := _hptrs_subvector(ap, kc, k - 1);
                    col_kp1 := _hptrs_subvector(ap, kc + k, k - 1);
                    _hptrs_lacgv(nrhs, *row_k, ldb, MT);
                    _hptrs_gemv(#char "C", k - 1, nrhs, neg_one_t, *btop, ldb, col_k, 1, one_t, *row_k, ldb, MT);
                    _hptrs_lacgv(nrhs, *row_k, ldb, MT);
                    _hptrs_lacgv(nrhs, *row_kp1, ldb, MT);
                    _hptrs_gemv(#char "C", k - 1, nrhs, neg_one_t, *btop, ldb, col_kp1, 1, one_t, *row_kp1, ldb, MT);
                    _hptrs_lacgv(nrhs, *row_kp1, ldb, MT);
                }

                kp = -fget(ipiv, k);
                if kp != k {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hptrs_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hptrs_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                kc += 2*k + 1;
                k += 2;
            }
        }
    } else {
        k = 1;
        kc = 1;
        while true {
            if k > n then break;

            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hptrs_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hptrs_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }

                if k < n {
                    lcol := _hptrs_subvector(ap, kc + 1, n - k);
                    brow := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    btail_raw := _hptrs_matrix_storage(b.*, ldb, k + 1, 1, n - k, nrhs);
                    btail := array2d_view(btail_raw, ldb, row_major=false);
                    _hptrs_geru(n - k, nrhs, neg_one_t, lcol, 1, brow, ldb, *btail, ldb, MT);
                }

                s = one(RT) / real(fget(ap, kc));
                row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                _hptrs_scal(nrhs, s, *row_k, ldb, RT, MT);
                kc += n - k + 1;
                k += 1;
            } else {
                kp = -fget(ipiv, k);
                if kp != k + 1 {
                    row_kp1 := _hptrs_matrix_row(b.*, ldb, k + 1, 1, nrhs);
                    row_kp := _hptrs_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hptrs_swap(nrhs, *row_kp1, ldb, *row_kp, ldb, MT);
                }

                if k < n - 1 {
                    col_k := _hptrs_subvector(ap, kc + 2, n - k - 1);
                    col_kp1 := _hptrs_subvector(ap, kc + (n - k + 2), n - k - 1);
                    brow_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    brow_kp1 := _hptrs_matrix_row(b.*, ldb, k + 1, 1, nrhs);
                    btail_raw := _hptrs_matrix_storage(b.*, ldb, k + 2, 1, n - k - 1, nrhs);
                    btail := array2d_view(btail_raw, ldb, row_major=false);
                    _hptrs_geru(n - k - 1, nrhs, neg_one_t, col_k, 1, brow_k, ldb, *btail, ldb, MT);
                    _hptrs_geru(n - k - 1, nrhs, neg_one_t, col_kp1, 1, brow_kp1, ldb, *btail, ldb, MT);
                }

                akm1k = fget(ap, kc + 1);
                akm1 = fget(ap, kc) / conjugate(akm1k);
                ak = fget(ap, kc + n - k + 1) / akm1k;
                denom = akm1*ak - one_t;
                for j:1..nrhs {
                    bkm1 = fget(b.*, k, j) / conjugate(akm1k);
                    bk = fget(b.*, k + 1, j) / akm1k;
                    fset(b, k, j, (ak*bkm1 - bk) / denom);
                    fset(b, k + 1, j, (akm1*bk - bkm1) / denom);
                }
                kc += 2*(n - k) + 1;
                k += 2;
            }
        }

        k = n;
        kc = n*(n + 1)/2 + 1;
        while true {
            if k < 1 then break;

            kc -= n - k + 1;
            if fget(ipiv, k) > 0 {
                if k < n {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    btail_raw := _hptrs_matrix_storage(b.*, ldb, k + 1, 1, n - k, nrhs);
                    btail := array2d_view(btail_raw, ldb, row_major=false);
                    lcol := _hptrs_subvector(ap, kc + 1, n - k);
                    _hptrs_lacgv(nrhs, *row_k, ldb, MT);
                    _hptrs_gemv(#char "C", n - k, nrhs, neg_one_t, *btail, ldb, lcol, 1, one_t, *row_k, ldb, MT);
                    _hptrs_lacgv(nrhs, *row_k, ldb, MT);
                }

                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hptrs_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hptrs_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k -= 1;
            } else {
                if k < n {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_km1 := _hptrs_matrix_row(b.*, ldb, k - 1, 1, nrhs);
                    btail_raw := _hptrs_matrix_storage(b.*, ldb, k + 1, 1, n - k, nrhs);
                    btail := array2d_view(btail_raw, ldb, row_major=false);
                    col_k := _hptrs_subvector(ap, kc + 1, n - k);
                    col_km1 := _hptrs_subvector(ap, kc - (n - k), n - k);
                    _hptrs_lacgv(nrhs, *row_k, ldb, MT);
                    _hptrs_gemv(#char "C", n - k, nrhs, neg_one_t, *btail, ldb, col_k, 1, one_t, *row_k, ldb, MT);
                    _hptrs_lacgv(nrhs, *row_k, ldb, MT);
                    _hptrs_lacgv(nrhs, *row_km1, ldb, MT);
                    _hptrs_gemv(#char "C", n - k, nrhs, neg_one_t, *btail, ldb, col_km1, 1, one_t, *row_km1, ldb, MT);
                    _hptrs_lacgv(nrhs, *row_km1, ldb, MT);
                }

                kp = -fget(ipiv, k);
                if kp != k {
                    row_k := _hptrs_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hptrs_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hptrs_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                kc -= n - k + 2;
                k -= 2;
            }
        }
    }
}
