// ? LAPACK *hetrs_rook()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrs_rook/zhetrs_rook
// *HETRS_ROOK solves a system of linear equations A*X = B for Hermitian matrices.

chetrs_rook :: #bake_arguments _hetrs_rook(RT=float32, MT=C32);
zhetrs_rook :: #bake_arguments _hetrs_rook(RT=float64, MT=C64);

#scope_file
_hetrs_rook_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrs_rook_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetrs_rook_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetrs_rook_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _hetrs_rook_raw_slice(a.data + (col - 1)*lda + (row - 1), ifx count < 0 then 0 else count);
}

_hetrs_rook_matrix_row :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    n := ifx count <= 0 then 0 else (count - 1)*lda + 1;
    return _hetrs_rook_raw_slice(a.data + (col - 1)*lda + (row - 1), n);
}

_hetrs_rook_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrs_rook_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hetrs_rook_scal :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_hetrs_rook_geru :: (
    $$m: int, $$n: int, alpha: $MT, x: $X/ArrayType1D, $$incx: int,
    y: $Y/ArrayType1D, $$incy: int, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT && A.data_type == MT;
} {
    #if MT == C32 { cgeru(m, n, alpha, x, incx, y, incy, a, lda); } else { zgeru(m, n, alpha, x, incx, y, incy, a, lda); }
}

_hetrs_rook_gemv :: (
    $$trans: u8, $$m: int, $$n: int, alpha: $MT, a: $A/ArrayType2D, $$lda: int,
    x: $X/ArrayType1D, $$incx: int, beta: $MT, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 {
        cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_hetrs_rook_lacgv :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { clacgv(n, x, incx); } else { zlacgv(n, x, incx); }
}

_hetrs_rook :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: $A/ArrayType2D, $$lda: int,
    ipiv: $IPIV/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    one_t  :: #run _hetrs_rook_complex_from_real(one(RT), RT, MT);
    neg_one_t :: #run _hetrs_rook_complex_from_real(-one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHETRS_ROOK" } else { "ZHETRS_ROOK" };

    upper : bool = ---;
    j, k, kp : int = ---;
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
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ldb < max(1, n) {
        info.* = -8;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 || nrhs == 0 {
        return;
    }

    if upper {
        k = n;
        while true {
            if k < 1 then break;

            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }

                if k > 1 {
                    a_col := _hetrs_rook_matrix_column(a, lda, 1, k, k - 1);
                    b_row := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    b_top_raw := _hetrs_rook_matrix_storage(b.*, ldb, 1, 1, k - 1, nrhs);
                    b_top := array2d_view(b_top_raw, ldb, row_major=false);
                    _hetrs_rook_geru(k - 1, nrhs, neg_one_t, a_col, 1, b_row, ldb, *b_top, ldb, MT);
                }

                s = one(RT) / real(fget(a, k, k));
                row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                _hetrs_rook_scal(nrhs, s, *row_k, ldb, RT, MT);
                k -= 1;
            } else {
                kp = -fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }

                kp = -fget(ipiv, k - 1);
                if kp != k - 1 {
                    row_km1 := _hetrs_rook_matrix_row(b.*, ldb, k - 1, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_km1, ldb, *row_kp, ldb, MT);
                }

                if k > 2 {
                    a_col_k := _hetrs_rook_matrix_column(a, lda, 1, k, k - 2);
                    a_col_km1 := _hetrs_rook_matrix_column(a, lda, 1, k - 1, k - 2);
                    b_row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    b_row_km1 := _hetrs_rook_matrix_row(b.*, ldb, k - 1, 1, nrhs);
                    b_top_raw := _hetrs_rook_matrix_storage(b.*, ldb, 1, 1, k - 2, nrhs);
                    b_top := array2d_view(b_top_raw, ldb, row_major=false);
                    _hetrs_rook_geru(k - 2, nrhs, neg_one_t, a_col_k, 1, b_row_k, ldb, *b_top, ldb, MT);
                    _hetrs_rook_geru(k - 2, nrhs, neg_one_t, a_col_km1, 1, b_row_km1, ldb, *b_top, ldb, MT);
                }

                akm1k = fget(a, k - 1, k);
                akm1 = fget(a, k - 1, k - 1) / akm1k;
                ak = fget(a, k, k) / conjugate(akm1k);
                denom = akm1*ak - one_t;
                for j:1..nrhs {
                    bkm1 = fget(b.*, k - 1, j) / akm1k;
                    bk = fget(b.*, k, j) / conjugate(akm1k);
                    fset(b, k - 1, j, (ak*bkm1 - bk) / denom);
                    fset(b, k, j, (akm1*bk - bkm1) / denom);
                }
                k -= 2;
            }
        }

        k = 1;
        while true {
            if k > n then break;

            if fget(ipiv, k) > 0 {
                if k > 1 {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    b_top_raw := _hetrs_rook_matrix_storage(b.*, ldb, 1, 1, k - 1, nrhs);
                    b_top := array2d_view(b_top_raw, ldb, row_major=false);
                    a_col := _hetrs_rook_matrix_column(a, lda, 1, k, k - 1);
                    _hetrs_rook_lacgv(nrhs, *row_k, ldb, MT);
                    _hetrs_rook_gemv(#char "C", k - 1, nrhs, neg_one_t, *b_top, ldb, a_col, 1, one_t, *row_k, ldb, MT);
                    _hetrs_rook_lacgv(nrhs, *row_k, ldb, MT);
                }

                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k += 1;
            } else {
                if k > 1 {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp1 := _hetrs_rook_matrix_row(b.*, ldb, k + 1, 1, nrhs);
                    b_top_raw := _hetrs_rook_matrix_storage(b.*, ldb, 1, 1, k - 1, nrhs);
                    b_top := array2d_view(b_top_raw, ldb, row_major=false);
                    a_col_k := _hetrs_rook_matrix_column(a, lda, 1, k, k - 1);
                    a_col_kp1 := _hetrs_rook_matrix_column(a, lda, 1, k + 1, k - 1);
                    _hetrs_rook_lacgv(nrhs, *row_k, ldb, MT);
                    _hetrs_rook_gemv(#char "C", k - 1, nrhs, neg_one_t, *b_top, ldb, a_col_k, 1, one_t, *row_k, ldb, MT);
                    _hetrs_rook_lacgv(nrhs, *row_k, ldb, MT);
                    _hetrs_rook_lacgv(nrhs, *row_kp1, ldb, MT);
                    _hetrs_rook_gemv(#char "C", k - 1, nrhs, neg_one_t, *b_top, ldb, a_col_kp1, 1, one_t, *row_kp1, ldb, MT);
                    _hetrs_rook_lacgv(nrhs, *row_kp1, ldb, MT);
                }

                kp = -fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }

                kp = -fget(ipiv, k + 1);
                if kp != k + 1 {
                    row_kp1 := _hetrs_rook_matrix_row(b.*, ldb, k + 1, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_kp1, ldb, *row_kp, ldb, MT);
                }
                k += 2;
            }
        }
    } else {
        k = 1;
        while true {
            if k > n then break;

            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }

                if k < n {
                    a_col := _hetrs_rook_matrix_column(a, lda, k + 1, k, n - k);
                    b_row := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    b_tail_raw := _hetrs_rook_matrix_storage(b.*, ldb, k + 1, 1, n - k, nrhs);
                    b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
                    _hetrs_rook_geru(n - k, nrhs, neg_one_t, a_col, 1, b_row, ldb, *b_tail, ldb, MT);
                }

                s = one(RT) / real(fget(a, k, k));
                row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                _hetrs_rook_scal(nrhs, s, *row_k, ldb, RT, MT);
                k += 1;
            } else {
                kp = -fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }

                kp = -fget(ipiv, k + 1);
                if kp != k + 1 {
                    row_kp1 := _hetrs_rook_matrix_row(b.*, ldb, k + 1, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_kp1, ldb, *row_kp, ldb, MT);
                }

                if k < n - 1 {
                    a_col_k := _hetrs_rook_matrix_column(a, lda, k + 2, k, n - k - 1);
                    a_col_kp1 := _hetrs_rook_matrix_column(a, lda, k + 2, k + 1, n - k - 1);
                    b_row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    b_row_kp1 := _hetrs_rook_matrix_row(b.*, ldb, k + 1, 1, nrhs);
                    b_tail_raw := _hetrs_rook_matrix_storage(b.*, ldb, k + 2, 1, n - k - 1, nrhs);
                    b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
                    _hetrs_rook_geru(n - k - 1, nrhs, neg_one_t, a_col_k, 1, b_row_k, ldb, *b_tail, ldb, MT);
                    _hetrs_rook_geru(n - k - 1, nrhs, neg_one_t, a_col_kp1, 1, b_row_kp1, ldb, *b_tail, ldb, MT);
                }

                akm1k = fget(a, k + 1, k);
                akm1 = fget(a, k, k) / conjugate(akm1k);
                ak = fget(a, k + 1, k + 1) / akm1k;
                denom = akm1*ak - one_t;
                for j:1..nrhs {
                    bkm1 = fget(b.*, k, j) / conjugate(akm1k);
                    bk = fget(b.*, k + 1, j) / akm1k;
                    fset(b, k, j, (ak*bkm1 - bk) / denom);
                    fset(b, k + 1, j, (akm1*bk - bkm1) / denom);
                }
                k += 2;
            }
        }

        k = n;
        while true {
            if k < 1 then break;

            if fget(ipiv, k) > 0 {
                if k < n {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    b_tail_raw := _hetrs_rook_matrix_storage(b.*, ldb, k + 1, 1, n - k, nrhs);
                    b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
                    a_col := _hetrs_rook_matrix_column(a, lda, k + 1, k, n - k);
                    _hetrs_rook_lacgv(nrhs, *row_k, ldb, MT);
                    _hetrs_rook_gemv(#char "C", n - k, nrhs, neg_one_t, *b_tail, ldb, a_col, 1, one_t, *row_k, ldb, MT);
                    _hetrs_rook_lacgv(nrhs, *row_k, ldb, MT);
                }

                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k -= 1;
            } else {
                if k < n {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_km1 := _hetrs_rook_matrix_row(b.*, ldb, k - 1, 1, nrhs);
                    b_tail_raw := _hetrs_rook_matrix_storage(b.*, ldb, k + 1, 1, n - k, nrhs);
                    b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
                    a_col_k := _hetrs_rook_matrix_column(a, lda, k + 1, k, n - k);
                    a_col_km1 := _hetrs_rook_matrix_column(a, lda, k + 1, k - 1, n - k);
                    _hetrs_rook_lacgv(nrhs, *row_k, ldb, MT);
                    _hetrs_rook_gemv(#char "C", n - k, nrhs, neg_one_t, *b_tail, ldb, a_col_k, 1, one_t, *row_k, ldb, MT);
                    _hetrs_rook_lacgv(nrhs, *row_k, ldb, MT);
                    _hetrs_rook_lacgv(nrhs, *row_km1, ldb, MT);
                    _hetrs_rook_gemv(#char "C", n - k, nrhs, neg_one_t, *b_tail, ldb, a_col_km1, 1, one_t, *row_km1, ldb, MT);
                    _hetrs_rook_lacgv(nrhs, *row_km1, ldb, MT);
                }

                kp = -fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_rook_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }

                kp = -fget(ipiv, k - 1);
                if kp != k - 1 {
                    row_km1 := _hetrs_rook_matrix_row(b.*, ldb, k - 1, 1, nrhs);
                    row_kp := _hetrs_rook_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_rook_swap(nrhs, *row_km1, ldb, *row_kp, ldb, MT);
                }
                k -= 2;
            }
        }
    }
}
