// ? LAPACK *hetrs_3()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrs_3/zhetrs_3
// *HETRS_3 solves a system of linear equations A*X = B for Hermitian matrices.

chetrs_3 :: #bake_arguments _hetrs_3(RT=float32, MT=C32);
zhetrs_3 :: #bake_arguments _hetrs_3(RT=float64, MT=C64);

#scope_file
_hetrs_3_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrs_3_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetrs_3_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetrs_3_matrix_row :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    n := ifx count <= 0 then 0 else (count - 1)*lda + 1;
    return _hetrs_3_raw_slice(a.data + (col - 1)*lda + (row - 1), n);
}

_hetrs_3_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrs_3_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hetrs_3_scal :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_hetrs_3_trsm :: (
    $$side: u8, $$uplo: u8, $$transa: u8, $$diag: u8, $$m: int, $$n: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_hetrs_3 :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: $A/ArrayType2D, $$lda: int,
    e: $E/ArrayType1D, ipiv: $IPIV/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && E.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    one_t :: #run _hetrs_3_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHETRS_3" } else { "ZHETRS_3" };

    upper : bool = ---;
    i, j, k, kp : int = ---;
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
        info.* = -9;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 || nrhs == 0 {
        return;
    }

    b_full_raw := _hetrs_3_matrix_storage(b.*, ldb, 1, 1, n, nrhs);
    b_full := array2d_view(b_full_raw, ldb, row_major=false);
    a_full_raw := _hetrs_3_matrix_storage(a, lda, 1, 1, n, n);
    a_full := array2d_view(a_full_raw, lda, row_major=false);

    if upper {
        for <k:n..1 {
            kp = abs(fget(ipiv, k));
            if kp != k {
                row_k := _hetrs_3_matrix_row(b.*, ldb, k, 1, nrhs);
                row_kp := _hetrs_3_matrix_row(b.*, ldb, kp, 1, nrhs);
                _hetrs_3_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
            }
        }

        _hetrs_3_trsm(#char "L", #char "U", #char "N", #char "U", n, nrhs, one_t, *a_full, lda, *b_full, ldb, MT);

        i = n;
        while true {
            if i < 1 then break;
            if fget(ipiv, i) > 0 {
                s = one(RT) / real(fget(a, i, i));
                row_i := _hetrs_3_matrix_row(b.*, ldb, i, 1, nrhs);
                _hetrs_3_scal(nrhs, s, *row_i, ldb, RT, MT);
            } else if i > 1 {
                akm1k = fget(e, i);
                akm1 = fget(a, i - 1, i - 1) / akm1k;
                ak = fget(a, i, i) / conjugate(akm1k);
                denom = akm1*ak - one_t;
                for j:1..nrhs {
                    bkm1 = fget(b.*, i - 1, j) / akm1k;
                    bk = fget(b.*, i, j) / conjugate(akm1k);
                    fset(b, i - 1, j, (ak*bkm1 - bk) / denom);
                    fset(b, i, j, (akm1*bk - bkm1) / denom);
                }
                i -= 1;
            }
            i -= 1;
        }

        _hetrs_3_trsm(#char "L", #char "U", #char "C", #char "U", n, nrhs, one_t, *a_full, lda, *b_full, ldb, MT);

        for k:1..n {
            kp = abs(fget(ipiv, k));
            if kp != k {
                row_k := _hetrs_3_matrix_row(b.*, ldb, k, 1, nrhs);
                row_kp := _hetrs_3_matrix_row(b.*, ldb, kp, 1, nrhs);
                _hetrs_3_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
            }
        }
    } else {
        for k:1..n {
            kp = abs(fget(ipiv, k));
            if kp != k {
                row_k := _hetrs_3_matrix_row(b.*, ldb, k, 1, nrhs);
                row_kp := _hetrs_3_matrix_row(b.*, ldb, kp, 1, nrhs);
                _hetrs_3_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
            }
        }

        _hetrs_3_trsm(#char "L", #char "L", #char "N", #char "U", n, nrhs, one_t, *a_full, lda, *b_full, ldb, MT);

        i = 1;
        while true {
            if i > n then break;
            if fget(ipiv, i) > 0 {
                s = one(RT) / real(fget(a, i, i));
                row_i := _hetrs_3_matrix_row(b.*, ldb, i, 1, nrhs);
                _hetrs_3_scal(nrhs, s, *row_i, ldb, RT, MT);
            } else if i < n {
                akm1k = fget(e, i);
                akm1 = fget(a, i, i) / conjugate(akm1k);
                ak = fget(a, i + 1, i + 1) / akm1k;
                denom = akm1*ak - one_t;
                for j:1..nrhs {
                    bkm1 = fget(b.*, i, j) / conjugate(akm1k);
                    bk = fget(b.*, i + 1, j) / akm1k;
                    fset(b, i, j, (ak*bkm1 - bk) / denom);
                    fset(b, i + 1, j, (akm1*bk - bkm1) / denom);
                }
                i += 1;
            }
            i += 1;
        }

        _hetrs_3_trsm(#char "L", #char "L", #char "C", #char "U", n, nrhs, one_t, *a_full, lda, *b_full, ldb, MT);

        for <k:n..1 {
            kp = abs(fget(ipiv, k));
            if kp != k {
                row_k := _hetrs_3_matrix_row(b.*, ldb, k, 1, nrhs);
                row_kp := _hetrs_3_matrix_row(b.*, ldb, kp, 1, nrhs);
                _hetrs_3_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
            }
        }
    }
}
