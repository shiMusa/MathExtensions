// ? LAPACK *hetrs2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrs2/zhetrs2
// *HETRS2 solves a system of linear equations A*X = B for Hermitian matrices.

chetrs2 :: #bake_arguments _hetrs2(RT=float32, MT=C32);
zhetrs2 :: #bake_arguments _hetrs2(RT=float64, MT=C64);

#scope_file
_hetrs2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrs2_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetrs2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetrs2_matrix_row :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    n := ifx count <= 0 then 0 else (count - 1)*lda + 1;
    return _hetrs2_raw_slice(a.data + (col - 1)*lda + (row - 1), n);
}

_hetrs2_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetrs2_raw_slice(x.data + (offset - 1), count);
}

_hetrs2_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrs2_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hetrs2_scal :: (
    $$n: int, alpha: $RT, x: *$X/ArrayType1D, $$incx: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT) && X.data_type == MT;
} {
    #if MT == C32 { csscal(n, alpha, x, incx); } else { zdscal(n, alpha, x, incx); }
}

_hetrs2_trsm :: (
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

_hetrs2_syconv :: (
    $$uplo: u8, $$way: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: $IPIV/ArrayType1D, work: *$WORK/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    #if MT == C32 {
        csyconv(uplo, way, n, a, lda, ipiv, work, info);
    } else {
        zsyconv(uplo, way, n, a, lda, ipiv, work, info);
    }
}

_hetrs2 :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: $IPIV/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int,
    work: *$WORK/ArrayType1D, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && WORK.data_type == MT;
} {
    one_t :: #run _hetrs2_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHETRS2" } else { "ZHETRS2" };

    upper : bool = ---;
    i, iinfo, j, k, kp : int = ---;
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

    _hetrs2_syconv(uplo, #char "C", n, a, lda, ipiv, work, *iinfo, MT);

    a_full_raw := _hetrs2_matrix_storage(a.*, lda, 1, 1, n, n);
    a_full := array2d_view(a_full_raw, lda, row_major=false);
    b_full_raw := _hetrs2_matrix_storage(b.*, ldb, 1, 1, n, nrhs);
    b_full := array2d_view(b_full_raw, ldb, row_major=false);

    if upper {
        k = n;
        while true {
            if k < 1 then break;
            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs2_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs2_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs2_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k -= 1;
            } else {
                kp = -fget(ipiv, k);
                if kp == -fget(ipiv, k - 1) {
                    row_km1 := _hetrs2_matrix_row(b.*, ldb, k - 1, 1, nrhs);
                    row_kp := _hetrs2_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs2_swap(nrhs, *row_km1, ldb, *row_kp, ldb, MT);
                }
                k -= 2;
            }
        }

        _hetrs2_trsm(#char "L", #char "U", #char "N", #char "U", n, nrhs, one_t, *a_full, lda, *b_full, ldb, MT);

        i = n;
        while true {
            if i < 1 then break;
            if fget(ipiv, i) > 0 {
                s = one(RT) / real(fget(a.*, i, i));
                row_i := _hetrs2_matrix_row(b.*, ldb, i, 1, nrhs);
                _hetrs2_scal(nrhs, s, *row_i, ldb, RT, MT);
            } else if i > 1 {
                if fget(ipiv, i - 1) == fget(ipiv, i) {
                    akm1k = fget(work, i);
                    akm1 = fget(a.*, i - 1, i - 1) / akm1k;
                    ak = fget(a.*, i, i) / conjugate(akm1k);
                    denom = akm1*ak - one_t;
                    for j:1..nrhs {
                        bkm1 = fget(b.*, i - 1, j) / akm1k;
                        bk = fget(b.*, i, j) / conjugate(akm1k);
                        fset(b, i - 1, j, (ak*bkm1 - bk) / denom);
                        fset(b, i, j, (akm1*bk - bkm1) / denom);
                    }
                    i -= 1;
                }
            }
            i -= 1;
        }

        _hetrs2_trsm(#char "L", #char "U", #char "C", #char "U", n, nrhs, one_t, *a_full, lda, *b_full, ldb, MT);

        k = 1;
        while true {
            if k > n then break;
            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs2_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs2_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs2_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k += 1;
            } else {
                kp = -fget(ipiv, k);
                if k < n && kp == -fget(ipiv, k + 1) {
                    row_k := _hetrs2_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs2_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs2_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
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
                    row_k := _hetrs2_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs2_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs2_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k += 1;
            } else {
                kp = -fget(ipiv, k + 1);
                if kp == -fget(ipiv, k) {
                    row_kp1 := _hetrs2_matrix_row(b.*, ldb, k + 1, 1, nrhs);
                    row_kp := _hetrs2_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs2_swap(nrhs, *row_kp1, ldb, *row_kp, ldb, MT);
                }
                k += 2;
            }
        }

        _hetrs2_trsm(#char "L", #char "L", #char "N", #char "U", n, nrhs, one_t, *a_full, lda, *b_full, ldb, MT);

        i = 1;
        while true {
            if i > n then break;
            if fget(ipiv, i) > 0 {
                s = one(RT) / real(fget(a.*, i, i));
                row_i := _hetrs2_matrix_row(b.*, ldb, i, 1, nrhs);
                _hetrs2_scal(nrhs, s, *row_i, ldb, RT, MT);
            } else {
                akm1k = fget(work, i);
                akm1 = fget(a.*, i, i) / conjugate(akm1k);
                ak = fget(a.*, i + 1, i + 1) / akm1k;
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

        _hetrs2_trsm(#char "L", #char "L", #char "C", #char "U", n, nrhs, one_t, *a_full, lda, *b_full, ldb, MT);

        k = n;
        while true {
            if k < 1 then break;
            if fget(ipiv, k) > 0 {
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs2_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs2_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs2_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k -= 1;
            } else {
                kp = -fget(ipiv, k);
                if k > 1 && kp == -fget(ipiv, k - 1) {
                    row_k := _hetrs2_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs2_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs2_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k -= 2;
            }
        }
    }

    _hetrs2_syconv(uplo, #char "R", n, a, lda, ipiv, work, *iinfo, MT);
}
