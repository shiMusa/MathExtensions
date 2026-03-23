// ? LAPACK *hetrs_aa()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrs_aa/zhetrs_aa
// *HETRS_AA solves a system of linear equations A*X = B for Hermitian matrices.

chetrs_aa :: #bake_arguments _hetrs_aa(RT=float32, MT=C32);
zhetrs_aa :: #bake_arguments _hetrs_aa(RT=float64, MT=C64);

#scope_file
_hetrs_aa_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrs_aa_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetrs_aa_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetrs_aa_matrix_row :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    n := ifx count <= 0 then 0 else (count - 1)*lda + 1;
    return _hetrs_aa_raw_slice(a.data + (col - 1)*lda + (row - 1), n);
}

_hetrs_aa_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetrs_aa_raw_slice(x.data + (offset - 1), count);
}

_hetrs_aa_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrs_aa_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetrs_aa_swap :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, y: *$Y/ArrayType1D, $$incy: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT && Y.data_type == MT;
} {
    #if MT == C32 { cswap(n, x, incx, y, incy); } else { zswap(n, x, incx, y, incy); }
}

_hetrs_aa_trsm :: (
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

_hetrs_aa_lacgv :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 { clacgv(n, x, incx); } else { zlacgv(n, x, incx); }
}

_hetrs_aa_gtsv :: (
    $$n: int, $$nrhs: int, dl: *$DL/ArrayType1D, d: *$D/ArrayType1D, du: *$DU/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           DL.data_type == MT && D.data_type == MT && DU.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        cgtsv(n, nrhs, dl, d, du, b, ldb, info);
    } else {
        zgtsv(n, nrhs, dl, d, du, b, ldb, info);
    }
}

_hetrs_aa :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: $A/ArrayType2D, $$lda: int,
    ipiv: $IPIV/ArrayType1D, b: *$B/ArrayType2D, $$ldb: int,
    work: *$WORK/ArrayType1D, $$lwork: int, info: *int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && IPIV.data_type == int &&
           B.data_type == MT && WORK.data_type == MT;
} {
    one_t :: #run _hetrs_aa_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHETRS_AA" } else { "ZHETRS_AA" };

    lquery, upper : bool = ---;
    k, kp, lwkmin : int = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;
    if min(n, nrhs) == 0 {
        lwkmin = 1;
    } else {
        lwkmin = 3*n - 2;
    }

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
    } else if lwork < lwkmin && !lquery {
        info.* = -10;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        fset(work, 1, _hetrs_aa_complex_from_real(_hetrs_aa_roundup_lwork(ncast(RT, lwkmin), RT), RT, MT));
        return;
    }

    if min(n, nrhs) == 0 {
        return;
    }

    if upper {
        if n > 1 {
            k = 1;
            while true {
                if k > n then break;
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_aa_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_aa_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_aa_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k += 1;
            }

            a_u_raw := _hetrs_aa_matrix_storage(a, lda, 1, 2, n - 1, n - 1);
            a_u := array2d_view(a_u_raw, lda, row_major=false);
            b_tail_raw := _hetrs_aa_matrix_storage(b.*, ldb, 2, 1, n - 1, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _hetrs_aa_trsm(#char "L", #char "U", #char "C", #char "U", n - 1, nrhs, one_t, *a_u, lda, *b_tail, ldb, MT);
        }

        for k:1..n {
            fset(work, n + k - 1, fget(a, k, k));
        }
        if n > 1 {
            for k:1..n - 1 {
                fset(work, 2*n + k - 1, fget(a, k, k + 1));
                fset(work, k, fget(a, k, k + 1));
            }
            work_dl := _hetrs_aa_vector_slice(work, 1, n - 1);
            _hetrs_aa_lacgv(n - 1, *work_dl, 1, MT);
        }
        work_dl := _hetrs_aa_vector_slice(work, 1, max(0, n - 1));
        work_d := _hetrs_aa_vector_slice(work, n, n);
        work_du := _hetrs_aa_vector_slice(work, 2*n, max(0, n - 1));
        b_full_raw := _hetrs_aa_matrix_storage(b.*, ldb, 1, 1, n, nrhs);
        b_full := array2d_view(b_full_raw, ldb, row_major=false);
        _hetrs_aa_gtsv(n, nrhs, *work_dl, *work_d, *work_du, *b_full, ldb, info, MT);

        if n > 1 {
            a_u_raw := _hetrs_aa_matrix_storage(a, lda, 1, 2, n - 1, n - 1);
            a_u := array2d_view(a_u_raw, lda, row_major=false);
            b_tail_raw := _hetrs_aa_matrix_storage(b.*, ldb, 2, 1, n - 1, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _hetrs_aa_trsm(#char "L", #char "U", #char "N", #char "U", n - 1, nrhs, one_t, *a_u, lda, *b_tail, ldb, MT);

            k = n;
            while true {
                if k < 1 then break;
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_aa_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_aa_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_aa_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k -= 1;
            }
        }
    } else {
        if n > 1 {
            k = 1;
            while true {
                if k > n then break;
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_aa_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_aa_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_aa_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k += 1;
            }

            a_l_raw := _hetrs_aa_matrix_storage(a, lda, 2, 1, n - 1, n - 1);
            a_l := array2d_view(a_l_raw, lda, row_major=false);
            b_tail_raw := _hetrs_aa_matrix_storage(b.*, ldb, 2, 1, n - 1, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _hetrs_aa_trsm(#char "L", #char "L", #char "N", #char "U", n - 1, nrhs, one_t, *a_l, lda, *b_tail, ldb, MT);
        }

        for k:1..n {
            fset(work, n + k - 1, fget(a, k, k));
        }
        if n > 1 {
            for k:1..n - 1 {
                fset(work, k, fget(a, k + 1, k));
                fset(work, 2*n + k - 1, fget(a, k + 1, k));
            }
            work_du := _hetrs_aa_vector_slice(work, 2*n, n - 1);
            _hetrs_aa_lacgv(n - 1, *work_du, 1, MT);
        }
        work_dl := _hetrs_aa_vector_slice(work, 1, max(0, n - 1));
        work_d := _hetrs_aa_vector_slice(work, n, n);
        work_du := _hetrs_aa_vector_slice(work, 2*n, max(0, n - 1));
        b_full_raw := _hetrs_aa_matrix_storage(b.*, ldb, 1, 1, n, nrhs);
        b_full := array2d_view(b_full_raw, ldb, row_major=false);
        _hetrs_aa_gtsv(n, nrhs, *work_dl, *work_d, *work_du, *b_full, ldb, info, MT);

        if n > 1 {
            a_l_raw := _hetrs_aa_matrix_storage(a, lda, 2, 1, n - 1, n - 1);
            a_l := array2d_view(a_l_raw, lda, row_major=false);
            b_tail_raw := _hetrs_aa_matrix_storage(b.*, ldb, 2, 1, n - 1, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _hetrs_aa_trsm(#char "L", #char "L", #char "C", #char "U", n - 1, nrhs, one_t, *a_l, lda, *b_tail, ldb, MT);

            k = n;
            while true {
                if k < 1 then break;
                kp = fget(ipiv, k);
                if kp != k {
                    row_k := _hetrs_aa_matrix_row(b.*, ldb, k, 1, nrhs);
                    row_kp := _hetrs_aa_matrix_row(b.*, ldb, kp, 1, nrhs);
                    _hetrs_aa_swap(nrhs, *row_k, ldb, *row_kp, ldb, MT);
                }
                k -= 1;
            }
        }
    }
}
