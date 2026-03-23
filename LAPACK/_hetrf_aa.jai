// ? LAPACK *hetrf_aa()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrf_aa/zhetrf_aa
// *HETRF_AA computes the blocked Aasen factorization of a complex
// Hermitian matrix.

chetrf_aa :: #bake_arguments _hetrf_aa(RT=float32, MT=C32);
zhetrf_aa :: #bake_arguments _hetrf_aa(RT=float64, MT=C64);

#scope_file
_hetrf_aa_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrf_aa_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetrf_aa_raw_slice(x.data + (offset - 1), count);
}

_hetrf_aa_matrix_storage :: inline (data: *$T, $$ld: int, $$rows: int, $$cols: int) -> [] T {
    return _hetrf_aa_raw_slice(data, ifx rows <= 0 || cols <= 0 then 0 else ld*(cols - 1) + rows);
}

_hetrf_aa_matrix_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    return _hetrf_aa_matrix_storage(a.data + (col - 1)*lda + (row - 1), lda, rows, cols);
}

_hetrf_aa_matrix_column :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    return _hetrf_aa_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetrf_aa_matrix_row :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _hetrf_aa_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_hetrf_aa_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrf_aa_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetrf_aa_lahef :: (
    $$uplo: u8, $$j1: int, $$m: int, $$nb: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, h: *$H/ArrayType2D, $$ldh: int, work: *$WORK/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && IPIV.data_type == int &&
           H.data_type == MT && WORK.data_type == MT;
} {
    #if MT == C32 {
        clahef_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
    } else {
        zlahef_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
    }
}

_hetrf_aa_copy :: (
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

_hetrf_aa_swap :: (
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

_hetrf_aa_scal :: (
    $$n: int, alpha: $MT, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        cscal(n, alpha, x, incx);
    } else {
        zscal(n, alpha, x, incx);
    }
}

_hetrf_aa_gemm :: (
    $$transa: u8, $$transb: u8,
    $$m: int, $$n: int, $$k: int,
    alpha: $MT, a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int,
    beta: MT, c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_hetrf_aa :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    ipiv: *$IPIV/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && IPIV.data_type == int && WORK.data_type == MT;
} {
    zero_t :: #run zero(MT);
    one_t  :: #run _hetrf_aa_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHETRF_AA" } else { "ZHETRF_AA" };

    lquery, upper : bool = ---;
    j, jb, j1, j2, j3, k1, k2, lwjmin, lwjopt, mj, nb, nj : int = ---;
    alpha : MT = ---;

    nb = ilaenv(1, ifx MT == C32 then "CHETRF_AA" else "ZHETRF_AA", cast(*u8) *uplo, n, -1, -1, -1);

    info.* = 0;
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;
    if n <= 1 {
        lwjmin = 1;
        lwjopt = 1;
    } else {
        lwjmin = 2*n;
        lwjopt = (nb + 1)*n;
    }

    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if lwork < lwjmin && !lquery {
        info.* = -7;
    }

    if info.* == 0 {
        fset(work, 1, _hetrf_aa_complex_from_real(_hetrf_aa_roundup_lwork(ncast(RT, lwjopt), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    if n == 0 {
        return;
    }
    fset(ipiv, 1, 1);
    if n == 1 {
        fset(a, 1, 1, _hetrf_aa_complex_from_real(real(fget(a.*, 1, 1)), RT, MT));
        return;
    }

    if lwork < (1 + nb)*n {
        nb = (lwork - n) / n;
    }

    work_h_raw := _hetrf_aa_matrix_storage(work.data, n, n, max(1, nb));
    work_h := array2d_view(work_h_raw, n, row_major=false);
    work_aux := _hetrf_aa_vector_slice(work, n*nb + 1, max(0, lwork - n*nb));

    if upper {
        first_row := _hetrf_aa_matrix_row(a.*, lda, 1, 1, n);
        work_col1 := _hetrf_aa_vector_slice(work, 1, n);
        _hetrf_aa_copy(n, *first_row, lda, *work_col1, 1, MT);

        j = 0;
        while true {
            if j >= n then break;

            j1 = j + 1;
            jb = min(n - j1 + 1, nb);
            k1 = max(1, j) - j;

            a_panel_raw := _hetrf_aa_matrix_view(a.*, lda, max(1, j), j + 1, n - max(1, j) + 1, n - j);
            a_panel := array2d_view(a_panel_raw, lda, row_major=false);
            ipiv_panel := _hetrf_aa_vector_slice(ipiv, j + 1, n - j);
            _hetrf_aa_lahef(uplo, 2 - k1, n - j, jb, *a_panel, lda, *ipiv_panel, *work_h, n, *work_aux, MT);

            for j2:j + 2..min(n, j + jb + 1) {
                fset(ipiv, j2, fget(ipiv.*, j2) + j);
                if j2 != fget(ipiv.*, j2) && (j1 - k1) > 2 {
                    col_j2 := _hetrf_aa_matrix_column(a.*, lda, 1, j2, j1 - k1 - 2);
                    col_pv := _hetrf_aa_matrix_column(a.*, lda, 1, fget(ipiv.*, j2), j1 - k1 - 2);
                    _hetrf_aa_swap(j1 - k1 - 2, *col_j2, 1, *col_pv, 1, MT);
                }
            }
            j += jb;

            if j < n {
                if j1 > 1 || jb > 1 {
                    alpha = conjugate(fget(a.*, j, j + 1));
                    fset(a, j, j + 1, one_t);

                    h_tail := _hetrf_aa_vector_slice(work, (j + 1 - j1 + 1) + jb*n, n - j);
                    row_prev := _hetrf_aa_matrix_row(a.*, lda, j - 1, j + 1, n - j);
                    _hetrf_aa_copy(n - j, *row_prev, lda, *h_tail, 1, MT);
                    _hetrf_aa_scal(n - j, alpha, *h_tail, 1, MT);

                    if j1 > 1 {
                        k2 = 1;
                    } else {
                        k2 = 0;
                        jb -= 1;
                    }

                    for j2:j + 1..n;nb {
                        nj = min(nb, n - j2 + 1);
                        j3 = j2;
                        for <mj:1..nj - 1 {
                            a_col_raw := _hetrf_aa_matrix_view(a.*, lda, j1 - k2, j3, jb + 1, 1);
                            a_col := array2d_view(a_col_raw, lda, row_major=false);
                            h_blk_raw := _hetrf_aa_matrix_storage(work.data + ((j3 - j1 + 1) + k1*n - 1), n, mj, jb + 1);
                            h_blk := array2d_view(h_blk_raw, n, row_major=false);
                            c_blk_raw := _hetrf_aa_matrix_view(a.*, lda, j3, j3, 1, mj);
                            c_blk := array2d_view(c_blk_raw, lda, row_major=false);
                            _hetrf_aa_gemm(#char "C", #char "T", 1, mj, jb + 1, -one_t, *a_col, lda, *h_blk, n, one_t, *c_blk, lda, MT);
                            j3 += 1;
                        }

                        if j3 <= n {
                            a_blk_raw := _hetrf_aa_matrix_view(a.*, lda, j1 - k2, j2, jb + 1, nj);
                            a_blk := array2d_view(a_blk_raw, lda, row_major=false);
                            h_blk_raw := _hetrf_aa_matrix_storage(work.data + ((j3 - j1 + 1) + k1*n - 1), n, n - j3 + 1, jb + 1);
                            h_blk := array2d_view(h_blk_raw, n, row_major=false);
                            c_blk_raw := _hetrf_aa_matrix_view(a.*, lda, j2, j3, nj, n - j3 + 1);
                            c_blk := array2d_view(c_blk_raw, lda, row_major=false);
                            _hetrf_aa_gemm(#char "C", #char "T", nj, n - j3 + 1, jb + 1, -one_t, *a_blk, lda, *h_blk, n, one_t, *c_blk, lda, MT);
                        }
                    }

                    fset(a, j, j + 1, conjugate(alpha));
                }

                work_col1 = _hetrf_aa_vector_slice(work, 1, n - j);
                diag_row := _hetrf_aa_matrix_row(a.*, lda, j + 1, j + 1, n - j);
                _hetrf_aa_copy(n - j, *diag_row, lda, *work_col1, 1, MT);
            }
        }
    } else {
        first_col := _hetrf_aa_matrix_column(a.*, lda, 1, 1, n);
        work_col1 := _hetrf_aa_vector_slice(work, 1, n);
        _hetrf_aa_copy(n, *first_col, 1, *work_col1, 1, MT);

        j = 0;
        while true {
            if j >= n then break;

            j1 = j + 1;
            jb = min(n - j1 + 1, nb);
            k1 = max(1, j) - j;

            a_panel_raw := _hetrf_aa_matrix_view(a.*, lda, j + 1, max(1, j), n - j, n - max(1, j) + 1);
            a_panel := array2d_view(a_panel_raw, lda, row_major=false);
            ipiv_panel := _hetrf_aa_vector_slice(ipiv, j + 1, n - j);
            _hetrf_aa_lahef(uplo, 2 - k1, n - j, jb, *a_panel, lda, *ipiv_panel, *work_h, n, *work_aux, MT);

            for j2:j + 2..min(n, j + jb + 1) {
                fset(ipiv, j2, fget(ipiv.*, j2) + j);
                if j2 != fget(ipiv.*, j2) && (j1 - k1) > 2 {
                    row_j2 := _hetrf_aa_matrix_row(a.*, lda, j2, 1, j1 - k1 - 2);
                    row_pv := _hetrf_aa_matrix_row(a.*, lda, fget(ipiv.*, j2), 1, j1 - k1 - 2);
                    _hetrf_aa_swap(j1 - k1 - 2, *row_j2, lda, *row_pv, lda, MT);
                }
            }
            j += jb;

            if j < n {
                if j1 > 1 || jb > 1 {
                    alpha = conjugate(fget(a.*, j + 1, j));
                    fset(a, j + 1, j, one_t);

                    h_tail := _hetrf_aa_vector_slice(work, (j + 1 - j1 + 1) + jb*n, n - j);
                    col_prev := _hetrf_aa_matrix_column(a.*, lda, j + 1, j - 1, n - j);
                    _hetrf_aa_copy(n - j, *col_prev, 1, *h_tail, 1, MT);
                    _hetrf_aa_scal(n - j, alpha, *h_tail, 1, MT);

                    if j1 > 1 {
                        k2 = 1;
                    } else {
                        k2 = 0;
                        jb -= 1;
                    }

                    for j2:j + 1..n;nb {
                        nj = min(nb, n - j2 + 1);
                        j3 = j2;
                        for <mj:1..nj - 1 {
                            h_blk_raw := _hetrf_aa_matrix_storage(work.data + ((j3 - j1 + 1) + k1*n - 1), n, mj, jb + 1);
                            h_blk := array2d_view(h_blk_raw, n, row_major=false);
                            a_row_raw := _hetrf_aa_matrix_view(a.*, lda, j3, j1 - k2, 1, jb + 1);
                            a_row := array2d_view(a_row_raw, lda, row_major=false);
                            c_blk_raw := _hetrf_aa_matrix_view(a.*, lda, j3, j3, mj, 1);
                            c_blk := array2d_view(c_blk_raw, lda, row_major=false);
                            _hetrf_aa_gemm(#char "N", #char "C", mj, 1, jb + 1, -one_t, *h_blk, n, *a_row, lda, one_t, *c_blk, lda, MT);
                            j3 += 1;
                        }

                        if j3 <= n {
                            h_blk_raw := _hetrf_aa_matrix_storage(work.data + ((j3 - j1 + 1) + k1*n - 1), n, n - j3 + 1, jb + 1);
                            h_blk := array2d_view(h_blk_raw, n, row_major=false);
                            a_blk_raw := _hetrf_aa_matrix_view(a.*, lda, j2, j1 - k2, nj, jb + 1);
                            a_blk := array2d_view(a_blk_raw, lda, row_major=false);
                            c_blk_raw := _hetrf_aa_matrix_view(a.*, lda, j3, j2, n - j3 + 1, nj);
                            c_blk := array2d_view(c_blk_raw, lda, row_major=false);
                            _hetrf_aa_gemm(#char "N", #char "C", n - j3 + 1, nj, jb + 1, -one_t, *h_blk, n, *a_blk, lda, one_t, *c_blk, lda, MT);
                        }
                    }

                    fset(a, j + 1, j, conjugate(alpha));
                }

                work_col1 = _hetrf_aa_vector_slice(work, 1, n - j);
                diag_col := _hetrf_aa_matrix_column(a.*, lda, j + 1, j + 1, n - j);
                _hetrf_aa_copy(n - j, *diag_col, 1, *work_col1, 1, MT);
            }
        }
    }

    fset(work, 1, _hetrf_aa_complex_from_real(_hetrf_aa_roundup_lwork(ncast(RT, lwjopt), RT), RT, MT));
}
