// ? LAPACK *hetrf_aa_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrf_aa_2stage/zhetrf_aa_2stage
// *HETRF_AA_2STAGE computes the two-stage Aasen factorization of a
// complex Hermitian matrix.

chetrf_aa_2stage :: #bake_arguments _hetrf_aa_2stage(RT=float32, MT=C32);
zhetrf_aa_2stage :: #bake_arguments _hetrf_aa_2stage(RT=float64, MT=C64);

#scope_file
_hetrf_aa_2stage_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrf_aa_2stage_storage :: inline (
    data: *$T, $$ld: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] T {
    count := ifx rows <= 0 || cols <= 0 then 0 else ld*(cols - 1) + rows;
    return _hetrf_aa_2stage_raw_slice(data + (col - 1)*ld + (row - 1), count);
}

_hetrf_aa_2stage_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetrf_aa_2stage_raw_slice(x.data + (offset - 1), count);
}

_hetrf_aa_2stage_matrix_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    return _hetrf_aa_2stage_storage(a.data, lda, row, col, rows, cols);
}

_hetrf_aa_2stage_matrix_column :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    return _hetrf_aa_2stage_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetrf_aa_2stage_matrix_row :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _hetrf_aa_2stage_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_hetrf_aa_2stage_tb_view :: inline (
    tb: *$TB/ArrayType1D, $$ldtb: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] TB.data_type {
    return _hetrf_aa_2stage_storage(tb.data, ldtb, row, col, rows, cols);
}

_hetrf_aa_2stage_tb_index :: inline ($$ldtb: int, $$row: int, $$col: int) -> int {
    return row + (col - 1)*ldtb;
}

_hetrf_aa_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrf_aa_2stage_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetrf_aa_2stage_copy :: (
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

_hetrf_aa_2stage_lacgv :: (
    $$n: int, x: *$X/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && X.data_type == MT;
} {
    #if MT == C32 {
        clacgv(n, x, incx);
    } else {
        zlacgv(n, x, incx);
    }
}

_hetrf_aa_2stage_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_hetrf_aa_2stage_laset :: (
    $$uplo: u8, $$m: int, $$n: int, alpha: $MT, beta: MT, a: *$A/ArrayType2D, $$lda: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT;
} {
    #if MT == C32 {
        claset(uplo, m, n, alpha, beta, a, lda);
    } else {
        zlaset(uplo, m, n, alpha, beta, a, lda);
    }
}

_hetrf_aa_2stage_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int,
    alpha: $MT, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int,
    beta: MT, c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_hetrf_aa_2stage_getrf :: (
    $$m: int, $$n: int, a: *$A/ArrayType2D, $$lda: int, ipiv: *$IPIV/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && IPIV.data_type == int;
} {
    #if MT == C32 {
        cgetrf(m, n, a, lda, ipiv, info);
    } else {
        zgetrf(m, n, a, lda, ipiv, info);
    }
}

_hetrf_aa_2stage_hegst :: (
    $$itype: int, $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        chegst(itype, uplo, n, a, lda, b, ldb, info);
    } else {
        zhegst(itype, uplo, n, a, lda, b, ldb, info);
    }
}

_hetrf_aa_2stage_swap :: (
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

_hetrf_aa_2stage_trsm :: (
    $$side: u8, $$uplo: u8, $$transa: u8, $$diag: u8, $$m: int, $$n: int,
    alpha: $MT, a: *$A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_hetrf_aa_2stage_gbtrf :: (
    $$m: int, $$n: int, $$kl: int, $$ku: int, ab: *$AB/ArrayType2D, $$ldab: int,
    ipiv: *$IPIV/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AB.data_type == MT && IPIV.data_type == int;
} {
    #if MT == C32 {
        cgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    } else {
        zgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    }
}

_hetrf_aa_2stage :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    tb: *$TB/ArrayType1D, $$ltb: int, ipiv: *$IPIV/ArrayType1D,
    ipiv2: *$IPIV2/ArrayType1D, work: *$WORK/ArrayType1D, $$lwork: int,
    info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && TB.data_type == MT &&
           IPIV.data_type == int && IPIV2.data_type == int &&
           WORK.data_type == MT;
} {
    zero_t :: #run zero(MT);
    one_t  :: #run _hetrf_aa_2stage_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHETRF_AA_2STAGE" } else { "ZHETRF_AA_2STAGE" };

    i, j, k, i1, i2, td : int = ---;
    iinfo, jb, kb, ldtb, nt : int = ---;
    piv : MT = ---;
    tquery, upper, wquery : bool = ---;
    nb : int = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    wquery = lwork == -1;
    tquery = ltb == -1;
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if ltb < max(1, 4*n) && !tquery {
        info.* = -6;
    } else if lwork < max(1, n) && !wquery {
        info.* = -10;
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    nb = ilaenv(1, ifx MT == C32 then "CHETRF_AA_2STAGE" else "ZHETRF_AA_2STAGE", cast(*u8) *uplo, n, -1, -1, -1);
    if tquery {
        fset(tb, 1, _hetrf_aa_2stage_complex_from_real(_hetrf_aa_2stage_roundup_lwork(ncast(RT, max(1, (3*nb + 1)*n)), RT), RT, MT));
    }
    if wquery {
        fset(work, 1, _hetrf_aa_2stage_complex_from_real(_hetrf_aa_2stage_roundup_lwork(ncast(RT, max(1, n*nb)), RT), RT, MT));
    }
    if tquery || wquery {
        return;
    }

    if n == 0 {
        return;
    }

    ldtb = ltb / n;
    if ldtb < 3*nb + 1 {
        nb = (ldtb - 1) / 3;
    }
    if lwork < nb*n {
        nb = lwork / n;
    }

    nt = (n + nb - 1) / nb;
    td = 2*nb;
    kb = min(nb, n);

    for j:1..kb {
        fset(ipiv, j, j);
    }
    fset(tb, 1, _hetrf_aa_2stage_complex_from_real(ncast(RT, nb), RT, MT));

    if upper {
        for j:0..nt - 1 {
            kb = min(nb, n - j*nb);
            for i:1..j - 1 {
                if i == 1 {
                    jb = if i == j - 1 then nb + kb else 2*nb;
                    tb_block_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + 1, i*nb + 1, nb, jb);
                    tb_block := array2d_view(tb_block_raw, ldtb - 1, row_major=false);
                    a_block_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (i - 1)*nb + 1, j*nb + 1, jb, kb);
                    a_block := array2d_view(a_block_raw, lda, row_major=false);
                    w_block_raw := _hetrf_aa_2stage_storage(work.data, n, i*nb + 1, 1, nb, kb);
                    w_block := array2d_view(w_block_raw, n, row_major=false);
                    _hetrf_aa_2stage_gemm(#char "N", #char "N", nb, kb, jb, one_t, *tb_block, ldtb - 1, *a_block, lda, zero_t, *w_block, n, MT);
                } else {
                    jb = if i == j - 1 then 2*nb + kb else 3*nb;
                    tb_block_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + nb + 1, (i - 1)*nb + 1, nb, jb);
                    tb_block := array2d_view(tb_block_raw, ldtb - 1, row_major=false);
                    a_block_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (i - 2)*nb + 1, j*nb + 1, jb, kb);
                    a_block := array2d_view(a_block_raw, lda, row_major=false);
                    w_block_raw := _hetrf_aa_2stage_storage(work.data, n, i*nb + 1, 1, nb, kb);
                    w_block := array2d_view(w_block_raw, n, row_major=false);
                    _hetrf_aa_2stage_gemm(#char "N", #char "N", nb, kb, jb, one_t, *tb_block, ldtb - 1, *a_block, lda, zero_t, *w_block, n, MT);
                }
            }

            tb_diag_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + 1, j*nb + 1, kb, kb);
            tb_diag := array2d_view(tb_diag_raw, ldtb - 1, row_major=false);
            a_diag_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, j*nb + 1, kb, kb);
            a_diag := array2d_view(a_diag_raw, lda, row_major=false);
            _hetrf_aa_2stage_lacpy(#char "U", kb, kb, *a_diag, lda, *tb_diag, ldtb - 1, MT);

            if j > 1 {
                a1_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, 1, j*nb + 1, (j - 1)*nb, kb);
                a1 := array2d_view(a1_raw, lda, row_major=false);
                w1_raw := _hetrf_aa_2stage_storage(work.data, n, nb + 1, 1, (j - 1)*nb, kb);
                w1 := array2d_view(w1_raw, n, row_major=false);
                _hetrf_aa_2stage_gemm(#char "C", #char "N", kb, kb, (j - 1)*nb, -one_t, *a1, lda, *w1, n, one_t, *tb_diag, ldtb - 1, MT);

                a2_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j - 1)*nb + 1, j*nb + 1, kb, kb);
                a2 := array2d_view(a2_raw, lda, row_major=false);
                tb_off_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + nb + 1, (j - 1)*nb + 1, kb, nb);
                tb_off := array2d_view(tb_off_raw, ldtb - 1, row_major=false);
                work_top_raw := _hetrf_aa_2stage_storage(work.data, n, 1, 1, kb, nb);
                work_top := array2d_view(work_top_raw, n, row_major=false);
                _hetrf_aa_2stage_gemm(#char "C", #char "N", kb, nb, kb, one_t, *a2, lda, *tb_off, ldtb - 1, zero_t, *work_top, n, MT);

                a3_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j - 2)*nb + 1, j*nb + 1, nb, kb);
                a3 := array2d_view(a3_raw, lda, row_major=false);
                _hetrf_aa_2stage_gemm(#char "N", #char "N", kb, kb, nb, -one_t, *work_top, n, *a3, lda, one_t, *tb_diag, ldtb - 1, MT);
            }

            if j > 0 {
                a_prev_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j - 1)*nb + 1, j*nb + 1, kb, kb);
                a_prev := array2d_view(a_prev_raw, lda, row_major=false);
                _hetrf_aa_2stage_hegst(1, #char "U", kb, *tb_diag, ldtb - 1, *a_prev, lda, *info, MT);
                info.* = 0;
            }

            for i:1..kb {
                diag_index := _hetrf_aa_2stage_tb_index(ldtb, td + 1, j*nb + i);
                fset(tb, diag_index, _hetrf_aa_2stage_complex_from_real(real(fget(tb, diag_index)), RT, MT));
                for k:i + 1..kb {
                    dst_index := _hetrf_aa_2stage_tb_index(ldtb, td + (k - i) + 1, j*nb + i);
                    src_index := _hetrf_aa_2stage_tb_index(ldtb, td - (k - (i + 1)), j*nb + k);
                    fset(tb, dst_index, conjugate(fget(tb, src_index)));
                }
            }

            if j < nt - 1 {
                if j > 0 {
                    if j == 1 {
                        tb_h_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + 1, j*nb + 1, kb, kb);
                        tb_h := array2d_view(tb_h_raw, ldtb - 1, row_major=false);
                        a_h_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j - 1)*nb + 1, j*nb + 1, kb, kb);
                        a_h := array2d_view(a_h_raw, lda, row_major=false);
                        w_h_raw := _hetrf_aa_2stage_storage(work.data, n, j*nb + 1, 1, kb, kb);
                        w_h := array2d_view(w_h_raw, n, row_major=false);
                        _hetrf_aa_2stage_gemm(#char "N", #char "N", kb, kb, kb, one_t, *tb_h, ldtb - 1, *a_h, lda, zero_t, *w_h, n, MT);
                    } else {
                        tb_h_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + nb + 1, (j - 1)*nb + 1, kb, nb + kb);
                        tb_h := array2d_view(tb_h_raw, ldtb - 1, row_major=false);
                        a_h_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j - 2)*nb + 1, j*nb + 1, nb + kb, kb);
                        a_h := array2d_view(a_h_raw, lda, row_major=false);
                        w_h_raw := _hetrf_aa_2stage_storage(work.data, n, j*nb + 1, 1, kb, kb);
                        w_h := array2d_view(w_h_raw, n, row_major=false);
                        _hetrf_aa_2stage_gemm(#char "N", #char "N", kb, kb, nb + kb, one_t, *tb_h, ldtb - 1, *a_h, lda, zero_t, *w_h, n, MT);
                    }

                    w_prev_raw := _hetrf_aa_2stage_storage(work.data, n, nb + 1, 1, j*nb, n - (j + 1)*nb);
                    w_prev := array2d_view(w_prev_raw, n, row_major=false);
                    a_trail_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, 1, (j + 1)*nb + 1, j*nb, n - (j + 1)*nb);
                    a_trail := array2d_view(a_trail_raw, lda, row_major=false);
                    a_upd_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (j + 1)*nb + 1, nb, n - (j + 1)*nb);
                    a_upd := array2d_view(a_upd_raw, lda, row_major=false);
                    _hetrf_aa_2stage_gemm(#char "C", #char "N", nb, n - (j + 1)*nb, j*nb, -one_t, *w_prev, n, *a_trail, lda, one_t, *a_upd, lda, MT);
                }

                for k:1..nb {
                    count := n - (j + 1)*nb;
                    a_row := _hetrf_aa_2stage_matrix_row(a.*, lda, j*nb + k, (j + 1)*nb + 1, count);
                    work_col := _hetrf_aa_2stage_vector_slice(work, 1 + (k - 1)*n, count);
                    _hetrf_aa_2stage_copy(count, *a_row, lda, *work_col, 1, MT);
                }

                work_panel_raw := _hetrf_aa_2stage_storage(work.data, n, 1, 1, n - (j + 1)*nb, nb);
                work_panel := array2d_view(work_panel_raw, n, row_major=false);
                ipiv_panel := _hetrf_aa_2stage_vector_slice(ipiv, (j + 1)*nb + 1, n - (j + 1)*nb);
                _hetrf_aa_2stage_getrf(n - (j + 1)*nb, nb, *work_panel, n, *ipiv_panel, *info, MT);
                info.* = 0;

                for k:1..nb {
                    count := n - k - (j + 1)*nb;
                    work_col := _hetrf_aa_2stage_vector_slice(work, k + 1 + (k - 1)*n, count);
                    a_col := _hetrf_aa_2stage_matrix_column(a.*, lda, j*nb + k, (j + 1)*nb + k + 1, count);
                    _hetrf_aa_2stage_copy(count, *work_col, 1, *a_col, lda, MT);

                    work_head := _hetrf_aa_2stage_vector_slice(work, 1 + (k - 1)*n, k);
                    _hetrf_aa_2stage_lacgv(k, *work_head, 1, MT);
                }

                kb = min(nb, n - (j + 1)*nb);
                tb_sub_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + nb + 1, j*nb + 1, kb, nb);
                tb_sub := array2d_view(tb_sub_raw, ldtb - 1, row_major=false);
                _hetrf_aa_2stage_laset(#char "F", kb, nb, zero_t, zero_t, *tb_sub, ldtb - 1, MT);

                work_upper_raw := _hetrf_aa_2stage_storage(work.data, n, 1, 1, kb, nb);
                work_upper := array2d_view(work_upper_raw, n, row_major=false);
                _hetrf_aa_2stage_lacpy(#char "U", kb, nb, *work_upper, n, *tb_sub, ldtb - 1, MT);

                if j > 0 {
                    a_tri_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j - 1)*nb + 1, j*nb + 1, nb, kb);
                    a_tri := array2d_view(a_tri_raw, lda, row_major=false);
                    _hetrf_aa_2stage_trsm(#char "R", #char "U", #char "N", #char "U", kb, nb, one_t, *a_tri, lda, *tb_sub, ldtb - 1, MT);
                }

                for k:1..nb {
                    for i:1..kb {
                        dst_index := _hetrf_aa_2stage_tb_index(ldtb, td - nb + k - i + 1, j*nb + nb + i);
                        src_index := _hetrf_aa_2stage_tb_index(ldtb, td + nb + i - k + 1, j*nb + k);
                        fset(tb, dst_index, conjugate(fget(tb, src_index)));
                    }
                }

                a_unit_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (j + 1)*nb + 1, kb, nb);
                a_unit := array2d_view(a_unit_raw, lda, row_major=false);
                _hetrf_aa_2stage_laset(#char "L", kb, nb, zero_t, one_t, *a_unit, lda, MT);

                for k:1..kb {
                    fset(ipiv, (j + 1)*nb + k, fget(ipiv.*, (j + 1)*nb + k) + (j + 1)*nb);
                    i1 = (j + 1)*nb + k;
                    i2 = fget(ipiv.*, (j + 1)*nb + k);
                    if i1 != i2 {
                        prev_col1 := _hetrf_aa_2stage_matrix_column(a.*, lda, (j + 1)*nb + 1, i1, k - 1);
                        prev_col2 := _hetrf_aa_2stage_matrix_column(a.*, lda, (j + 1)*nb + 1, i2, k - 1);
                        _hetrf_aa_2stage_swap(k - 1, *prev_col1, 1, *prev_col2, 1, MT);

                        if i2 > i1 + 1 {
                            row_seg := _hetrf_aa_2stage_matrix_row(a.*, lda, i1, i1 + 1, i2 - i1 - 1);
                            col_seg := _hetrf_aa_2stage_matrix_column(a.*, lda, i1 + 1, i2, i2 - i1 - 1);
                            _hetrf_aa_2stage_swap(i2 - i1 - 1, *row_seg, lda, *col_seg, 1, MT);
                            _hetrf_aa_2stage_lacgv(i2 - i1 - 1, *col_seg, 1, MT);
                        }

                        row_diag := _hetrf_aa_2stage_matrix_row(a.*, lda, i1, i1 + 1, i2 - i1);
                        _hetrf_aa_2stage_lacgv(i2 - i1, *row_diag, lda, MT);

                        if i2 < n {
                            row_tail1 := _hetrf_aa_2stage_matrix_row(a.*, lda, i1, i2 + 1, n - i2);
                            row_tail2 := _hetrf_aa_2stage_matrix_row(a.*, lda, i2, i2 + 1, n - i2);
                            _hetrf_aa_2stage_swap(n - i2, *row_tail1, lda, *row_tail2, lda, MT);
                        }

                        piv = fget(a.*, i1, i1);
                        fset(a, i1, i1, fget(a.*, i2, i2));
                        fset(a, i2, i2, piv);

                        if j > 0 {
                            prev_top1 := _hetrf_aa_2stage_matrix_column(a.*, lda, 1, i1, j*nb);
                            prev_top2 := _hetrf_aa_2stage_matrix_column(a.*, lda, 1, i2, j*nb);
                            _hetrf_aa_2stage_swap(j*nb, *prev_top1, 1, *prev_top2, 1, MT);
                        }
                    }
                }
            }
        }
    } else {
        for j:0..nt - 1 {
            kb = min(nb, n - j*nb);
            for i:1..j - 1 {
                if i == 1 {
                    jb = if i == j - 1 then nb + kb else 2*nb;
                    tb_block_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + 1, i*nb + 1, nb, jb);
                    tb_block := array2d_view(tb_block_raw, ldtb - 1, row_major=false);
                    a_block_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (i - 1)*nb + 1, kb, jb);
                    a_block := array2d_view(a_block_raw, lda, row_major=false);
                    w_block_raw := _hetrf_aa_2stage_storage(work.data, n, i*nb + 1, 1, nb, kb);
                    w_block := array2d_view(w_block_raw, n, row_major=false);
                    _hetrf_aa_2stage_gemm(#char "N", #char "C", nb, kb, jb, one_t, *tb_block, ldtb - 1, *a_block, lda, zero_t, *w_block, n, MT);
                } else {
                    jb = if i == j - 1 then 2*nb + kb else 3*nb;
                    tb_block_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + nb + 1, (i - 1)*nb + 1, nb, jb);
                    tb_block := array2d_view(tb_block_raw, ldtb - 1, row_major=false);
                    a_block_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (i - 2)*nb + 1, kb, jb);
                    a_block := array2d_view(a_block_raw, lda, row_major=false);
                    w_block_raw := _hetrf_aa_2stage_storage(work.data, n, i*nb + 1, 1, nb, kb);
                    w_block := array2d_view(w_block_raw, n, row_major=false);
                    _hetrf_aa_2stage_gemm(#char "N", #char "C", nb, kb, jb, one_t, *tb_block, ldtb - 1, *a_block, lda, zero_t, *w_block, n, MT);
                }
            }

            tb_diag_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + 1, j*nb + 1, kb, kb);
            tb_diag := array2d_view(tb_diag_raw, ldtb - 1, row_major=false);
            a_diag_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, j*nb + 1, kb, kb);
            a_diag := array2d_view(a_diag_raw, lda, row_major=false);
            _hetrf_aa_2stage_lacpy(#char "L", kb, kb, *a_diag, lda, *tb_diag, ldtb - 1, MT);

            if j > 1 {
                a1_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, 1, kb, (j - 1)*nb);
                a1 := array2d_view(a1_raw, lda, row_major=false);
                w1_raw := _hetrf_aa_2stage_storage(work.data, n, nb + 1, 1, (j - 1)*nb, kb);
                w1 := array2d_view(w1_raw, n, row_major=false);
                _hetrf_aa_2stage_gemm(#char "N", #char "N", kb, kb, (j - 1)*nb, -one_t, *a1, lda, *w1, n, one_t, *tb_diag, ldtb - 1, MT);

                a2_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (j - 1)*nb + 1, kb, kb);
                a2 := array2d_view(a2_raw, lda, row_major=false);
                tb_off_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + nb + 1, (j - 1)*nb + 1, kb, nb);
                tb_off := array2d_view(tb_off_raw, ldtb - 1, row_major=false);
                work_top_raw := _hetrf_aa_2stage_storage(work.data, n, 1, 1, kb, nb);
                work_top := array2d_view(work_top_raw, n, row_major=false);
                _hetrf_aa_2stage_gemm(#char "N", #char "N", kb, nb, kb, one_t, *a2, lda, *tb_off, ldtb - 1, zero_t, *work_top, n, MT);

                a3_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (j - 2)*nb + 1, kb, nb);
                a3 := array2d_view(a3_raw, lda, row_major=false);
                _hetrf_aa_2stage_gemm(#char "N", #char "C", kb, kb, nb, -one_t, *work_top, n, *a3, lda, one_t, *tb_diag, ldtb - 1, MT);
            }

            if j > 0 {
                a_prev_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (j - 1)*nb + 1, kb, kb);
                a_prev := array2d_view(a_prev_raw, lda, row_major=false);
                _hetrf_aa_2stage_hegst(1, #char "L", kb, *tb_diag, ldtb - 1, *a_prev, lda, *info, MT);
                info.* = 0;
            }

            for i:1..kb {
                diag_index := _hetrf_aa_2stage_tb_index(ldtb, td + 1, j*nb + i);
                fset(tb, diag_index, _hetrf_aa_2stage_complex_from_real(real(fget(tb, diag_index)), RT, MT));
                for k:i + 1..kb {
                    dst_index := _hetrf_aa_2stage_tb_index(ldtb, td - (k - (i + 1)), j*nb + k);
                    src_index := _hetrf_aa_2stage_tb_index(ldtb, td + (k - i) + 1, j*nb + i);
                    fset(tb, dst_index, conjugate(fget(tb, src_index)));
                }
            }

            if j < nt - 1 {
                if j > 0 {
                    if j == 1 {
                        tb_h_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + 1, j*nb + 1, kb, kb);
                        tb_h := array2d_view(tb_h_raw, ldtb - 1, row_major=false);
                        a_h_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (j - 1)*nb + 1, kb, kb);
                        a_h := array2d_view(a_h_raw, lda, row_major=false);
                        w_h_raw := _hetrf_aa_2stage_storage(work.data, n, j*nb + 1, 1, kb, kb);
                        w_h := array2d_view(w_h_raw, n, row_major=false);
                        _hetrf_aa_2stage_gemm(#char "N", #char "C", kb, kb, kb, one_t, *tb_h, ldtb - 1, *a_h, lda, zero_t, *w_h, n, MT);
                    } else {
                        tb_h_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + nb + 1, (j - 1)*nb + 1, kb, nb + kb);
                        tb_h := array2d_view(tb_h_raw, ldtb - 1, row_major=false);
                        a_h_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (j - 2)*nb + 1, kb, nb + kb);
                        a_h := array2d_view(a_h_raw, lda, row_major=false);
                        w_h_raw := _hetrf_aa_2stage_storage(work.data, n, j*nb + 1, 1, kb, kb);
                        w_h := array2d_view(w_h_raw, n, row_major=false);
                        _hetrf_aa_2stage_gemm(#char "N", #char "C", kb, kb, nb + kb, one_t, *tb_h, ldtb - 1, *a_h, lda, zero_t, *w_h, n, MT);
                    }

                    a_prev_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j + 1)*nb + 1, 1, n - (j + 1)*nb, j*nb);
                    a_prev := array2d_view(a_prev_raw, lda, row_major=false);
                    w_prev_raw := _hetrf_aa_2stage_storage(work.data, n, nb + 1, 1, j*nb, nb);
                    w_prev := array2d_view(w_prev_raw, n, row_major=false);
                    a_upd_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j + 1)*nb + 1, j*nb + 1, n - (j + 1)*nb, nb);
                    a_upd := array2d_view(a_upd_raw, lda, row_major=false);
                    _hetrf_aa_2stage_gemm(#char "N", #char "N", n - (j + 1)*nb, nb, j*nb, -one_t, *a_prev, lda, *w_prev, n, one_t, *a_upd, lda, MT);
                }

                panel_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j + 1)*nb + 1, j*nb + 1, n - (j + 1)*nb, nb);
                panel := array2d_view(panel_raw, lda, row_major=false);
                ipiv_panel := _hetrf_aa_2stage_vector_slice(ipiv, (j + 1)*nb + 1, n - (j + 1)*nb);
                _hetrf_aa_2stage_getrf(n - (j + 1)*nb, nb, *panel, lda, *ipiv_panel, *info, MT);
                info.* = 0;

                kb = min(nb, n - (j + 1)*nb);
                tb_sub_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, td + nb + 1, j*nb + 1, kb, nb);
                tb_sub := array2d_view(tb_sub_raw, ldtb - 1, row_major=false);
                _hetrf_aa_2stage_laset(#char "F", kb, nb, zero_t, zero_t, *tb_sub, ldtb - 1, MT);
                panel_upper_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j + 1)*nb + 1, j*nb + 1, kb, nb);
                panel_upper := array2d_view(panel_upper_raw, lda, row_major=false);
                _hetrf_aa_2stage_lacpy(#char "U", kb, nb, *panel_upper, lda, *tb_sub, ldtb - 1, MT);

                if j > 0 {
                    a_tri_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, j*nb + 1, (j - 1)*nb + 1, kb, nb);
                    a_tri := array2d_view(a_tri_raw, lda, row_major=false);
                    _hetrf_aa_2stage_trsm(#char "R", #char "L", #char "C", #char "U", kb, nb, one_t, *a_tri, lda, *tb_sub, ldtb - 1, MT);
                }

                for k:1..nb {
                    for i:1..kb {
                        dst_index := _hetrf_aa_2stage_tb_index(ldtb, td - nb + k - i + 1, j*nb + nb + i);
                        src_index := _hetrf_aa_2stage_tb_index(ldtb, td + nb + i - k + 1, j*nb + k);
                        fset(tb, dst_index, conjugate(fget(tb, src_index)));
                    }
                }

                a_unit_raw := _hetrf_aa_2stage_matrix_view(a.*, lda, (j + 1)*nb + 1, j*nb + 1, kb, nb);
                a_unit := array2d_view(a_unit_raw, lda, row_major=false);
                _hetrf_aa_2stage_laset(#char "U", kb, nb, zero_t, one_t, *a_unit, lda, MT);

                for k:1..kb {
                    fset(ipiv, (j + 1)*nb + k, fget(ipiv.*, (j + 1)*nb + k) + (j + 1)*nb);
                    i1 = (j + 1)*nb + k;
                    i2 = fget(ipiv.*, (j + 1)*nb + k);
                    if i1 != i2 {
                        prev_row1 := _hetrf_aa_2stage_matrix_row(a.*, lda, i1, (j + 1)*nb + 1, k - 1);
                        prev_row2 := _hetrf_aa_2stage_matrix_row(a.*, lda, i2, (j + 1)*nb + 1, k - 1);
                        _hetrf_aa_2stage_swap(k - 1, *prev_row1, lda, *prev_row2, lda, MT);

                        if i2 > i1 + 1 {
                            col_seg := _hetrf_aa_2stage_matrix_column(a.*, lda, i1 + 1, i1, i2 - i1 - 1);
                            row_seg := _hetrf_aa_2stage_matrix_row(a.*, lda, i2, i1 + 1, i2 - i1 - 1);
                            _hetrf_aa_2stage_swap(i2 - i1 - 1, *col_seg, 1, *row_seg, lda, MT);
                            _hetrf_aa_2stage_lacgv(i2 - i1 - 1, *row_seg, lda, MT);
                        }

                        col_diag := _hetrf_aa_2stage_matrix_column(a.*, lda, i1 + 1, i1, i2 - i1);
                        _hetrf_aa_2stage_lacgv(i2 - i1, *col_diag, 1, MT);

                        if i2 < n {
                            col_tail1 := _hetrf_aa_2stage_matrix_column(a.*, lda, i2 + 1, i1, n - i2);
                            col_tail2 := _hetrf_aa_2stage_matrix_column(a.*, lda, i2 + 1, i2, n - i2);
                            _hetrf_aa_2stage_swap(n - i2, *col_tail1, 1, *col_tail2, 1, MT);
                        }

                        piv = fget(a.*, i1, i1);
                        fset(a, i1, i1, fget(a.*, i2, i2));
                        fset(a, i2, i2, piv);

                        if j > 0 {
                            prev_left1 := _hetrf_aa_2stage_matrix_row(a.*, lda, i1, 1, j*nb);
                            prev_left2 := _hetrf_aa_2stage_matrix_row(a.*, lda, i2, 1, j*nb);
                            _hetrf_aa_2stage_swap(j*nb, *prev_left1, lda, *prev_left2, lda, MT);
                        }
                    }
                }
            }
        }
    }

    tb_band_raw := _hetrf_aa_2stage_tb_view(tb, ldtb, 1, 1, ldtb, n);
    tb_band := array2d_view(tb_band_raw, ldtb, row_major=false);
    _hetrf_aa_2stage_gbtrf(n, n, nb, nb, *tb_band, ldtb, ipiv2, info, MT);
}
