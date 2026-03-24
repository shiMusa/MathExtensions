// ? LAPACK *labrd()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slabrd/dlabrd/clabrd/zlabrd
// *LABRD reduces the first NB rows and columns of a general matrix
// to bidiagonal form.

slabrd :: #bake_arguments _labrd(RT=float32, MT=float32);
dlabrd :: #bake_arguments _labrd(RT=float64, MT=float64);
clabrd :: #bake_arguments _labrd(RT=float32, MT=C32);
zlabrd :: #bake_arguments _labrd(RT=float64, MT=C64);

#scope_file
_labrd_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_labrd_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _labrd_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_labrd_matrix_row :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _labrd_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_labrd_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _labrd_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_labrd_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_labrd_vector_entry_ptr :: inline (x: *$X/ArrayType1D, $$index: int) -> *X.data_type {
    return x.data + (index - 1);
}

_labrd_real_part :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return x;
    } else {
        return real(x);
    }
}

_labrd_store_real :: inline (
    a: *$A/ArrayType2D,
    $$row: int,
    $$col: int,
    value: $RT,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT;
} {
    fset(a, row, col, ncast(MT, value));
}

_labrd_gemv :: (
    $$trans: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    x: $X/ArrayType1D,
    $$incx: int,
    beta: MT,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == float64 {
        dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else if MT == C32 {
        cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    } else {
        zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
}

_labrd_lacgv :: (
    $$n: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == C32 {
        clacgv(n, x, incx);
    } else if MT == C64 {
        zlacgv(n, x, incx);
    }
}

_labrd_larfg :: (
    $$n: int,
    alpha: *$MT,
    x: *$X/ArrayType1D,
    $$incx: int,
    tau: *$TAU/ArrayType1D,
    $$tau_index: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           TAU.data_type == MT;
} {
    tau_ptr := _labrd_vector_entry_ptr(tau, tau_index);
    #if MT == float32 {
        slarfg(n, alpha, x, incx, tau_ptr);
    } else if MT == float64 {
        dlarfg(n, alpha, x, incx, tau_ptr);
    } else if MT == C32 {
        clarfg(n, alpha, x, incx, tau_ptr);
    } else {
        zlarfg(n, alpha, x, incx, tau_ptr);
    }
}

_labrd_scal :: (
    $$n: int,
    alpha: $MT,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        sscal(n, alpha, x, incx);
    } else if MT == float64 {
        dscal(n, alpha, x, incx);
    } else if MT == C32 {
        cscal(n, alpha, x, incx);
    } else {
        zscal(n, alpha, x, incx);
    }
}

_labrd :: (
    $$m: int,
    $$n: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    tauq: *$TAUQ/ArrayType1D,
    taup: *$TAUP/ArrayType1D,
    x: *$X/ArrayType2D,
    $$ldx: int,
    y: *$Y/ArrayType2D,
    $$ldy: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           A.data_type == MT &&
           D.data_type == RT &&
           E.data_type == RT &&
           TAUQ.data_type == MT &&
           TAUP.data_type == MT &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    zero_t :: #run zero(MT);
    one_t  :: #run one(MT);
    minus_one_t :: #run -one(MT);
    trans_conj :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };

    alpha : MT = ---;
    i : int = ---;

    if m <= 0 || n <= 0 then return;

    if m >= n {
        for i:1..nb {
            a_i1_row := _labrd_matrix_row(a.*, lda, i, 1, i - 1);
            y_i1_row := _labrd_matrix_row(y.*, ldy, i, 1, i - 1);
            a_i_col_tail_raw := _labrd_matrix_column(a.*, lda, i, i, m - i + 1);

            #if MT == C32 || MT == C64 {
                _labrd_lacgv(i - 1, *y_i1_row, ldy, MT);
            }
            a_i1_block_raw := _labrd_matrix_storage(a.*, lda, i, 1, m - i + 1, i - 1);
            a_i1_block := array2d_view(a_i1_block_raw, lda, row_major=false);
            _labrd_gemv(#char "N", m - i + 1, i - 1, minus_one_t, a_i1_block, lda, y_i1_row, ldy, one_t, *a_i_col_tail_raw, 1, MT);
            #if MT == C32 || MT == C64 {
                _labrd_lacgv(i - 1, *y_i1_row, ldy, MT);
            }

            x_i1_row := _labrd_matrix_row(x.*, ldx, i, 1, i - 1);
            a_1i_col := _labrd_matrix_column(a.*, lda, 1, i, i - 1);
            x_i1_block_raw := _labrd_matrix_storage(x.*, ldx, i, 1, m - i + 1, i - 1);
            x_i1_block := array2d_view(x_i1_block_raw, ldx, row_major=false);
            _labrd_gemv(#char "N", m - i + 1, i - 1, minus_one_t, x_i1_block, ldx, a_1i_col, 1, one_t, *a_i_col_tail_raw, 1, MT);

            #if MT == C32 || MT == C64 {
                alpha = fget(a.*, i, i);
                col_tail := _labrd_matrix_column(a.*, lda, min(i + 1, m), i, max(0, m - i));
                _labrd_larfg(m - i + 1, *alpha, *col_tail, 1, tauq, i, MT);
                fset(d, i, _labrd_real_part(alpha, RT, MT));
            } else {
                col_tail := _labrd_matrix_column(a.*, lda, min(i + 1, m), i, max(0, m - i));
                _labrd_larfg(m - i + 1, _labrd_matrix_entry_ptr(a, lda, i, i), *col_tail, 1, tauq, i, MT);
                fset(d, i, fget(a.*, i, i));
            }

            if i < n {
                fset(a, i, i, one_t);

                a_i_ip1_block_raw := _labrd_matrix_storage(a.*, lda, i, i + 1, m - i + 1, n - i);
                a_i_ip1_block := array2d_view(a_i_ip1_block_raw, lda, row_major=false);
                y_ip1_i_col := _labrd_matrix_column(y.*, ldy, i + 1, i, n - i);
                _labrd_gemv(trans_conj, m - i + 1, n - i, one_t, a_i_ip1_block, lda, a_i_col_tail_raw, 1, zero_t, *y_ip1_i_col, 1, MT);

                a_i_1_block_raw := _labrd_matrix_storage(a.*, lda, i, 1, m - i + 1, i - 1);
                a_i_1_block := array2d_view(a_i_1_block_raw, lda, row_major=false);
                y_1_i_col := _labrd_matrix_column(y.*, ldy, 1, i, i - 1);
                _labrd_gemv(trans_conj, m - i + 1, i - 1, one_t, a_i_1_block, lda, a_i_col_tail_raw, 1, zero_t, *y_1_i_col, 1, MT);

                y_ip1_1_block_raw := _labrd_matrix_storage(y.*, ldy, i + 1, 1, n - i, i - 1);
                y_ip1_1_block := array2d_view(y_ip1_1_block_raw, ldy, row_major=false);
                _labrd_gemv(#char "N", n - i, i - 1, minus_one_t, y_ip1_1_block, ldy, y_1_i_col, 1, one_t, *y_ip1_i_col, 1, MT);

                x_i_1_block_raw := _labrd_matrix_storage(x.*, ldx, i, 1, m - i + 1, i - 1);
                x_i_1_block := array2d_view(x_i_1_block_raw, ldx, row_major=false);
                _labrd_gemv(trans_conj, m - i + 1, i - 1, one_t, x_i_1_block, ldx, a_i_col_tail_raw, 1, zero_t, *y_1_i_col, 1, MT);

                a_1_ip1_block_raw := _labrd_matrix_storage(a.*, lda, 1, i + 1, i - 1, n - i);
                a_1_ip1_block := array2d_view(a_1_ip1_block_raw, lda, row_major=false);
                _labrd_gemv(trans_conj, i - 1, n - i, minus_one_t, a_1_ip1_block, lda, y_1_i_col, 1, one_t, *y_ip1_i_col, 1, MT);

                _labrd_scal(n - i, fget(tauq.*, i), y_ip1_i_col, 1, MT);

                a_i_ip1_row := _labrd_matrix_row(a.*, lda, i, i + 1, n - i);
                #if MT == C32 || MT == C64 {
                    a_i_1_row := _labrd_matrix_row(a.*, lda, i, 1, i);
                    _labrd_lacgv(n - i, *a_i_ip1_row, lda, MT);
                    _labrd_lacgv(i, *a_i_1_row, lda, MT);
                }

                y_ip1_1_wide_raw := _labrd_matrix_storage(y.*, ldy, i + 1, 1, n - i, i);
                y_ip1_1_wide := array2d_view(y_ip1_1_wide_raw, ldy, row_major=false);
                a_i_1_row_for_y := _labrd_matrix_row(a.*, lda, i, 1, i);
                _labrd_gemv(#char "N", n - i, i, minus_one_t, y_ip1_1_wide, ldy, a_i_1_row_for_y, lda, one_t, *a_i_ip1_row, lda, MT);

                #if MT == C32 || MT == C64 {
                    _labrd_lacgv(i, *a_i_1_row, lda, MT);
                    x_i_1_row_conj := _labrd_matrix_row(x.*, ldx, i, 1, i - 1);
                    _labrd_lacgv(i - 1, *x_i_1_row_conj, ldx, MT);
                }
                a_1_ip1_block2_raw := _labrd_matrix_storage(a.*, lda, 1, i + 1, i - 1, n - i);
                a_1_ip1_block2 := array2d_view(a_1_ip1_block2_raw, lda, row_major=false);
                x_i_1_row := _labrd_matrix_row(x.*, ldx, i, 1, i - 1);
                _labrd_gemv(trans_conj, i - 1, n - i, minus_one_t, a_1_ip1_block2, lda, x_i_1_row, ldx, one_t, *a_i_ip1_row, lda, MT);
                #if MT == C32 || MT == C64 {
                    _labrd_lacgv(i - 1, *x_i_1_row, ldx, MT);
                }

                #if MT == C32 || MT == C64 {
                    alpha = fget(a.*, i, i + 1);
                    row_tail := _labrd_matrix_row(a.*, lda, i, min(i + 2, n), max(0, n - i - 1));
                    _labrd_larfg(n - i, *alpha, *row_tail, lda, taup, i, MT);
                    fset(e, i, _labrd_real_part(alpha, RT, MT));
                } else {
                    row_tail := _labrd_matrix_row(a.*, lda, i, min(i + 2, n), max(0, n - i - 1));
                    _labrd_larfg(n - i, _labrd_matrix_entry_ptr(a, lda, i, i + 1), *row_tail, lda, taup, i, MT);
                    fset(e, i, fget(a.*, i, i + 1));
                }
                fset(a, i, i + 1, one_t);

                a_ip1_ip1_block_raw := _labrd_matrix_storage(a.*, lda, i + 1, i + 1, m - i, n - i);
                a_ip1_ip1_block := array2d_view(a_ip1_ip1_block_raw, lda, row_major=false);
                a_i_ip1_row2 := _labrd_matrix_row(a.*, lda, i, i + 1, n - i);
                x_ip1_i_col := _labrd_matrix_column(x.*, ldx, i + 1, i, m - i);
                _labrd_gemv(#char "N", m - i, n - i, one_t, a_ip1_ip1_block, lda, a_i_ip1_row2, lda, zero_t, *x_ip1_i_col, 1, MT);

                y_ip1_1_wide2_raw := _labrd_matrix_storage(y.*, ldy, i + 1, 1, n - i, i);
                y_ip1_1_wide2 := array2d_view(y_ip1_1_wide2_raw, ldy, row_major=false);
                x_1_i_col := _labrd_matrix_column(x.*, ldx, 1, i, i);
                _labrd_gemv(trans_conj, n - i, i, one_t, y_ip1_1_wide2, ldy, a_i_ip1_row2, lda, zero_t, *x_1_i_col, 1, MT);

                a_ip1_1_block_raw := _labrd_matrix_storage(a.*, lda, i + 1, 1, m - i, i);
                a_ip1_1_block := array2d_view(a_ip1_1_block_raw, lda, row_major=false);
                _labrd_gemv(#char "N", m - i, i, minus_one_t, a_ip1_1_block, lda, x_1_i_col, 1, one_t, *x_ip1_i_col, 1, MT);

                a_1_ip1_short_raw := _labrd_matrix_storage(a.*, lda, 1, i + 1, i - 1, n - i);
                a_1_ip1_short := array2d_view(a_1_ip1_short_raw, lda, row_major=false);
                x_1_i_short := _labrd_matrix_column(x.*, ldx, 1, i, i - 1);
                _labrd_gemv(#char "N", i - 1, n - i, one_t, a_1_ip1_short, lda, a_i_ip1_row2, lda, zero_t, *x_1_i_short, 1, MT);

                x_ip1_1_short_raw := _labrd_matrix_storage(x.*, ldx, i + 1, 1, m - i, i - 1);
                x_ip1_1_short := array2d_view(x_ip1_1_short_raw, ldx, row_major=false);
                _labrd_gemv(#char "N", m - i, i - 1, minus_one_t, x_ip1_1_short, ldx, x_1_i_short, 1, one_t, *x_ip1_i_col, 1, MT);

                _labrd_scal(m - i, fget(taup.*, i), x_ip1_i_col, 1, MT);

                #if MT == C32 || MT == C64 {
                    _labrd_lacgv(n - i, *a_i_ip1_row2, lda, MT);
                }
            }
        }
    } else {
        for i:1..nb {
            a_i_row := _labrd_matrix_row(a.*, lda, i, i, n - i + 1);
            #if MT == C32 || MT == C64 {
                a_i_1_row_conj := _labrd_matrix_row(a.*, lda, i, 1, i - 1);
                _labrd_lacgv(n - i + 1, *a_i_row, lda, MT);
                _labrd_lacgv(i - 1, *a_i_1_row_conj, lda, MT);
            }

            y_i_1_block_raw := _labrd_matrix_storage(y.*, ldy, i, 1, n - i + 1, i - 1);
            y_i_1_block := array2d_view(y_i_1_block_raw, ldy, row_major=false);
            a_i_1_row := _labrd_matrix_row(a.*, lda, i, 1, i - 1);
            _labrd_gemv(#char "N", n - i + 1, i - 1, minus_one_t, y_i_1_block, ldy, a_i_1_row, lda, one_t, *a_i_row, lda, MT);

            #if MT == C32 || MT == C64 {
                _labrd_lacgv(i - 1, *a_i_1_row, lda, MT);
                x_i_1_row_conj := _labrd_matrix_row(x.*, ldx, i, 1, i - 1);
                _labrd_lacgv(i - 1, *x_i_1_row_conj, ldx, MT);
            }
            a_1_i_block_raw := _labrd_matrix_storage(a.*, lda, 1, i, i - 1, n - i + 1);
            a_1_i_block := array2d_view(a_1_i_block_raw, lda, row_major=false);
            x_i_1_row := _labrd_matrix_row(x.*, ldx, i, 1, i - 1);
            _labrd_gemv(trans_conj, i - 1, n - i + 1, minus_one_t, a_1_i_block, lda, x_i_1_row, ldx, one_t, *a_i_row, lda, MT);
            #if MT == C32 || MT == C64 {
                _labrd_lacgv(i - 1, *x_i_1_row, ldx, MT);
            }

            #if MT == C32 || MT == C64 {
                alpha = fget(a.*, i, i);
                row_tail := _labrd_matrix_row(a.*, lda, i, min(i + 1, n), max(0, n - i));
                _labrd_larfg(n - i + 1, *alpha, *row_tail, lda, taup, i, MT);
                fset(d, i, _labrd_real_part(alpha, RT, MT));
            } else {
                row_tail := _labrd_matrix_row(a.*, lda, i, min(i + 1, n), max(0, n - i));
                _labrd_larfg(n - i + 1, _labrd_matrix_entry_ptr(a, lda, i, i), *row_tail, lda, taup, i, MT);
                fset(d, i, fget(a.*, i, i));
            }

            if i < m {
                fset(a, i, i, one_t);

                a_ip1_i_block_raw := _labrd_matrix_storage(a.*, lda, i + 1, i, m - i, n - i + 1);
                a_ip1_i_block := array2d_view(a_ip1_i_block_raw, lda, row_major=false);
                x_ip1_i_col := _labrd_matrix_column(x.*, ldx, i + 1, i, m - i);
                _labrd_gemv(#char "N", m - i, n - i + 1, one_t, a_ip1_i_block, lda, a_i_row, lda, zero_t, *x_ip1_i_col, 1, MT);

                y_i_1_block2_raw := _labrd_matrix_storage(y.*, ldy, i, 1, n - i + 1, i - 1);
                y_i_1_block2 := array2d_view(y_i_1_block2_raw, ldy, row_major=false);
                x_1_i_col := _labrd_matrix_column(x.*, ldx, 1, i, i - 1);
                _labrd_gemv(trans_conj, n - i + 1, i - 1, one_t, y_i_1_block2, ldy, a_i_row, lda, zero_t, *x_1_i_col, 1, MT);

                a_ip1_1_block_raw := _labrd_matrix_storage(a.*, lda, i + 1, 1, m - i, i - 1);
                a_ip1_1_block := array2d_view(a_ip1_1_block_raw, lda, row_major=false);
                _labrd_gemv(#char "N", m - i, i - 1, minus_one_t, a_ip1_1_block, lda, x_1_i_col, 1, one_t, *x_ip1_i_col, 1, MT);

                a_1_i_block2_raw := _labrd_matrix_storage(a.*, lda, 1, i, i - 1, n - i + 1);
                a_1_i_block2 := array2d_view(a_1_i_block2_raw, lda, row_major=false);
                x_1_i_col2 := _labrd_matrix_column(x.*, ldx, 1, i, i - 1);
                _labrd_gemv(#char "N", i - 1, n - i + 1, one_t, a_1_i_block2, lda, a_i_row, lda, zero_t, *x_1_i_col2, 1, MT);

                x_ip1_1_block_raw := _labrd_matrix_storage(x.*, ldx, i + 1, 1, m - i, i - 1);
                x_ip1_1_block := array2d_view(x_ip1_1_block_raw, ldx, row_major=false);
                _labrd_gemv(#char "N", m - i, i - 1, minus_one_t, x_ip1_1_block, ldx, x_1_i_col2, 1, one_t, *x_ip1_i_col, 1, MT);

                _labrd_scal(m - i, fget(taup.*, i), x_ip1_i_col, 1, MT);

                #if MT == C32 || MT == C64 {
                    _labrd_lacgv(n - i + 1, *a_i_row, lda, MT);
                }

                y_i_1_row := _labrd_matrix_row(y.*, ldy, i, 1, i - 1);
                #if MT == C32 || MT == C64 {
                    _labrd_lacgv(i - 1, *y_i_1_row, ldy, MT);
                }
                a_ip1_1_block2_raw := _labrd_matrix_storage(a.*, lda, i + 1, 1, m - i, i - 1);
                a_ip1_1_block2 := array2d_view(a_ip1_1_block2_raw, lda, row_major=false);
                a_ip1_i_col := _labrd_matrix_column(a.*, lda, i + 1, i, m - i);
                _labrd_gemv(#char "N", m - i, i - 1, minus_one_t, a_ip1_1_block2, lda, y_i_1_row, ldy, one_t, *a_ip1_i_col, 1, MT);
                #if MT == C32 || MT == C64 {
                    _labrd_lacgv(i - 1, *y_i_1_row, ldy, MT);
                }

                x_ip1_1_wide_raw := _labrd_matrix_storage(x.*, ldx, i + 1, 1, m - i, i);
                x_ip1_1_wide := array2d_view(x_ip1_1_wide_raw, ldx, row_major=false);
                a_1_i_col := _labrd_matrix_column(a.*, lda, 1, i, i);
                _labrd_gemv(#char "N", m - i, i, minus_one_t, x_ip1_1_wide, ldx, a_1_i_col, 1, one_t, *a_ip1_i_col, 1, MT);

                #if MT == C32 || MT == C64 {
                    alpha = fget(a.*, i + 1, i);
                    col_tail := _labrd_matrix_column(a.*, lda, min(i + 2, m), i, max(0, m - i - 1));
                    _labrd_larfg(m - i, *alpha, *col_tail, 1, tauq, i, MT);
                    fset(e, i, _labrd_real_part(alpha, RT, MT));
                } else {
                    col_tail := _labrd_matrix_column(a.*, lda, min(i + 2, m), i, max(0, m - i - 1));
                    _labrd_larfg(m - i, _labrd_matrix_entry_ptr(a, lda, i + 1, i), *col_tail, 1, tauq, i, MT);
                    fset(e, i, fget(a.*, i + 1, i));
                }
                fset(a, i + 1, i, one_t);

                a_ip1_ip1_block_raw := _labrd_matrix_storage(a.*, lda, i + 1, i + 1, m - i, n - i);
                a_ip1_ip1_block := array2d_view(a_ip1_ip1_block_raw, lda, row_major=false);
                a_ip1_i_col2 := _labrd_matrix_column(a.*, lda, i + 1, i, m - i);
                y_ip1_i_col := _labrd_matrix_column(y.*, ldy, i + 1, i, n - i);
                _labrd_gemv(trans_conj, m - i, n - i, one_t, a_ip1_ip1_block, lda, a_ip1_i_col2, 1, zero_t, *y_ip1_i_col, 1, MT);

                a_ip1_1_block3_raw := _labrd_matrix_storage(a.*, lda, i + 1, 1, m - i, i - 1);
                a_ip1_1_block3 := array2d_view(a_ip1_1_block3_raw, lda, row_major=false);
                y_1_i_col := _labrd_matrix_column(y.*, ldy, 1, i, i - 1);
                _labrd_gemv(trans_conj, m - i, i - 1, one_t, a_ip1_1_block3, lda, a_ip1_i_col2, 1, zero_t, *y_1_i_col, 1, MT);

                y_ip1_1_block3_raw := _labrd_matrix_storage(y.*, ldy, i + 1, 1, n - i, i - 1);
                y_ip1_1_block3 := array2d_view(y_ip1_1_block3_raw, ldy, row_major=false);
                _labrd_gemv(#char "N", n - i, i - 1, minus_one_t, y_ip1_1_block3, ldy, y_1_i_col, 1, one_t, *y_ip1_i_col, 1, MT);

                x_ip1_1_wide2_raw := _labrd_matrix_storage(x.*, ldx, i + 1, 1, m - i, i);
                x_ip1_1_wide2 := array2d_view(x_ip1_1_wide2_raw, ldx, row_major=false);
                _labrd_gemv(trans_conj, m - i, i, one_t, x_ip1_1_wide2, ldx, a_ip1_i_col2, 1, zero_t, *y_1_i_col, 1, MT);

                a_1_ip1_block_raw := _labrd_matrix_storage(a.*, lda, 1, i + 1, i, n - i);
                a_1_ip1_block := array2d_view(a_1_ip1_block_raw, lda, row_major=false);
                _labrd_gemv(trans_conj, i, n - i, minus_one_t, a_1_ip1_block, lda, y_1_i_col, 1, one_t, *y_ip1_i_col, 1, MT);

                _labrd_scal(n - i, fget(tauq.*, i), y_ip1_i_col, 1, MT);
            } else {
                #if MT == C32 || MT == C64 {
                    _labrd_lacgv(n - i + 1, *a_i_row, lda, MT);
                }
            }
        }
    }
}
