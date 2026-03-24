// ? LAPACK *lahr2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slahr2/dlahr2/clahr2/zlahr2
// *LAHR2 reduces the specified number of first columns of a general
// matrix and returns the auxiliary block-reflector data.

slahr2 :: #bake_arguments _lahr2(MT=float32);
dlahr2 :: #bake_arguments _lahr2(MT=float64);
clahr2 :: #bake_arguments _lahr2(MT=C32);
zlahr2 :: #bake_arguments _lahr2(MT=C64);

#scope_file
_lahr2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_lahr2_vector_entry_ptr :: inline (x: *$X/ArrayType1D, $$index: int) -> *X.data_type {
    return x.data + (index - 1);
}

_lahr2_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_lahr2_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _lahr2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lahr2_matrix_row :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _lahr2_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_lahr2_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _lahr2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_lahr2_axpy :: (
    $$n: int,
    alpha: $MT,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        saxpy(n, alpha, x, incx, y, incy);
    } else if MT == float64 {
        daxpy(n, alpha, x, incx, y, incy);
    } else if MT == C32 {
        caxpy(n, alpha, x, incx, y, incy);
    } else {
        zaxpy(n, alpha, x, incx, y, incy);
    }
}

_lahr2_copy :: (
    $$n: int,
    x: $X/ArrayType1D,
    $$incx: int,
    y: *$Y/ArrayType1D,
    $$incy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT &&
           Y.data_type == MT;
} {
    #if MT == float32 {
        scopy(n, x, incx, y, incy);
    } else if MT == float64 {
        dcopy(n, x, incx, y, incy);
    } else if MT == C32 {
        ccopy(n, x, incx, y, incy);
    } else {
        zcopy(n, x, incx, y, incy);
    }
}

_lahr2_gemm :: (
    $$transa: u8,
    $$transb: u8,
    $$m: int,
    $$n: int,
    $$k: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: MT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
} {
    #if MT == float32 {
        sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == float64 {
        dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if MT == C32 {
        cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_lahr2_gemv :: (
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

_lahr2_lacgv :: (
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

_lahr2_lacpy :: (
    $$uplo: u8,
    $$m: int,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == float32 {
        slacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == float64 {
        dlacpy(uplo, m, n, a, lda, b, ldb);
    } else if MT == C32 {
        clacpy(uplo, m, n, a, lda, b, ldb);
    } else {
        zlacpy(uplo, m, n, a, lda, b, ldb);
    }
}

_lahr2_larfg :: (
    $$n: int,
    alpha: *$MT,
    x: *$X/ArrayType1D,
    $$incx: int,
    tau: *$MT,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           X.data_type == MT;
} {
    #if MT == float32 {
        slarfg(n, alpha, x, incx, tau);
    } else if MT == float64 {
        dlarfg(n, alpha, x, incx, tau);
    } else if MT == C32 {
        clarfg(n, alpha, x, incx, tau);
    } else {
        zlarfg(n, alpha, x, incx, tau);
    }
}

_lahr2_scal :: (
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

_lahr2_trmm :: (
    $$side: u8,
    $$uplo: u8,
    $$transa: u8,
    $$diag: u8,
    $$m: int,
    $$n: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == float32 {
        strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == float64 {
        dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else if MT == C32 {
        ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_lahr2_trmv :: (
    $$uplo: u8,
    $$trans: u8,
    $$diag: u8,
    $$n: int,
    a: $A/ArrayType2D,
    $$lda: int,
    x: *$X/ArrayType1D,
    $$incx: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           X.data_type == MT;
} {
    #if MT == float32 {
        strmv(uplo, trans, diag, n, a, lda, x, incx);
    } else if MT == float64 {
        dtrmv(uplo, trans, diag, n, a, lda, x, incx);
    } else if MT == C32 {
        ctrmv(uplo, trans, diag, n, a, lda, x, incx);
    } else {
        ztrmv(uplo, trans, diag, n, a, lda, x, incx);
    }
}

_lahr2 :: (
    $$n: int,
    $$k: int,
    $$nb: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    tau: *$TAU/ArrayType1D,
    t: *$T/ArrayType2D,
    $$ldt: int,
    y: *$Y/ArrayType2D,
    $$ldy: int,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT &&
           TAU.data_type == MT &&
           T.data_type == MT &&
           Y.data_type == MT;
} {
    zero_t :: #run zero(MT);
    one_t :: #run one(MT);
    minus_one_t :: #run -one_t;
    trans_conj :: #run if MT == float32 || MT == float64 { #char "T" } else { #char "C" };

    ei : MT = ---;
    i : int = ---;

    if n <= 1 || nb <= 0 then return;

    for i:1..nb {
        if i > 1 {
            a_kip1_row := _lahr2_matrix_row(a.*, lda, k + i - 1, 1, i - 1);
            a_kp1_i_head := _lahr2_matrix_column(a.*, lda, k + 1, i, i - 1);
            a_kp1_i_tail := _lahr2_matrix_column(a.*, lda, k + 1, i, n - k);
            a_kp1_1_raw := _lahr2_matrix_storage(a.*, lda, k + 1, 1, i - 1, i - 1);
            a_kp1_1 := array2d_view(a_kp1_1_raw, lda, row_major=false);
            a_ki_1_raw := _lahr2_matrix_storage(a.*, lda, k + i, 1, n - k - i + 1, i - 1);
            a_ki_1 := array2d_view(a_ki_1_raw, lda, row_major=false);
            a_ki_i := _lahr2_matrix_column(a.*, lda, k + i, i, n - k - i + 1);
            t_1_nb := _lahr2_matrix_column(t.*, ldt, 1, nb, i - 1);
            t_top_raw := _lahr2_matrix_storage(t.*, ldt, 1, 1, i - 1, i - 1);
            t_top := array2d_view(t_top_raw, ldt, row_major=false);
            y_kp1_1_raw := _lahr2_matrix_storage(y.*, ldy, k + 1, 1, n - k, i - 1);
            y_kp1_1 := array2d_view(y_kp1_1_raw, ldy, row_major=false);

            #if MT == C32 || MT == C64 {
                _lahr2_lacgv(i - 1, *a_kip1_row, lda, MT);
            }
            _lahr2_gemv(#char "N", n - k, i - 1, minus_one_t, y_kp1_1, ldy,
                        a_kip1_row, lda, one_t, *a_kp1_i_tail, 1, MT);
            #if MT == C32 || MT == C64 {
                _lahr2_lacgv(i - 1, *a_kip1_row, lda, MT);
            }

            _lahr2_copy(i - 1, a_kp1_i_head, 1, *t_1_nb, 1, MT);
            _lahr2_trmv(#char "L", trans_conj, #char "U", i - 1, a_kp1_1, lda, *t_1_nb, 1, MT);
            _lahr2_gemv(trans_conj, n - k - i + 1, i - 1, one_t, a_ki_1, lda,
                        a_ki_i, 1, one_t, *t_1_nb, 1, MT);
            _lahr2_trmv(#char "U", trans_conj, #char "N", i - 1, t_top, ldt, *t_1_nb, 1, MT);
            _lahr2_gemv(#char "N", n - k - i + 1, i - 1, minus_one_t, a_ki_1, lda,
                        t_1_nb, 1, one_t, *a_ki_i, 1, MT);
            _lahr2_trmv(#char "L", #char "N", #char "U", i - 1, a_kp1_1, lda, *t_1_nb, 1, MT);
            _lahr2_axpy(i - 1, minus_one_t, t_1_nb, 1, *a_kp1_i_head, 1, MT);

            fset(a, k + i - 1, i - 1, ei);
        }

        alpha_ptr := _lahr2_matrix_entry_ptr(a, lda, k + i, i);
        col_tail := _lahr2_matrix_column(a.*, lda, min(k + i + 1, n), i, max(0, n - k - i));
        _lahr2_larfg(n - k - i + 1, alpha_ptr, *col_tail, 1, _lahr2_vector_entry_ptr(tau, i), MT);
        ei = fget(a.*, k + i, i);
        fset(a, k + i, i, one_t);

        a_kp1_ip1_raw := _lahr2_matrix_storage(a.*, lda, k + 1, i + 1, n - k, n - k - i + 1);
        a_kp1_ip1 := array2d_view(a_kp1_ip1_raw, lda, row_major=false);
        a_ki_i := _lahr2_matrix_column(a.*, lda, k + i, i, n - k - i + 1);
        y_kp1_i := _lahr2_matrix_column(y.*, ldy, k + 1, i, n - k);
        _lahr2_gemv(#char "N", n - k, n - k - i + 1, one_t, a_kp1_ip1, lda,
                    a_ki_i, 1, zero_t, *y_kp1_i, 1, MT);

        a_ki_1_raw := _lahr2_matrix_storage(a.*, lda, k + i, 1, n - k - i + 1, i - 1);
        a_ki_1 := array2d_view(a_ki_1_raw, lda, row_major=false);
        t_1_i := _lahr2_matrix_column(t.*, ldt, 1, i, i - 1);
        _lahr2_gemv(trans_conj, n - k - i + 1, i - 1, one_t, a_ki_1, lda,
                    a_ki_i, 1, zero_t, *t_1_i, 1, MT);

        y_kp1_1_raw := _lahr2_matrix_storage(y.*, ldy, k + 1, 1, n - k, i - 1);
        y_kp1_1 := array2d_view(y_kp1_1_raw, ldy, row_major=false);
        _lahr2_gemv(#char "N", n - k, i - 1, minus_one_t, y_kp1_1, ldy,
                    t_1_i, 1, one_t, *y_kp1_i, 1, MT);
        _lahr2_scal(n - k, fget(tau.*, i), *y_kp1_i, 1, MT);

        _lahr2_scal(i - 1, -fget(tau.*, i), *t_1_i, 1, MT);
        t_top_raw := _lahr2_matrix_storage(t.*, ldt, 1, 1, i - 1, i - 1);
        t_top := array2d_view(t_top_raw, ldt, row_major=false);
        _lahr2_trmv(#char "U", #char "N", #char "N", i - 1, t_top, ldt, *t_1_i, 1, MT);
        fset(t, i, i, fget(tau.*, i));
    }

    fset(a, k + nb, nb, ei);

    y_top_raw := _lahr2_matrix_storage(y.*, ldy, 1, 1, k, nb);
    y_top := array2d_view(y_top_raw, ldy, row_major=false);
    a_top_raw := _lahr2_matrix_storage(a.*, lda, 1, 2, k, nb);
    a_top := array2d_view(a_top_raw, lda, row_major=false);
    _lahr2_lacpy(#char "A", k, nb, a_top, lda, *y_top, ldy, MT);

    a_v_raw := _lahr2_matrix_storage(a.*, lda, k + 1, 1, nb, nb);
    a_v := array2d_view(a_v_raw, lda, row_major=false);
    _lahr2_trmm(#char "R", #char "L", #char "N", #char "U", k, nb, one_t, a_v, lda, *y_top, ldy, MT);

    if n > k + nb {
        a_top_rest_raw := _lahr2_matrix_storage(a.*, lda, 1, 2 + nb, k, n - k - nb);
        a_top_rest := array2d_view(a_top_rest_raw, lda, row_major=false);
        a_bottom_raw := _lahr2_matrix_storage(a.*, lda, k + 1 + nb, 1, n - k - nb, nb);
        a_bottom := array2d_view(a_bottom_raw, lda, row_major=false);
        _lahr2_gemm(#char "N", #char "N", k, nb, n - k - nb, one_t,
                    a_top_rest, lda, a_bottom, lda, one_t, *y_top, ldy, MT);
    }

    t_nb_raw := _lahr2_matrix_storage(t.*, ldt, 1, 1, nb, nb);
    t_nb := array2d_view(t_nb_raw, ldt, row_major=false);
    _lahr2_trmm(#char "R", #char "U", #char "N", #char "N", k, nb, one_t, t_nb, ldt, *y_top, ldy, MT);
}
