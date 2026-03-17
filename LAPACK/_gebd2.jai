// ? LAPACK *gebd2()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment sgebd2/cgebd2/dgebd2/zgebd2
// *GEBD2 reduces a general matrix to bidiagonal form
// using an unblocked algorithm.

sgebd2 :: #bake_arguments _gebd2(RT=float32, MT=float32);
dgebd2 :: #bake_arguments _gebd2(RT=float64, MT=float64);
cgebd2 :: #bake_arguments _gebd2(RT=float32, MT=C32);
zgebd2 :: #bake_arguments _gebd2(RT=float64, MT=C64);

#scope_file
_gebd2_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = count;
    return view;
}

_gebd2_matrix_column :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    return _gebd2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gebd2_matrix_row :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$count: int
) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _gebd2_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_gebd2_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _gebd2_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_gebd2_matrix_entry_ptr :: inline (
    a: *$A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int
) -> *A.data_type {
    return a.data + (col - 1)*lda + (row - 1);
}

_gebd2_vector_entry_ptr :: inline (x: *$X/ArrayType1D, $$index: int) -> *X.data_type {
    return x.data + (index - 1);
}

_gebd2_real_part :: inline (x: $MT, $RT: Type, $MT: Type) -> RT #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT));
} {
    #if MT == float32 || MT == float64 {
        return x;
    } else {
        return real(x);
    }
}

_gebd2_store_real :: inline (
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

_gebd2_larfg :: (
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

_gebd2_larf1f :: (
    $$side: u8,
    $$m: int,
    $$n: int,
    v: $V/ArrayType1D,
    $$incv: int,
    tau: $MT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    work: *$WORK/ArrayType1D,
    $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           V.data_type == MT &&
           C.data_type == MT &&
           WORK.data_type == MT;
} {
    #if MT == float32 {
        slarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == float64 {
        dlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else if MT == C32 {
        clarf1f(side, m, n, v, incv, tau, c, ldc, work);
    } else {
        zlarf1f(side, m, n, v, incv, tau, c, ldc, work);
    }
}

_gebd2_lacgv :: (
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

_gebd2 :: (
    $$m: int,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    d: *$D/ArrayType1D,
    e: *$E/ArrayType1D,
    tauq: *$TAUQ/ArrayType1D,
    taup: *$TAUP/ArrayType1D,
    work: *$WORK/ArrayType1D,
    info: *int,
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
           WORK.data_type == MT;
} {
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == float32 {
        "SGEBD2"
    } else if MT == float64 {
        "DGEBD2"
    } else if MT == C32 {
        "CGEBD2"
    } else {
        "ZGEBD2"
    };

    alpha : MT = ---;
    i : int = ---;

    info.* = 0;
    if m < 0 {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, m) {
        info.* = -4;
    }
    if info.* < 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if m >= n {
        for i:1..n {
            #if MT == C32 || MT == C64 {
                alpha = fget(a.*, i, i);
                col_tail := _gebd2_matrix_column(a.*, lda, min(i + 1, m), i, max(0, m - i));
                _gebd2_larfg(m - i + 1, *alpha, *col_tail, 1, _gebd2_vector_entry_ptr(tauq, i), MT);
                fset(d, i, _gebd2_real_part(alpha, RT, MT));
            } else {
                col_tail := _gebd2_matrix_column(a.*, lda, min(i + 1, m), i, max(0, m - i));
                _gebd2_larfg(m - i + 1, _gebd2_matrix_entry_ptr(a, lda, i, i), *col_tail, 1,
                             _gebd2_vector_entry_ptr(tauq, i), MT);
                fset(d, i, _gebd2_real_part(fget(a.*, i, i), RT, MT));
            }

            if i < n {
                v := _gebd2_matrix_column(a.*, lda, i, i, m - i + 1);
                c_raw := _gebd2_matrix_storage(a.*, lda, i, i + 1, m - i + 1, n - i);
                c_view := array2d_view(c_raw, lda, row_major=false);
                _gebd2_larf1f(#char "L", m - i + 1, n - i, v, 1, conjugate(fget(tauq.*, i)),
                              *c_view, lda, work, MT);
            }

            #if MT == C32 || MT == C64 {
                _gebd2_store_real(a, i, i, fget(d.*, i), RT, MT);
            }

            if i < n {
                #if MT == C32 || MT == C64 {
                    row_vec := _gebd2_matrix_row(a.*, lda, i, i + 1, n - i);
                    _gebd2_lacgv(n - i, *row_vec, lda, MT);
                    alpha = fget(a.*, i, i + 1);
                    row_tail := _gebd2_matrix_row(a.*, lda, i, min(i + 2, n), max(0, n - i - 1));
                    _gebd2_larfg(n - i, *alpha, *row_tail, lda, _gebd2_vector_entry_ptr(taup, i), MT);
                    fset(e, i, _gebd2_real_part(alpha, RT, MT));

                    c_raw := _gebd2_matrix_storage(a.*, lda, i + 1, i + 1, m - i, n - i);
                    c_view := array2d_view(c_raw, lda, row_major=false);
                    _gebd2_larf1f(#char "R", m - i, n - i, row_vec, lda, fget(taup.*, i),
                                  *c_view, lda, work, MT);
                    _gebd2_lacgv(n - i, *row_vec, lda, MT);
                    _gebd2_store_real(a, i, i + 1, fget(e.*, i), RT, MT);
                } else {
                    row_vec := _gebd2_matrix_row(a.*, lda, i, i + 1, n - i);
                    row_tail := _gebd2_matrix_row(a.*, lda, i, min(i + 2, n), max(0, n - i - 1));
                    _gebd2_larfg(n - i, _gebd2_matrix_entry_ptr(a, lda, i, i + 1), *row_tail, lda,
                                 _gebd2_vector_entry_ptr(taup, i), MT);
                    fset(e, i, _gebd2_real_part(fget(a.*, i, i + 1), RT, MT));

                    c_raw := _gebd2_matrix_storage(a.*, lda, i + 1, i + 1, m - i, n - i);
                    c_view := array2d_view(c_raw, lda, row_major=false);
                    _gebd2_larf1f(#char "R", m - i, n - i, row_vec, lda, fget(taup.*, i),
                                  *c_view, lda, work, MT);
                }
            } else {
                fset(taup, i, zero_t);
            }
        }
    } else {
        for i:1..m {
            #if MT == C32 || MT == C64 {
                row_vec := _gebd2_matrix_row(a.*, lda, i, i, n - i + 1);
                _gebd2_lacgv(n - i + 1, *row_vec, lda, MT);
                alpha = fget(a.*, i, i);
                row_tail := _gebd2_matrix_row(a.*, lda, i, min(i + 1, n), max(0, n - i));
                _gebd2_larfg(n - i + 1, *alpha, *row_tail, lda, _gebd2_vector_entry_ptr(taup, i), MT);
                fset(d, i, _gebd2_real_part(alpha, RT, MT));
            } else {
                row_tail := _gebd2_matrix_row(a.*, lda, i, min(i + 1, n), max(0, n - i));
                _gebd2_larfg(n - i + 1, _gebd2_matrix_entry_ptr(a, lda, i, i), *row_tail, lda,
                             _gebd2_vector_entry_ptr(taup, i), MT);
                fset(d, i, _gebd2_real_part(fget(a.*, i, i), RT, MT));
            }

            if i < m {
                v := _gebd2_matrix_row(a.*, lda, i, i, n - i + 1);
                c_raw := _gebd2_matrix_storage(a.*, lda, i + 1, i, m - i, n - i + 1);
                c_view := array2d_view(c_raw, lda, row_major=false);
                _gebd2_larf1f(#char "R", m - i, n - i + 1, v, lda, fget(taup.*, i),
                              *c_view, lda, work, MT);
            }

            #if MT == C32 || MT == C64 {
                row_vec := _gebd2_matrix_row(a.*, lda, i, i, n - i + 1);
                _gebd2_lacgv(n - i + 1, *row_vec, lda, MT);
                _gebd2_store_real(a, i, i, fget(d.*, i), RT, MT);
            }

            if i < m {
                #if MT == C32 || MT == C64 {
                    alpha = fget(a.*, i + 1, i);
                    col_tail := _gebd2_matrix_column(a.*, lda, min(i + 2, m), i, max(0, m - i - 1));
                    _gebd2_larfg(m - i, *alpha, *col_tail, 1, _gebd2_vector_entry_ptr(tauq, i), MT);
                    fset(e, i, _gebd2_real_part(alpha, RT, MT));

                    v := _gebd2_matrix_column(a.*, lda, i + 1, i, m - i);
                    c_raw := _gebd2_matrix_storage(a.*, lda, i + 1, i + 1, m - i, n - i);
                    c_view := array2d_view(c_raw, lda, row_major=false);
                    _gebd2_larf1f(#char "L", m - i, n - i, v, 1, conjugate(fget(tauq.*, i)),
                                  *c_view, lda, work, MT);
                    _gebd2_store_real(a, i + 1, i, fget(e.*, i), RT, MT);
                } else {
                    v := _gebd2_matrix_column(a.*, lda, i + 1, i, m - i);
                    col_tail := _gebd2_matrix_column(a.*, lda, min(i + 2, m), i, max(0, m - i - 1));
                    _gebd2_larfg(m - i, _gebd2_matrix_entry_ptr(a, lda, i + 1, i), *col_tail, 1,
                                 _gebd2_vector_entry_ptr(tauq, i), MT);
                    fset(e, i, _gebd2_real_part(fget(a.*, i + 1, i), RT, MT));

                    c_raw := _gebd2_matrix_storage(a.*, lda, i + 1, i + 1, m - i, n - i);
                    c_view := array2d_view(c_raw, lda, row_major=false);
                    _gebd2_larf1f(#char "L", m - i, n - i, v, 1, fget(tauq.*, i),
                                  *c_view, lda, work, MT);
                }
            } else {
                fset(tauq, i, zero_t);
            }
        }
    }
}
