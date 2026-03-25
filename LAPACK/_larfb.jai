// ? LAPACK *larfb()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarfb/dlarfb/clarfb/zlarfb
// *LARFB applies a block reflector to a matrix.

slarfb :: #bake_arguments _larfb(RT=float32, MT=float32);
dlarfb :: #bake_arguments _larfb(RT=float64, MT=float64);
clarfb :: #bake_arguments _larfb(RT=float32, MT=C32);
zlarfb :: #bake_arguments _larfb(RT=float64, MT=C64);

#scope_file
_larfb_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_larfb_matrix_storage :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _larfb_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_larfb_matrix_storage_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _larfb_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_larfb_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int, alpha: $MT,
    a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int, beta: $MT,
    c: *$C/ArrayType2D, $$ldc: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT && C.data_type == MT;
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

_larfb_trmm :: (
    $$side: u8, $$uplo: u8, $$transa: u8, $$diag: u8, $$m: int, $$n: int, alpha: $MT,
    a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT;
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

_larfb_copy_row_to_work :: (
    c: *$C/ArrayType2D, $$row: int, work: *$WORK/ArrayType2D, $$col: int, $$n: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           C.data_type == MT && WORK.data_type == MT;
} {
    for i:1..n {
        val := fget(c.*, row, i);
        #if MT == C32 || MT == C64 {
            val = conjugate(val);
        }
        fset(work, i, col, val);
    }
}

_larfb_copy_col_to_work :: (
    c: *$C/ArrayType2D, $$col: int, work: *$WORK/ArrayType2D, $$wcol: int, $$m: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           C.data_type == MT && WORK.data_type == MT;
} {
    for i:1..m {
        fset(work, i, wcol, fget(c.*, i, col));
    }
}

_larfb_left_subtract_value :: inline (x: $MT, $MT: Type) -> MT #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    #if MT == C32 || MT == C64 {
        return conjugate(x);
    } else {
        return x;
    }
}

_larfb :: (
    $$side: u8, $$trans: u8, $$direct: u8, $$storev: u8, $$m: int, $$n: int, $$k: int,
    v: $V/ArrayType2D, $$ldv: int, t: $T/ArrayType2D, $$ldt: int, c: *$C/ArrayType2D, $$ldc: int,
    work: *$WORK/ArrayType2D, $$ldwork: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT && T.data_type == MT && C.data_type == MT && WORK.data_type == MT;
} {
    one_t :: #run one(MT);
    neg_one_t : MT = -one_t;
    transh :: u8 = if MT == C32 || MT == C64 { #char "C" } else { #char "T" };
    transt : u8 = if lsame(trans, #char "N") { transh } else { #char "N" };

    if m <= 0 || n <= 0 || k <= 0 then return;

    t_raw := _larfb_matrix_storage_view(t, ldt, 1, 1, k, k);
    t_view := array2d_view(t_raw, ldt, row_major=false);

    if lsame(storev, #char "C") {
        if lsame(direct, #char "F") {
            if lsame(side, #char "L") {
                for j:1..k {
                    _larfb_copy_row_to_work(c, j, work, j, n, MT);
                }

                v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, k, k);
                v1 := array2d_view(v1_raw, ldv, row_major=false);
                _larfb_trmm(#char "R", #char "L", #char "N", #char "U", n, k, one_t, v1.*, ldv, work, ldwork, MT);

                if m > k {
                    c2_raw := _larfb_matrix_storage(c, ldc, k + 1, 1, m - k, n);
                    c2 := array2d_view(c2_raw, ldc, row_major=false);
                    v2_raw := _larfb_matrix_storage_view(v, ldv, k + 1, 1, m - k, k);
                    v2 := array2d_view(v2_raw, ldv, row_major=false);
                    _larfb_gemm(transh, #char "N", n, k, m - k, one_t, c2.*, ldc, v2.*, ldv, one_t, work, ldwork, MT);
                }

                _larfb_trmm(#char "R", #char "U", transt, #char "N", n, k, one_t, t_view.*, ldt, work, ldwork, MT);

                if m > k {
                    c2_raw := _larfb_matrix_storage(c, ldc, k + 1, 1, m - k, n);
                    c2 := array2d_view(c2_raw, ldc, row_major=false);
                    v2_raw := _larfb_matrix_storage_view(v, ldv, k + 1, 1, m - k, k);
                    v2 := array2d_view(v2_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", transh, m - k, n, k, neg_one_t, v2.*, ldv, work.*, ldwork, one_t, *c2, ldc, MT);
                }

                _larfb_trmm(#char "R", #char "L", transh, #char "U", n, k, one_t, v1.*, ldv, work, ldwork, MT);

                for j:1..k {
                    for i:1..n {
                        fset(c, j, i, fget(c.*, j, i) - _larfb_left_subtract_value(fget(work.*, i, j), MT));
                    }
                }
            } else if lsame(side, #char "R") {
                for j:1..k {
                    _larfb_copy_col_to_work(c, j, work, j, m, MT);
                }

                v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, k, k);
                v1 := array2d_view(v1_raw, ldv, row_major=false);
                _larfb_trmm(#char "R", #char "L", #char "N", #char "U", m, k, one_t, v1.*, ldv, work, ldwork, MT);

                if n > k {
                    c2_raw := _larfb_matrix_storage(c, ldc, 1, k + 1, m, n - k);
                    c2 := array2d_view(c2_raw, ldc, row_major=false);
                    v2_raw := _larfb_matrix_storage_view(v, ldv, k + 1, 1, n - k, k);
                    v2 := array2d_view(v2_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", #char "N", m, k, n - k, one_t, c2.*, ldc, v2.*, ldv, one_t, work, ldwork, MT);
                }

                _larfb_trmm(#char "R", #char "U", trans, #char "N", m, k, one_t, t_view.*, ldt, work, ldwork, MT);

                if n > k {
                    c2_raw := _larfb_matrix_storage(c, ldc, 1, k + 1, m, n - k);
                    c2 := array2d_view(c2_raw, ldc, row_major=false);
                    v2_raw := _larfb_matrix_storage_view(v, ldv, k + 1, 1, n - k, k);
                    v2 := array2d_view(v2_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", transh, m, n - k, k, neg_one_t, work.*, ldwork, v2.*, ldv, one_t, *c2, ldc, MT);
                }

                _larfb_trmm(#char "R", #char "L", transh, #char "U", m, k, one_t, v1.*, ldv, work, ldwork, MT);

                for j:1..k {
                    for i:1..m {
                        fset(c, i, j, fget(c.*, i, j) - fget(work.*, i, j));
                    }
                }
            }
        } else {
            if lsame(side, #char "L") {
                for j:1..k {
                    _larfb_copy_row_to_work(c, m - k + j, work, j, n, MT);
                }

                v2_raw := _larfb_matrix_storage_view(v, ldv, m - k + 1, 1, k, k);
                v2 := array2d_view(v2_raw, ldv, row_major=false);
                _larfb_trmm(#char "R", #char "U", #char "N", #char "U", n, k, one_t, v2.*, ldv, work, ldwork, MT);

                if m > k {
                    c1_raw := _larfb_matrix_storage(c, ldc, 1, 1, m - k, n);
                    c1 := array2d_view(c1_raw, ldc, row_major=false);
                    v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, m - k, k);
                    v1 := array2d_view(v1_raw, ldv, row_major=false);
                    _larfb_gemm(transh, #char "N", n, k, m - k, one_t, c1.*, ldc, v1.*, ldv, one_t, work, ldwork, MT);
                }

                _larfb_trmm(#char "R", #char "L", transt, #char "N", n, k, one_t, t_view.*, ldt, work, ldwork, MT);

                if m > k {
                    c1_raw := _larfb_matrix_storage(c, ldc, 1, 1, m - k, n);
                    c1 := array2d_view(c1_raw, ldc, row_major=false);
                    v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, m - k, k);
                    v1 := array2d_view(v1_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", transh, m - k, n, k, neg_one_t, v1.*, ldv, work.*, ldwork, one_t, *c1, ldc, MT);
                }

                _larfb_trmm(#char "R", #char "U", transh, #char "U", n, k, one_t, v2.*, ldv, work, ldwork, MT);

                for j:1..k {
                    for i:1..n {
                        fset(c, m - k + j, i, fget(c.*, m - k + j, i) - _larfb_left_subtract_value(fget(work.*, i, j), MT));
                    }
                }
            } else if lsame(side, #char "R") {
                for j:1..k {
                    _larfb_copy_col_to_work(c, n - k + j, work, j, m, MT);
                }

                v2_raw := _larfb_matrix_storage_view(v, ldv, n - k + 1, 1, k, k);
                v2 := array2d_view(v2_raw, ldv, row_major=false);
                _larfb_trmm(#char "R", #char "U", #char "N", #char "U", m, k, one_t, v2.*, ldv, work, ldwork, MT);

                if n > k {
                    c1_raw := _larfb_matrix_storage(c, ldc, 1, 1, m, n - k);
                    c1 := array2d_view(c1_raw, ldc, row_major=false);
                    v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, n - k, k);
                    v1 := array2d_view(v1_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", #char "N", m, k, n - k, one_t, c1.*, ldc, v1.*, ldv, one_t, work, ldwork, MT);
                }

                _larfb_trmm(#char "R", #char "L", trans, #char "N", m, k, one_t, t_view.*, ldt, work, ldwork, MT);

                if n > k {
                    c1_raw := _larfb_matrix_storage(c, ldc, 1, 1, m, n - k);
                    c1 := array2d_view(c1_raw, ldc, row_major=false);
                    v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, n - k, k);
                    v1 := array2d_view(v1_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", transh, m, n - k, k, neg_one_t, work.*, ldwork, v1.*, ldv, one_t, *c1, ldc, MT);
                }

                _larfb_trmm(#char "R", #char "U", transh, #char "U", m, k, one_t, v2.*, ldv, work, ldwork, MT);

                for j:1..k {
                    for i:1..m {
                        fset(c, i, n - k + j, fget(c.*, i, n - k + j) - fget(work.*, i, j));
                    }
                }
            }
        }
    } else if lsame(storev, #char "R") {
        if lsame(direct, #char "F") {
            if lsame(side, #char "L") {
                for j:1..k {
                    _larfb_copy_row_to_work(c, j, work, j, n, MT);
                }

                v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, k, k);
                v1 := array2d_view(v1_raw, ldv, row_major=false);
                _larfb_trmm(#char "R", #char "U", transh, #char "U", n, k, one_t, v1.*, ldv, work, ldwork, MT);

                if m > k {
                    c2_raw := _larfb_matrix_storage(c, ldc, k + 1, 1, m - k, n);
                    c2 := array2d_view(c2_raw, ldc, row_major=false);
                    v2_raw := _larfb_matrix_storage_view(v, ldv, 1, k + 1, k, m - k);
                    v2 := array2d_view(v2_raw, ldv, row_major=false);
                    _larfb_gemm(transh, transh, n, k, m - k, one_t, c2.*, ldc, v2.*, ldv, one_t, work, ldwork, MT);
                }

                _larfb_trmm(#char "R", #char "U", transt, #char "N", n, k, one_t, t_view.*, ldt, work, ldwork, MT);

                if m > k {
                    c2_raw := _larfb_matrix_storage(c, ldc, k + 1, 1, m - k, n);
                    c2 := array2d_view(c2_raw, ldc, row_major=false);
                    v2_raw := _larfb_matrix_storage_view(v, ldv, 1, k + 1, k, m - k);
                    v2 := array2d_view(v2_raw, ldv, row_major=false);
                    _larfb_gemm(transh, transh, m - k, n, k, neg_one_t, v2.*, ldv, work.*, ldwork, one_t, *c2, ldc, MT);
                }

                _larfb_trmm(#char "R", #char "U", #char "N", #char "U", n, k, one_t, v1.*, ldv, work, ldwork, MT);

                for j:1..k {
                    for i:1..n {
                        fset(c, j, i, fget(c.*, j, i) - _larfb_left_subtract_value(fget(work.*, i, j), MT));
                    }
                }
            } else if lsame(side, #char "R") {
                for j:1..k {
                    _larfb_copy_col_to_work(c, j, work, j, m, MT);
                }

                v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, k, k);
                v1 := array2d_view(v1_raw, ldv, row_major=false);
                _larfb_trmm(#char "R", #char "U", transh, #char "U", m, k, one_t, v1.*, ldv, work, ldwork, MT);

                if n > k {
                    c2_raw := _larfb_matrix_storage(c, ldc, 1, k + 1, m, n - k);
                    c2 := array2d_view(c2_raw, ldc, row_major=false);
                    v2_raw := _larfb_matrix_storage_view(v, ldv, 1, k + 1, k, n - k);
                    v2 := array2d_view(v2_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", transh, m, k, n - k, one_t, c2.*, ldc, v2.*, ldv, one_t, work, ldwork, MT);
                }

                _larfb_trmm(#char "R", #char "U", trans, #char "N", m, k, one_t, t_view.*, ldt, work, ldwork, MT);

                if n > k {
                    c2_raw := _larfb_matrix_storage(c, ldc, 1, k + 1, m, n - k);
                    c2 := array2d_view(c2_raw, ldc, row_major=false);
                    v2_raw := _larfb_matrix_storage_view(v, ldv, 1, k + 1, k, n - k);
                    v2 := array2d_view(v2_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", #char "N", m, n - k, k, neg_one_t, work.*, ldwork, v2.*, ldv, one_t, *c2, ldc, MT);
                }

                _larfb_trmm(#char "R", #char "U", #char "N", #char "U", m, k, one_t, v1.*, ldv, work, ldwork, MT);

                for j:1..k {
                    for i:1..m {
                        fset(c, i, j, fget(c.*, i, j) - fget(work.*, i, j));
                    }
                }
            }
        } else {
            if lsame(side, #char "L") {
                for j:1..k {
                    _larfb_copy_row_to_work(c, m - k + j, work, j, n, MT);
                }

                v2_raw := _larfb_matrix_storage_view(v, ldv, 1, m - k + 1, k, k);
                v2 := array2d_view(v2_raw, ldv, row_major=false);
                _larfb_trmm(#char "R", #char "L", transh, #char "U", n, k, one_t, v2.*, ldv, work, ldwork, MT);

                if m > k {
                    c1_raw := _larfb_matrix_storage(c, ldc, 1, 1, m - k, n);
                    c1 := array2d_view(c1_raw, ldc, row_major=false);
                    v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, k, m - k);
                    v1 := array2d_view(v1_raw, ldv, row_major=false);
                    _larfb_gemm(transh, transh, n, k, m - k, one_t, c1.*, ldc, v1.*, ldv, one_t, work, ldwork, MT);
                }

                _larfb_trmm(#char "R", #char "L", transt, #char "N", n, k, one_t, t_view.*, ldt, work, ldwork, MT);

                if m > k {
                    c1_raw := _larfb_matrix_storage(c, ldc, 1, 1, m - k, n);
                    c1 := array2d_view(c1_raw, ldc, row_major=false);
                    v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, k, m - k);
                    v1 := array2d_view(v1_raw, ldv, row_major=false);
                    _larfb_gemm(transh, transh, m - k, n, k, neg_one_t, v1.*, ldv, work.*, ldwork, one_t, *c1, ldc, MT);
                }

                _larfb_trmm(#char "R", #char "L", #char "N", #char "U", n, k, one_t, v2.*, ldv, work, ldwork, MT);

                for j:1..k {
                    for i:1..n {
                        fset(c, m - k + j, i, fget(c.*, m - k + j, i) - _larfb_left_subtract_value(fget(work.*, i, j), MT));
                    }
                }
            } else if lsame(side, #char "R") {
                for j:1..k {
                    _larfb_copy_col_to_work(c, n - k + j, work, j, m, MT);
                }

                v2_raw := _larfb_matrix_storage_view(v, ldv, 1, n - k + 1, k, k);
                v2 := array2d_view(v2_raw, ldv, row_major=false);
                _larfb_trmm(#char "R", #char "L", transh, #char "U", m, k, one_t, v2.*, ldv, work, ldwork, MT);

                if n > k {
                    c1_raw := _larfb_matrix_storage(c, ldc, 1, 1, m, n - k);
                    c1 := array2d_view(c1_raw, ldc, row_major=false);
                    v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, k, n - k);
                    v1 := array2d_view(v1_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", transh, m, k, n - k, one_t, c1.*, ldc, v1.*, ldv, one_t, work, ldwork, MT);
                }

                _larfb_trmm(#char "R", #char "L", trans, #char "N", m, k, one_t, t_view.*, ldt, work, ldwork, MT);

                if n > k {
                    c1_raw := _larfb_matrix_storage(c, ldc, 1, 1, m, n - k);
                    c1 := array2d_view(c1_raw, ldc, row_major=false);
                    v1_raw := _larfb_matrix_storage_view(v, ldv, 1, 1, k, n - k);
                    v1 := array2d_view(v1_raw, ldv, row_major=false);
                    _larfb_gemm(#char "N", #char "N", m, n - k, k, neg_one_t, work.*, ldwork, v1.*, ldv, one_t, *c1, ldc, MT);
                }

                _larfb_trmm(#char "R", #char "L", #char "N", #char "U", m, k, one_t, v2.*, ldv, work, ldwork, MT);

                for j:1..k {
                    for i:1..m {
                        fset(c, i, n - k + j, fget(c.*, i, n - k + j) - fget(work.*, i, j));
                    }
                }
            }
        }
    }
}
