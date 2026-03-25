// ? LAPACK *larfb_gett()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarfb_gett/dlarfb_gett/clarfb_gett/zlarfb_gett
// *LARFB_GETT is an auxiliary block-reflector helper.

slarfb_gett :: #bake_arguments _larfb_gett(RT=float32, MT=float32);
dlarfb_gett :: #bake_arguments _larfb_gett(RT=float64, MT=float64);
clarfb_gett :: #bake_arguments _larfb_gett(RT=float32, MT=C32);
zlarfb_gett :: #bake_arguments _larfb_gett(RT=float64, MT=C64);

#scope_file
_larfb_gett_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_larfb_gett_matrix_storage :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _larfb_gett_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_larfb_gett_matrix_storage_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _larfb_gett_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_larfb_gett_gemm :: (
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

_larfb_gett_trmm :: (
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

_larfb_gett :: (
    $$ident: u8, $$m: int, $$n: int, $$k: int, t: $T/ArrayType2D, $$ldt: int, a: *$A/ArrayType2D, $$lda: int,
    b: *$B/ArrayType2D, $$ldb: int, work: *$WORK/ArrayType2D, $$ldwork: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           T.data_type == MT && A.data_type == MT && B.data_type == MT && WORK.data_type == MT;
} {
    one_t :: #run one(MT);
    zero_t :: #run zero(MT);
    neg_one_t : MT = -one_t;
    lnotident := !lsame(ident, #char "I");
    trans_h :: u8 = if MT == C32 || MT == C64 { #char "C" } else { #char "T" };

    if m < 0 || n <= 0 || k == 0 || k > n then return;

    if n > k {
        for j:1..n - k {
            for i:1..k {
                fset(work, i, j, fget(a.*, i, k + j));
            }
        }

        w2_raw := _larfb_gett_matrix_storage(work, ldwork, 1, 1, k, n - k);
        w2 := array2d_view(w2_raw, ldwork, row_major=false);

        if lnotident {
            a1_raw := _larfb_gett_matrix_storage(a, lda, 1, 1, k, k);
            a1 := array2d_view(a1_raw, lda, row_major=false);
            _larfb_gett_trmm(#char "L", #char "L", trans_h, #char "U", k, n - k, one_t, a1.*, lda, *w2, ldwork, MT);
        }

        if m > 0 {
            b1_raw := _larfb_gett_matrix_storage(b, ldb, 1, 1, m, k);
            b1 := array2d_view(b1_raw, ldb, row_major=false);
            b2_raw := _larfb_gett_matrix_storage(b, ldb, 1, k + 1, m, n - k);
            b2 := array2d_view(b2_raw, ldb, row_major=false);
            _larfb_gett_gemm(trans_h, #char "N", k, n - k, m, one_t, b1.*, ldb, b2.*, ldb, one_t, *w2, ldwork, MT);
        }

        t_raw := _larfb_gett_matrix_storage_view(t, ldt, 1, 1, k, k);
        t_view := array2d_view(t_raw, ldt, row_major=false);
        _larfb_gett_trmm(#char "L", #char "U", #char "N", #char "N", k, n - k, one_t, t_view.*, ldt, *w2, ldwork, MT);

        if m > 0 {
            b1_raw := _larfb_gett_matrix_storage(b, ldb, 1, 1, m, k);
            b1 := array2d_view(b1_raw, ldb, row_major=false);
            b2_raw := _larfb_gett_matrix_storage(b, ldb, 1, k + 1, m, n - k);
            b2 := array2d_view(b2_raw, ldb, row_major=false);
            _larfb_gett_gemm(#char "N", #char "N", m, n - k, k, neg_one_t, b1.*, ldb, w2.*, ldwork, one_t, *b2, ldb, MT);
        }

        if lnotident {
            a1_raw := _larfb_gett_matrix_storage(a, lda, 1, 1, k, k);
            a1 := array2d_view(a1_raw, lda, row_major=false);
            _larfb_gett_trmm(#char "L", #char "L", #char "N", #char "U", k, n - k, one_t, a1.*, lda, *w2, ldwork, MT);
        }

        for j:1..n - k {
            for i:1..k {
                fset(a, i, k + j, fget(a.*, i, k + j) - fget(work.*, i, j));
            }
        }
    }

    for j:1..k {
        for i:1..j {
            fset(work, i, j, fget(a.*, i, j));
        }
    }
    for j:1..k - 1 {
        for i:j + 1..k {
            fset(work, i, j, zero_t);
        }
    }

    w1_raw := _larfb_gett_matrix_storage(work, ldwork, 1, 1, k, k);
    w1 := array2d_view(w1_raw, ldwork, row_major=false);

    if lnotident {
        a1_raw := _larfb_gett_matrix_storage(a, lda, 1, 1, k, k);
        a1 := array2d_view(a1_raw, lda, row_major=false);
        _larfb_gett_trmm(#char "L", #char "L", trans_h, #char "U", k, k, one_t, a1.*, lda, *w1, ldwork, MT);
    }

    t_raw2 := _larfb_gett_matrix_storage_view(t, ldt, 1, 1, k, k);
    t_view2 := array2d_view(t_raw2, ldt, row_major=false);
    _larfb_gett_trmm(#char "L", #char "U", #char "N", #char "N", k, k, one_t, t_view2.*, ldt, *w1, ldwork, MT);

    if m > 0 {
        _larfb_gett_trmm(#char "R", #char "U", #char "N", #char "N", m, k, neg_one_t, w1.*, ldwork, b, ldb, MT);
    }

    if lnotident {
        a1_raw := _larfb_gett_matrix_storage(a, lda, 1, 1, k, k);
        a1 := array2d_view(a1_raw, lda, row_major=false);
        _larfb_gett_trmm(#char "L", #char "L", #char "N", #char "U", k, k, one_t, a1.*, lda, *w1, ldwork, MT);

        for j:1..k - 1 {
            for i:j + 1..k {
                fset(a, i, j, -fget(work.*, i, j));
            }
        }
    }

    for j:1..k {
        for i:1..j {
            fset(a, i, j, fget(a.*, i, j) - fget(work.*, i, j));
        }
    }
}
