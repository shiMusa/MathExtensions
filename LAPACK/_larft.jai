// ? LAPACK *larft()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment slarft/dlarft/clarft/zlarft
// *LARFT forms the triangular factor of a block reflector.

slarft :: #bake_arguments _larft(RT=float32, MT=float32);
dlarft :: #bake_arguments _larft(RT=float64, MT=float64);
clarft :: #bake_arguments _larft(RT=float32, MT=C32);
zlarft :: #bake_arguments _larft(RT=float64, MT=C64);

#scope_file
_larft_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = max(0, count);
    return view;
}

_larft_matrix_storage :: inline (
    a: *$A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _larft_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_larft_vector_storage :: inline (x: $X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _larft_raw_slice(x.data + offset, count);
}

_larft_trmm :: (
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

_larft_gemm :: (
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

_larft_lacpy :: (
    $$uplo: u8, $$m: int, $$n: int, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == float32 || MT == float64 || MT == C32 || MT == C64) &&
           A.data_type == MT && B.data_type == MT;
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

_larft_conj_if_needed :: inline (x: $MT, $MT: Type) -> MT #modify {
    return MT == float32 || MT == float64 || MT == C32 || MT == C64;
} {
    #if MT == C32 || MT == C64 {
        return conjugate(x);
    } else {
        return x;
    }
}

_larft :: (
    $$direct: u8, $$storev: u8, $$n: int, $$k: int, v: $V/ArrayType2D, $$ldv: int,
    tau: $TAU/ArrayType1D, t: *$T/ArrayType2D, $$ldt: int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           (MT == RT || MT == complex_with_inner_type(RT)) &&
           V.data_type == MT && TAU.data_type == MT && T.data_type == MT;
} {
    one_t     :: #run one(MT);
    zero_t    :: #run zero(MT);
    neg_one_t : MT = -one_t;

    colv, dirf, lq, ql, qr : bool = ---;
    l : int = ---;

    if n == 0 || k == 0 then return;

    if n == 1 || k == 1 {
        fset(t, 1, 1, fget(tau, 1));
        return;
    }

    l = k/2;

    dirf = lsame(direct, #char "F");
    colv = lsame(storev, #char "C");
    qr = dirf && colv;
    lq = dirf && !colv;
    ql = !dirf && colv;

    if qr {
        v11_raw := _larft_matrix_storage(v, ldv, 1, 1, n, l);
        v11 := array2d_view(v11_raw, ldv, row_major=false);
        tau1 := _larft_vector_storage(tau, 0, l);
        t11_raw := _larft_matrix_storage(t, ldt, 1, 1, l, l);
        t11 := array2d_view(t11_raw, ldt, row_major=false);
        _larft(direct, storev, n, l, v11.*, ldv, tau1, *t11, ldt, RT, MT);

        v22_raw := _larft_matrix_storage(v, ldv, l + 1, l + 1, n - l, k - l);
        v22 := array2d_view(v22_raw, ldv, row_major=false);
        tau2 := _larft_vector_storage(tau, l, k - l);
        t22_raw := _larft_matrix_storage(t, ldt, l + 1, l + 1, k - l, k - l);
        t22 := array2d_view(t22_raw, ldt, row_major=false);
        _larft(direct, storev, n - l, k - l, v22.*, ldv, tau2, *t22, ldt, RT, MT);

        for j:1..l {
            for i:1..k - l {
                fset(t, j, l + i, _larft_conj_if_needed(fget(v, l + i, j), MT));
            }
        }

        t12_raw := _larft_matrix_storage(t, ldt, 1, l + 1, l, k - l);
        t12 := array2d_view(t12_raw, ldt, row_major=false);
        v22_tri_raw := _larft_matrix_storage(v, ldv, l + 1, l + 1, k - l, k - l);
        v22_tri := array2d_view(v22_tri_raw, ldv, row_major=false);
        _larft_trmm(#char "R", #char "L", #char "N", #char "U", l, k - l, one_t, v22_tri.*, ldv, *t12, ldt, MT);

        if n > k {
            v31_raw := _larft_matrix_storage(v, ldv, k + 1, 1, n - k, l);
            v31 := array2d_view(v31_raw, ldv, row_major=false);
            v32_raw := _larft_matrix_storage(v, ldv, k + 1, l + 1, n - k, k - l);
            v32 := array2d_view(v32_raw, ldv, row_major=false);
            left_trans :: u8 = if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
            _larft_gemm(left_trans, #char "N", l, k - l, n - k, one_t, v31.*, ldv, v32.*, ldv, one_t, *t12, ldt, MT);
        }

        _larft_trmm(#char "L", #char "U", #char "N", #char "N", l, k - l, neg_one_t, t11.*, ldt, *t12, ldt, MT);
        _larft_trmm(#char "R", #char "U", #char "N", #char "N", l, k - l, one_t, t22.*, ldt, *t12, ldt, MT);
    } else if lq {
        v11_raw := _larft_matrix_storage(v, ldv, 1, 1, l, n);
        v11 := array2d_view(v11_raw, ldv, row_major=false);
        tau1 := _larft_vector_storage(tau, 0, l);
        t11_raw := _larft_matrix_storage(t, ldt, 1, 1, l, l);
        t11 := array2d_view(t11_raw, ldt, row_major=false);
        _larft(direct, storev, n, l, v11.*, ldv, tau1, *t11, ldt, RT, MT);

        v22_raw := _larft_matrix_storage(v, ldv, l + 1, l + 1, k - l, n - l);
        v22 := array2d_view(v22_raw, ldv, row_major=false);
        tau2 := _larft_vector_storage(tau, l, k - l);
        t22_raw := _larft_matrix_storage(t, ldt, l + 1, l + 1, k - l, k - l);
        t22 := array2d_view(t22_raw, ldt, row_major=false);
        _larft(direct, storev, n - l, k - l, v22.*, ldv, tau2, *t22, ldt, RT, MT);

        t12_raw := _larft_matrix_storage(t, ldt, 1, l + 1, l, k - l);
        t12 := array2d_view(t12_raw, ldt, row_major=false);
        v12_raw := _larft_matrix_storage(v, ldv, 1, l + 1, l, k - l);
        v12 := array2d_view(v12_raw, ldv, row_major=false);
        _larft_lacpy(#char "A", l, k - l, v12.*, ldv, *t12, ldt, MT);

        v22_tri_raw := _larft_matrix_storage(v, ldv, l + 1, l + 1, k - l, k - l);
        v22_tri := array2d_view(v22_tri_raw, ldv, row_major=false);
        right_trans :: u8 = if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
        _larft_trmm(#char "R", #char "U", right_trans, #char "U", l, k - l, one_t, v22_tri.*, ldv, *t12, ldt, MT);

        if n > k {
            v13_raw := _larft_matrix_storage(v, ldv, 1, k + 1, l, n - k);
            v13 := array2d_view(v13_raw, ldv, row_major=false);
            v23_raw := _larft_matrix_storage(v, ldv, l + 1, k + 1, k - l, n - k);
            v23 := array2d_view(v23_raw, ldv, row_major=false);
            second_trans :: u8 = if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
            _larft_gemm(#char "N", second_trans, l, k - l, n - k, one_t, v13.*, ldv, v23.*, ldv, one_t, *t12, ldt, MT);
        }

        _larft_trmm(#char "L", #char "U", #char "N", #char "N", l, k - l, neg_one_t, t11.*, ldt, *t12, ldt, MT);
        _larft_trmm(#char "R", #char "U", #char "N", #char "N", l, k - l, one_t, t22.*, ldt, *t12, ldt, MT);
    } else if ql {
        v11_raw := _larft_matrix_storage(v, ldv, 1, 1, n - l, k - l);
        v11 := array2d_view(v11_raw, ldv, row_major=false);
        tau1 := _larft_vector_storage(tau, 0, k - l);
        t11_raw := _larft_matrix_storage(t, ldt, 1, 1, k - l, k - l);
        t11 := array2d_view(t11_raw, ldt, row_major=false);
        _larft(direct, storev, n - l, k - l, v11.*, ldv, tau1, *t11, ldt, RT, MT);

        v22_raw := _larft_matrix_storage(v, ldv, 1, k - l + 1, n, l);
        v22 := array2d_view(v22_raw, ldv, row_major=false);
        tau2 := _larft_vector_storage(tau, k - l, l);
        t22_raw := _larft_matrix_storage(t, ldt, k - l + 1, k - l + 1, l, l);
        t22 := array2d_view(t22_raw, ldt, row_major=false);
        _larft(direct, storev, n, l, v22.*, ldv, tau2, *t22, ldt, RT, MT);

        for j:1..k - l {
            for i:1..l {
                fset(t, k - l + i, j, _larft_conj_if_needed(fget(v, n - k + j, k - l + i), MT));
            }
        }

        t21_raw := _larft_matrix_storage(t, ldt, k - l + 1, 1, l, k - l);
        t21 := array2d_view(t21_raw, ldt, row_major=false);
        v21_tri_raw := _larft_matrix_storage(v, ldv, n - k + 1, 1, k - l, k - l);
        v21_tri := array2d_view(v21_tri_raw, ldv, row_major=false);
        _larft_trmm(#char "R", #char "U", #char "N", #char "U", l, k - l, one_t, v21_tri.*, ldv, *t21, ldt, MT);

        if n > k {
            v12_raw := _larft_matrix_storage(v, ldv, 1, k - l + 1, n - k, l);
            v12 := array2d_view(v12_raw, ldv, row_major=false);
            v11_top_raw := _larft_matrix_storage(v, ldv, 1, 1, n - k, k - l);
            v11_top := array2d_view(v11_top_raw, ldv, row_major=false);
            left_trans :: u8 = if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
            _larft_gemm(left_trans, #char "N", l, k - l, n - k, one_t, v12.*, ldv, v11_top.*, ldv, one_t, *t21, ldt, MT);
        }

        _larft_trmm(#char "L", #char "L", #char "N", #char "N", l, k - l, neg_one_t, t22.*, ldt, *t21, ldt, MT);
        _larft_trmm(#char "R", #char "L", #char "N", #char "N", l, k - l, one_t, t11.*, ldt, *t21, ldt, MT);
    } else {
        v11_raw := _larft_matrix_storage(v, ldv, 1, 1, k - l, n - l);
        v11 := array2d_view(v11_raw, ldv, row_major=false);
        tau1 := _larft_vector_storage(tau, 0, k - l);
        t11_raw := _larft_matrix_storage(t, ldt, 1, 1, k - l, k - l);
        t11 := array2d_view(t11_raw, ldt, row_major=false);
        _larft(direct, storev, n - l, k - l, v11.*, ldv, tau1, *t11, ldt, RT, MT);

        v22_raw := _larft_matrix_storage(v, ldv, k - l + 1, 1, l, n);
        v22 := array2d_view(v22_raw, ldv, row_major=false);
        tau2 := _larft_vector_storage(tau, k - l, l);
        t22_raw := _larft_matrix_storage(t, ldt, k - l + 1, k - l + 1, l, l);
        t22 := array2d_view(t22_raw, ldt, row_major=false);
        _larft(direct, storev, n, l, v22.*, ldv, tau2, *t22, ldt, RT, MT);

        t21_raw := _larft_matrix_storage(t, ldt, k - l + 1, 1, l, k - l);
        t21 := array2d_view(t21_raw, ldt, row_major=false);
        v22_block_raw := _larft_matrix_storage(v, ldv, k - l + 1, n - k + 1, l, k - l);
        v22_block := array2d_view(v22_block_raw, ldv, row_major=false);
        _larft_lacpy(#char "A", l, k - l, v22_block.*, ldv, *t21, ldt, MT);

        right_trans :: u8 = if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
        v12_tri_raw := _larft_matrix_storage(v, ldv, 1, n - k + 1, k - l, k - l);
        v12_tri := array2d_view(v12_tri_raw, ldv, row_major=false);
        _larft_trmm(#char "R", #char "L", right_trans, #char "U", l, k - l, one_t, v12_tri.*, ldv, *t21, ldt, MT);

        if n > k {
            v21_raw := _larft_matrix_storage(v, ldv, k - l + 1, 1, l, n - k);
            v21 := array2d_view(v21_raw, ldv, row_major=false);
            v11_left_raw := _larft_matrix_storage(v, ldv, 1, 1, k - l, n - k);
            v11_left := array2d_view(v11_left_raw, ldv, row_major=false);
            second_trans :: u8 = if MT == float32 || MT == float64 { #char "T" } else { #char "C" };
            _larft_gemm(#char "N", second_trans, l, k - l, n - k, one_t, v21.*, ldv, v11_left.*, ldv, one_t, *t21, ldt, MT);
        }

        _larft_trmm(#char "L", #char "L", #char "N", #char "N", l, k - l, neg_one_t, t22.*, ldt, *t21, ldt, MT);
        _larft_trmm(#char "R", #char "L", #char "N", #char "N", l, k - l, one_t, t11.*, ldt, *t21, ldt, MT);
    }
}
