// ? LAPACK *hegst()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chegst/zhegst
// *HEGST reduces a complex Hermitian-definite generalized eigenproblem
// to standard form.

chegst :: #bake_arguments _hegst(RT=float32, MT=C32);
zhegst :: #bake_arguments _hegst(RT=float64, MT=C64);

#scope_file
_hegst_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hegst_matrix_storage :: inline (
    a: $A/ArrayType2D,
    $$lda: int,
    $$row: int,
    $$col: int,
    $$m: int,
    $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hegst_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hegst_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hegst_hegs2 :: (
    $$itype: int,
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == C32 {
        chegs2(itype, uplo, n, a, lda, b, ldb, info);
    } else {
        zhegs2(itype, uplo, n, a, lda, b, ldb, info);
    }
}

_hegst_hemm :: (
    $$side: u8,
    $$uplo: u8,
    $$m: int,
    $$n: int,
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
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
} {
    #if MT == C32 {
        chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_hegst_her2k :: (
    $$uplo: u8,
    $$trans: u8,
    $$n: int,
    $$k: int,
    alpha: $MT,
    a: $A/ArrayType2D,
    $$lda: int,
    b: $B/ArrayType2D,
    $$ldb: int,
    beta: $RT,
    c: *$C/ArrayType2D,
    $$ldc: int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT &&
           C.data_type == MT;
} {
    #if MT == C32 {
        cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}

_hegst_trmm :: (
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
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == C32 {
        ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_hegst_trsm :: (
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
    return (MT == C32 || MT == C64) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    #if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_hegst :: (
    $$itype: int,
    $$uplo: u8,
    $$n: int,
    a: *$A/ArrayType2D,
    $$lda: int,
    b: *$B/ArrayType2D,
    $$ldb: int,
    info: *int,
    $RT: Type,
    $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT &&
           B.data_type == MT;
} {
    one_r  :: #run one(RT);
    cone_t :: #run _hegst_complex_from_real(one_r, RT, MT);
    half_t :: #run _hegst_complex_from_real(cast(RT) 0.5, RT, MT);
    routine_name :: #run if MT == C32 { "CHEGST" } else { "ZHEGST" };

    upper : bool = ---;
    k, kb, nb : int = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if itype < 1 || itype > 3 {
        info.* = -1;
    } else if !(upper || lsame(uplo, #char "L")) {
        info.* = -2;
    } else if n < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ldb < max(1, n) {
        info.* = -7;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 then return;

    nb = ilaenv(1, ifx MT == C32 then "CHEGST" else "ZHEGST", cast(*u8) *uplo, n, -1, -1, -1);
    if nb <= 1 || nb >= n {
        _hegst_hegs2(itype, uplo, n, a, lda, b, ldb, info, MT);
        return;
    }

    if itype == 1 {
        if upper {
            k = 1;
            while k <= n {
                kb = min(n - k + 1, nb);

                a_kk_raw := _hegst_matrix_storage(a.*, lda, k, k, kb, kb);
                b_kk_raw := _hegst_matrix_storage(b.*, ldb, k, k, kb, kb);
                a_kk := array2d_view(a_kk_raw, lda, row_major=false);
                b_kk := array2d_view(b_kk_raw, ldb, row_major=false);
                _hegst_hegs2(itype, uplo, kb, *a_kk, lda, *b_kk, ldb, info, MT);

                if k + kb <= n {
                    ntrail := n - k - kb + 1;
                    a_k12_raw := _hegst_matrix_storage(a.*, lda, k, k + kb, kb, ntrail);
                    b_k12_raw := _hegst_matrix_storage(b.*, ldb, k, k + kb, kb, ntrail);
                    a_22_raw  := _hegst_matrix_storage(a.*, lda, k + kb, k + kb, ntrail, ntrail);
                    b_22_raw  := _hegst_matrix_storage(b.*, ldb, k + kb, k + kb, ntrail, ntrail);
                    a_k12 := array2d_view(a_k12_raw, lda, row_major=false);
                    b_k12 := array2d_view(b_k12_raw, ldb, row_major=false);
                    a_22 := array2d_view(a_22_raw, lda, row_major=false);
                    b_22 := array2d_view(b_22_raw, ldb, row_major=false);

                    _hegst_trsm(#char "L", uplo, #char "C", #char "N", kb, ntrail, cone_t, b_kk, ldb, *a_k12, lda, MT);
                    _hegst_hemm(#char "L", uplo, kb, ntrail, -half_t, a_kk, lda, b_k12, ldb, cone_t, *a_k12, lda, MT);
                    _hegst_her2k(uplo, #char "C", ntrail, kb, -cone_t, a_k12, lda, b_k12, ldb, one_r, *a_22, lda, RT, MT);
                    _hegst_hemm(#char "L", uplo, kb, ntrail, -half_t, a_kk, lda, b_k12, ldb, cone_t, *a_k12, lda, MT);
                    _hegst_trsm(#char "R", uplo, #char "N", #char "N", kb, ntrail, cone_t, b_22, ldb, *a_k12, lda, MT);
                }

                k += nb;
            }
        } else {
            k = 1;
            while k <= n {
                kb = min(n - k + 1, nb);

                a_kk_raw := _hegst_matrix_storage(a.*, lda, k, k, kb, kb);
                b_kk_raw := _hegst_matrix_storage(b.*, ldb, k, k, kb, kb);
                a_kk := array2d_view(a_kk_raw, lda, row_major=false);
                b_kk := array2d_view(b_kk_raw, ldb, row_major=false);
                _hegst_hegs2(itype, uplo, kb, *a_kk, lda, *b_kk, ldb, info, MT);

                if k + kb <= n {
                    ntrail := n - k - kb + 1;
                    a_21_raw := _hegst_matrix_storage(a.*, lda, k + kb, k, ntrail, kb);
                    b_21_raw := _hegst_matrix_storage(b.*, ldb, k + kb, k, ntrail, kb);
                    a_22_raw := _hegst_matrix_storage(a.*, lda, k + kb, k + kb, ntrail, ntrail);
                    b_22_raw := _hegst_matrix_storage(b.*, ldb, k + kb, k + kb, ntrail, ntrail);
                    a_21 := array2d_view(a_21_raw, lda, row_major=false);
                    b_21 := array2d_view(b_21_raw, ldb, row_major=false);
                    a_22 := array2d_view(a_22_raw, lda, row_major=false);
                    b_22 := array2d_view(b_22_raw, ldb, row_major=false);

                    _hegst_trsm(#char "R", uplo, #char "C", #char "N", ntrail, kb, cone_t, b_kk, ldb, *a_21, lda, MT);
                    _hegst_hemm(#char "R", uplo, ntrail, kb, -half_t, a_kk, lda, b_21, ldb, cone_t, *a_21, lda, MT);
                    _hegst_her2k(uplo, #char "N", ntrail, kb, -cone_t, a_21, lda, b_21, ldb, one_r, *a_22, lda, RT, MT);
                    _hegst_hemm(#char "R", uplo, ntrail, kb, -half_t, a_kk, lda, b_21, ldb, cone_t, *a_21, lda, MT);
                    _hegst_trsm(#char "L", uplo, #char "N", #char "N", ntrail, kb, cone_t, b_22, ldb, *a_21, lda, MT);
                }

                k += nb;
            }
        }
    } else {
        if upper {
            k = 1;
            while k <= n {
                kb = min(n - k + 1, nb);

                a_kk_raw := _hegst_matrix_storage(a.*, lda, k, k, kb, kb);
                b_kk_raw := _hegst_matrix_storage(b.*, ldb, k, k, kb, kb);
                a_kk := array2d_view(a_kk_raw, lda, row_major=false);
                b_kk := array2d_view(b_kk_raw, ldb, row_major=false);

                if k > 1 {
                    kprev := k - 1;
                    a_1k_raw := _hegst_matrix_storage(a.*, lda, 1, k, kprev, kb);
                    b_1k_raw := _hegst_matrix_storage(b.*, ldb, 1, k, kprev, kb);
                    a_11_raw := _hegst_matrix_storage(a.*, lda, 1, 1, kprev, kprev);
                    b_11_raw := _hegst_matrix_storage(b.*, ldb, 1, 1, kprev, kprev);
                    a_1k := array2d_view(a_1k_raw, lda, row_major=false);
                    b_1k := array2d_view(b_1k_raw, ldb, row_major=false);
                    a_11 := array2d_view(a_11_raw, lda, row_major=false);
                    b_11 := array2d_view(b_11_raw, ldb, row_major=false);

                    _hegst_trmm(#char "L", uplo, #char "N", #char "N", kprev, kb, cone_t, b_11, ldb, *a_1k, lda, MT);
                    _hegst_hemm(#char "R", uplo, kprev, kb, half_t, a_kk, lda, b_1k, ldb, cone_t, *a_1k, lda, MT);
                    _hegst_her2k(uplo, #char "N", kprev, kb, cone_t, a_1k, lda, b_1k, ldb, one_r, *a_11, lda, RT, MT);
                    _hegst_hemm(#char "R", uplo, kprev, kb, half_t, a_kk, lda, b_1k, ldb, cone_t, *a_1k, lda, MT);
                    _hegst_trmm(#char "R", uplo, #char "C", #char "N", kprev, kb, cone_t, b_kk, ldb, *a_1k, lda, MT);
                }

                _hegst_hegs2(itype, uplo, kb, *a_kk, lda, *b_kk, ldb, info, MT);
                k += nb;
            }
        } else {
            k = 1;
            while k <= n {
                kb = min(n - k + 1, nb);

                a_kk_raw := _hegst_matrix_storage(a.*, lda, k, k, kb, kb);
                b_kk_raw := _hegst_matrix_storage(b.*, ldb, k, k, kb, kb);
                a_kk := array2d_view(a_kk_raw, lda, row_major=false);
                b_kk := array2d_view(b_kk_raw, ldb, row_major=false);

                if k > 1 {
                    kprev := k - 1;
                    a_k1_raw := _hegst_matrix_storage(a.*, lda, k, 1, kb, kprev);
                    b_k1_raw := _hegst_matrix_storage(b.*, ldb, k, 1, kb, kprev);
                    a_11_raw := _hegst_matrix_storage(a.*, lda, 1, 1, kprev, kprev);
                    b_11_raw := _hegst_matrix_storage(b.*, ldb, 1, 1, kprev, kprev);
                    a_k1 := array2d_view(a_k1_raw, lda, row_major=false);
                    b_k1 := array2d_view(b_k1_raw, ldb, row_major=false);
                    a_11 := array2d_view(a_11_raw, lda, row_major=false);
                    b_11 := array2d_view(b_11_raw, ldb, row_major=false);

                    _hegst_trmm(#char "R", uplo, #char "N", #char "N", kb, kprev, cone_t, b_11, ldb, *a_k1, lda, MT);
                    _hegst_hemm(#char "L", uplo, kb, kprev, half_t, a_kk, lda, b_k1, ldb, cone_t, *a_k1, lda, MT);
                    _hegst_her2k(uplo, #char "C", kprev, kb, cone_t, a_k1, lda, b_k1, ldb, one_r, *a_11, lda, RT, MT);
                    _hegst_hemm(#char "L", uplo, kb, kprev, half_t, a_kk, lda, b_k1, ldb, cone_t, *a_k1, lda, MT);
                    _hegst_trmm(#char "L", uplo, #char "C", #char "N", kb, kprev, cone_t, b_kk, ldb, *a_k1, lda, MT);
                }

                _hegst_hegs2(itype, uplo, kb, *a_kk, lda, *b_kk, ldb, info, MT);
                k += nb;
            }
        }
    }
}
