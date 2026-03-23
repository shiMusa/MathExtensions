// ? LAPACK *hfrk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chfrk/zhfrk
// *HFRK performs a Hermitian rank-k update for matrices in RFP format.

chfrk :: #bake_arguments _hfrk(RT=float32, MT=C32);
zhfrk :: #bake_arguments _hfrk(RT=float64, MT=C64);

#scope_file
_hfrk_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hfrk_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hfrk_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hfrk_rfp_storage :: inline (
    c: *$C/ArrayType1D, $$offset: int, $$ldc: int, $$m: int, $$n: int
) -> [] C.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*ldc + m;
    return _hfrk_raw_slice(c.data + (offset - 1), count);
}

_hfrk_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hfrk_herk :: (
    $$uplo: u8, $$trans: u8, $$n: int, $$k: int, alpha: $RT,
    a: $A/ArrayType2D, $$lda: int, beta: $RT, c: *$C/ArrayType2D, $$ldc: int,
    $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && C.data_type == MT;
} {
    #if MT == C32 {
        cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
    } else {
        zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
    }
}

_hfrk_gemm :: (
    $$transa: u8, $$transb: u8, $$m: int, $$n: int, $$k: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int, b: $B/ArrayType2D, $$ldb: int,
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

_hfrk :: (
    $$transr: u8, $$uplo: u8, $$trans: u8, $$n: int, $$k: int,
    alpha: $RT, a: $A/ArrayType2D, $$lda: int, beta: RT,
    c: *$C/ArrayType1D, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && C.data_type == MT;
} {
    zero_r :: #run zero(RT);
    one_r  :: #run one(RT);
    zero_t :: #run zero(MT);
    routine_name :: #run if MT == C32 { "CHFRK " } else { "ZHFRK " };

    lower, normaltransr, nisodd, notrans : bool = ---;
    info, j, nrowa, n1, n2, nk : int = ---;
    calpha, cbeta : MT = ---;

    info = 0;
    normaltransr = lsame(transr, #char "N");
    lower = lsame(uplo, #char "L");
    notrans = lsame(trans, #char "N");

    if notrans {
        nrowa = n;
    } else {
        nrowa = k;
    }

    if !(normaltransr || lsame(transr, #char "C")) {
        info = -1;
    } else if !(lower || lsame(uplo, #char "U")) {
        info = -2;
    } else if !(notrans || lsame(trans, #char "C")) {
        info = -3;
    } else if n < 0 {
        info = -4;
    } else if k < 0 {
        info = -5;
    } else if lda < max(1, nrowa) {
        info = -8;
    }
    if info != 0 {
        xerbla(routine_name, -info);
        return;
    }

    if n == 0 || (((alpha == zero_r) || (k == 0)) && (beta == one_r)) {
        return;
    }

    if alpha == zero_r && beta == zero_r {
        for j:1..((n*(n + 1)) / 2) {
            fset(c, j, zero_t);
        }
        return;
    }

    calpha = _hfrk_complex_from_real(alpha, RT, MT);
    cbeta = _hfrk_complex_from_real(beta, RT, MT);

    if (n % 2) == 0 {
        nisodd = false;
        nk = n / 2;
        n1 = 0;
        n2 = 0;
    } else {
        nisodd = true;
        if lower {
            n2 = n / 2;
            n1 = n - n2;
        } else {
            n1 = n / 2;
            n2 = n - n1;
        }
        nk = 0;
    }

    if nisodd {
        if normaltransr {
            if lower {
                if notrans {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, n1, k);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, 1, n, n1, n1);
                    c11 := array2d_view(c11_raw, n, row_major=false);
                    _hfrk_herk(#char "L", #char "N", n1, k, alpha, *a11, lda, beta, *c11, n, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, n1 + 1, 1, n2, k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, n + 1, n, n2, n2);
                    c22 := array2d_view(c22_raw, n, row_major=false);
                    _hfrk_herk(#char "U", #char "N", n2, k, alpha, *a22, lda, beta, *c22, n, RT, MT);

                    c21_raw := _hfrk_rfp_storage(c, n1 + 1, n, n2, n1);
                    c21 := array2d_view(c21_raw, n, row_major=false);
                    _hfrk_gemm(#char "N", #char "C", n2, n1, k, calpha, *a22, lda, *a11, lda, cbeta, *c21, n, MT);
                } else {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, k, n1);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, 1, n, n1, n1);
                    c11 := array2d_view(c11_raw, n, row_major=false);
                    _hfrk_herk(#char "L", #char "C", n1, k, alpha, *a11, lda, beta, *c11, n, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, 1, n1 + 1, k, n2);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, n + 1, n, n2, n2);
                    c22 := array2d_view(c22_raw, n, row_major=false);
                    _hfrk_herk(#char "U", #char "C", n2, k, alpha, *a22, lda, beta, *c22, n, RT, MT);

                    c21_raw := _hfrk_rfp_storage(c, n1 + 1, n, n2, n1);
                    c21 := array2d_view(c21_raw, n, row_major=false);
                    _hfrk_gemm(#char "C", #char "N", n2, n1, k, calpha, *a22, lda, *a11, lda, cbeta, *c21, n, MT);
                }
            } else {
                if notrans {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, n1, k);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, n2 + 1, n, n1, n1);
                    c11 := array2d_view(c11_raw, n, row_major=false);
                    _hfrk_herk(#char "L", #char "N", n1, k, alpha, *a11, lda, beta, *c11, n, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, n2, 1, n2, k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, n1 + 1, n, n2, n2);
                    c22 := array2d_view(c22_raw, n, row_major=false);
                    _hfrk_herk(#char "U", #char "N", n2, k, alpha, *a22, lda, beta, *c22, n, RT, MT);

                    c12_raw := _hfrk_rfp_storage(c, 1, n, n1, n2);
                    c12 := array2d_view(c12_raw, n, row_major=false);
                    _hfrk_gemm(#char "N", #char "C", n1, n2, k, calpha, *a11, lda, *a22, lda, cbeta, *c12, n, MT);
                } else {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, k, n1);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, n2 + 1, n, n1, n1);
                    c11 := array2d_view(c11_raw, n, row_major=false);
                    _hfrk_herk(#char "L", #char "C", n1, k, alpha, *a11, lda, beta, *c11, n, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, 1, n2, k, n2);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, n1 + 1, n, n2, n2);
                    c22 := array2d_view(c22_raw, n, row_major=false);
                    _hfrk_herk(#char "U", #char "C", n2, k, alpha, *a22, lda, beta, *c22, n, RT, MT);

                    c12_raw := _hfrk_rfp_storage(c, 1, n, n1, n2);
                    c12 := array2d_view(c12_raw, n, row_major=false);
                    _hfrk_gemm(#char "C", #char "N", n1, n2, k, calpha, *a11, lda, *a22, lda, cbeta, *c12, n, MT);
                }
            }
        } else {
            if lower {
                if notrans {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, n1, k);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, 1, n1, n1, n1);
                    c11 := array2d_view(c11_raw, n1, row_major=false);
                    _hfrk_herk(#char "U", #char "N", n1, k, alpha, *a11, lda, beta, *c11, n1, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, n1 + 1, 1, n2, k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, 2, n1, n2, n2);
                    c22 := array2d_view(c22_raw, n1, row_major=false);
                    _hfrk_herk(#char "L", #char "N", n2, k, alpha, *a22, lda, beta, *c22, n1, RT, MT);

                    c12_raw := _hfrk_rfp_storage(c, n1*n1 + 1, n1, n1, n2);
                    c12 := array2d_view(c12_raw, n1, row_major=false);
                    _hfrk_gemm(#char "N", #char "C", n1, n2, k, calpha, *a11, lda, *a22, lda, cbeta, *c12, n1, MT);
                } else {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, k, n1);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, 1, n1, n1, n1);
                    c11 := array2d_view(c11_raw, n1, row_major=false);
                    _hfrk_herk(#char "U", #char "C", n1, k, alpha, *a11, lda, beta, *c11, n1, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, 1, n1 + 1, k, n2);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, 2, n1, n2, n2);
                    c22 := array2d_view(c22_raw, n1, row_major=false);
                    _hfrk_herk(#char "L", #char "C", n2, k, alpha, *a22, lda, beta, *c22, n1, RT, MT);

                    c12_raw := _hfrk_rfp_storage(c, n1*n1 + 1, n1, n1, n2);
                    c12 := array2d_view(c12_raw, n1, row_major=false);
                    _hfrk_gemm(#char "C", #char "N", n1, n2, k, calpha, *a11, lda, *a22, lda, cbeta, *c12, n1, MT);
                }
            } else {
                if notrans {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, n1, k);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, n2*n2 + 1, n2, n1, n1);
                    c11 := array2d_view(c11_raw, n2, row_major=false);
                    _hfrk_herk(#char "U", #char "N", n1, k, alpha, *a11, lda, beta, *c11, n2, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, n1 + 1, 1, n2, k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, n1*n2 + 1, n2, n2, n2);
                    c22 := array2d_view(c22_raw, n2, row_major=false);
                    _hfrk_herk(#char "L", #char "N", n2, k, alpha, *a22, lda, beta, *c22, n2, RT, MT);

                    c21_raw := _hfrk_rfp_storage(c, 1, n2, n2, n1);
                    c21 := array2d_view(c21_raw, n2, row_major=false);
                    _hfrk_gemm(#char "N", #char "C", n2, n1, k, calpha, *a22, lda, *a11, lda, cbeta, *c21, n2, MT);
                } else {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, k, n1);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, n2*n2 + 1, n2, n1, n1);
                    c11 := array2d_view(c11_raw, n2, row_major=false);
                    _hfrk_herk(#char "U", #char "C", n1, k, alpha, *a11, lda, beta, *c11, n2, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, 1, n1 + 1, k, n2);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, n1*n2 + 1, n2, n2, n2);
                    c22 := array2d_view(c22_raw, n2, row_major=false);
                    _hfrk_herk(#char "L", #char "C", n2, k, alpha, *a22, lda, beta, *c22, n2, RT, MT);

                    c21_raw := _hfrk_rfp_storage(c, 1, n2, n2, n1);
                    c21 := array2d_view(c21_raw, n2, row_major=false);
                    _hfrk_gemm(#char "C", #char "N", n2, n1, k, calpha, *a22, lda, *a11, lda, cbeta, *c21, n2, MT);
                }
            }
        }
    } else {
        if normaltransr {
            if lower {
                if notrans {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, nk, k);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, 2, n + 1, nk, nk);
                    c11 := array2d_view(c11_raw, n + 1, row_major=false);
                    _hfrk_herk(#char "L", #char "N", nk, k, alpha, *a11, lda, beta, *c11, n + 1, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, nk + 1, 1, nk, k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, 1, n + 1, nk, nk);
                    c22 := array2d_view(c22_raw, n + 1, row_major=false);
                    _hfrk_herk(#char "U", #char "N", nk, k, alpha, *a22, lda, beta, *c22, n + 1, RT, MT);

                    c21_raw := _hfrk_rfp_storage(c, nk + 2, n + 1, nk, nk);
                    c21 := array2d_view(c21_raw, n + 1, row_major=false);
                    _hfrk_gemm(#char "N", #char "C", nk, nk, k, calpha, *a22, lda, *a11, lda, cbeta, *c21, n + 1, MT);
                } else {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, k, nk);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, 2, n + 1, nk, nk);
                    c11 := array2d_view(c11_raw, n + 1, row_major=false);
                    _hfrk_herk(#char "L", #char "C", nk, k, alpha, *a11, lda, beta, *c11, n + 1, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, 1, nk + 1, k, nk);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, 1, n + 1, nk, nk);
                    c22 := array2d_view(c22_raw, n + 1, row_major=false);
                    _hfrk_herk(#char "U", #char "C", nk, k, alpha, *a22, lda, beta, *c22, n + 1, RT, MT);

                    c21_raw := _hfrk_rfp_storage(c, nk + 2, n + 1, nk, nk);
                    c21 := array2d_view(c21_raw, n + 1, row_major=false);
                    _hfrk_gemm(#char "C", #char "N", nk, nk, k, calpha, *a22, lda, *a11, lda, cbeta, *c21, n + 1, MT);
                }
            } else {
                if notrans {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, nk, k);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, nk + 2, n + 1, nk, nk);
                    c11 := array2d_view(c11_raw, n + 1, row_major=false);
                    _hfrk_herk(#char "L", #char "N", nk, k, alpha, *a11, lda, beta, *c11, n + 1, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, nk + 1, 1, nk, k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, nk + 1, n + 1, nk, nk);
                    c22 := array2d_view(c22_raw, n + 1, row_major=false);
                    _hfrk_herk(#char "U", #char "N", nk, k, alpha, *a22, lda, beta, *c22, n + 1, RT, MT);

                    c12_raw := _hfrk_rfp_storage(c, 1, n + 1, nk, nk);
                    c12 := array2d_view(c12_raw, n + 1, row_major=false);
                    _hfrk_gemm(#char "N", #char "C", nk, nk, k, calpha, *a11, lda, *a22, lda, cbeta, *c12, n + 1, MT);
                } else {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, k, nk);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, nk + 2, n + 1, nk, nk);
                    c11 := array2d_view(c11_raw, n + 1, row_major=false);
                    _hfrk_herk(#char "L", #char "C", nk, k, alpha, *a11, lda, beta, *c11, n + 1, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, 1, nk + 1, k, nk);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, nk + 1, n + 1, nk, nk);
                    c22 := array2d_view(c22_raw, n + 1, row_major=false);
                    _hfrk_herk(#char "U", #char "C", nk, k, alpha, *a22, lda, beta, *c22, n + 1, RT, MT);

                    c12_raw := _hfrk_rfp_storage(c, 1, n + 1, nk, nk);
                    c12 := array2d_view(c12_raw, n + 1, row_major=false);
                    _hfrk_gemm(#char "C", #char "N", nk, nk, k, calpha, *a11, lda, *a22, lda, cbeta, *c12, n + 1, MT);
                }
            }
        } else {
            if lower {
                if notrans {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, nk, k);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, nk + 1, nk, nk, nk);
                    c11 := array2d_view(c11_raw, nk, row_major=false);
                    _hfrk_herk(#char "U", #char "N", nk, k, alpha, *a11, lda, beta, *c11, nk, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, nk + 1, 1, nk, k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, 1, nk, nk, nk);
                    c22 := array2d_view(c22_raw, nk, row_major=false);
                    _hfrk_herk(#char "L", #char "N", nk, k, alpha, *a22, lda, beta, *c22, nk, RT, MT);

                    c12_raw := _hfrk_rfp_storage(c, ((nk + 1)*nk) + 1, nk, nk, nk);
                    c12 := array2d_view(c12_raw, nk, row_major=false);
                    _hfrk_gemm(#char "N", #char "C", nk, nk, k, calpha, *a11, lda, *a22, lda, cbeta, *c12, nk, MT);
                } else {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, k, nk);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, nk + 1, nk, nk, nk);
                    c11 := array2d_view(c11_raw, nk, row_major=false);
                    _hfrk_herk(#char "U", #char "C", nk, k, alpha, *a11, lda, beta, *c11, nk, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, 1, nk + 1, k, nk);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, 1, nk, nk, nk);
                    c22 := array2d_view(c22_raw, nk, row_major=false);
                    _hfrk_herk(#char "L", #char "C", nk, k, alpha, *a22, lda, beta, *c22, nk, RT, MT);

                    c12_raw := _hfrk_rfp_storage(c, ((nk + 1)*nk) + 1, nk, nk, nk);
                    c12 := array2d_view(c12_raw, nk, row_major=false);
                    _hfrk_gemm(#char "C", #char "N", nk, nk, k, calpha, *a11, lda, *a22, lda, cbeta, *c12, nk, MT);
                }
            } else {
                if notrans {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, nk, k);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, nk*(nk + 1) + 1, nk, nk, nk);
                    c11 := array2d_view(c11_raw, nk, row_major=false);
                    _hfrk_herk(#char "U", #char "N", nk, k, alpha, *a11, lda, beta, *c11, nk, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, nk + 1, 1, nk, k);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, nk*nk + 1, nk, nk, nk);
                    c22 := array2d_view(c22_raw, nk, row_major=false);
                    _hfrk_herk(#char "L", #char "N", nk, k, alpha, *a22, lda, beta, *c22, nk, RT, MT);

                    c21_raw := _hfrk_rfp_storage(c, 1, nk, nk, nk);
                    c21 := array2d_view(c21_raw, nk, row_major=false);
                    _hfrk_gemm(#char "N", #char "C", nk, nk, k, calpha, *a22, lda, *a11, lda, cbeta, *c21, nk, MT);
                } else {
                    a11_raw := _hfrk_matrix_storage(a, lda, 1, 1, k, nk);
                    a11 := array2d_view(a11_raw, lda, row_major=false);
                    c11_raw := _hfrk_rfp_storage(c, nk*(nk + 1) + 1, nk, nk, nk);
                    c11 := array2d_view(c11_raw, nk, row_major=false);
                    _hfrk_herk(#char "U", #char "C", nk, k, alpha, *a11, lda, beta, *c11, nk, RT, MT);

                    a22_raw := _hfrk_matrix_storage(a, lda, 1, nk + 1, k, nk);
                    a22 := array2d_view(a22_raw, lda, row_major=false);
                    c22_raw := _hfrk_rfp_storage(c, nk*nk + 1, nk, nk, nk);
                    c22 := array2d_view(c22_raw, nk, row_major=false);
                    _hfrk_herk(#char "L", #char "C", nk, k, alpha, *a22, lda, beta, *c22, nk, RT, MT);

                    c21_raw := _hfrk_rfp_storage(c, 1, nk, nk, nk);
                    c21 := array2d_view(c21_raw, nk, row_major=false);
                    _hfrk_gemm(#char "C", #char "N", nk, nk, k, calpha, *a22, lda, *a11, lda, cbeta, *c21, nk, MT);
                }
            }
        }
    }
}
