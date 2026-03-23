// ? LAPACK *hetrf_rk()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrf_rk/zhetrf_rk
// *HETRF_RK computes the blocked bounded Bunch-Kaufman rook factorization
// of a complex Hermitian matrix with a separate E array.

chetrf_rk :: #bake_arguments _hetrf_rk(RT=float32, MT=C32);
zhetrf_rk :: #bake_arguments _hetrf_rk(RT=float64, MT=C64);

#scope_file
_hetrf_rk_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrf_rk_vector_slice :: inline (x: *$X/ArrayType1D, $$offset: int, $$count: int) -> [] X.data_type {
    return _hetrf_rk_raw_slice(x.data + (offset - 1), count);
}

_hetrf_rk_matrix_storage :: inline (data: *$T, $$ld: int, $$rows: int, $$cols: int) -> [] T {
    return _hetrf_rk_raw_slice(data, ifx rows <= 0 || cols <= 0 then 0 else ld*(cols - 1) + rows);
}

_hetrf_rk_matrix_view :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$rows: int, $$cols: int
) -> [] A.data_type {
    return _hetrf_rk_matrix_storage(a.data + (col - 1)*lda + (row - 1), lda, rows, cols);
}

_hetrf_rk_matrix_row :: inline (a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$count: int) -> [] A.data_type {
    stride_count := ifx count <= 0 then 0 else 1 + (count - 1)*lda;
    return _hetrf_rk_raw_slice(a.data + (col - 1)*lda + (row - 1), stride_count);
}

_hetrf_rk_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrf_rk_roundup_lwork :: ($x: $RT, $RT: Type) -> RT #modify {
    return RT == float32 || RT == float64;
} {
    #if RT == float32 { return sroundup_lwork(x); } else { return x; }
}

_hetrf_rk_lahef :: (
    $$uplo: u8, $$n: int, $$nb: int, kb: *int, a: *$A/ArrayType2D, $$lda: int,
    e: *$E/ArrayType1D, ipiv: *$IPIV/ArrayType1D, w: *$W/ArrayType2D,
    $$ldw: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && E.data_type == MT &&
           IPIV.data_type == int && W.data_type == MT;
} {
    #if MT == C32 {
        clahef_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
    } else {
        zlahef_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
    }
}

_hetrf_rk_hetf2 :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    e: *$E/ArrayType1D, ipiv: *$IPIV/ArrayType1D, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) &&
           A.data_type == MT && E.data_type == MT && IPIV.data_type == int;
} {
    #if MT == C32 {
        chetf2_rk(uplo, n, a, lda, e, ipiv, info);
    } else {
        zhetf2_rk(uplo, n, a, lda, e, ipiv, info);
    }
}

_hetrf_rk_swap :: (
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

_hetrf_rk :: (
    $$uplo: u8, $$n: int, a: *$A/ArrayType2D, $$lda: int,
    e: *$E/ArrayType1D, ipiv: *$IPIV/ArrayType1D, work: *$WORK/ArrayType1D,
    $$lwork: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && E.data_type == MT &&
           IPIV.data_type == int && WORK.data_type == MT;
} {
    routine_name :: #run if MT == C32 { "CHETRF_RK" } else { "ZHETRF_RK" };

    lquery, upper : bool = ---;
    i, iinfo, ip, iws, k, kb, ldwork, lwkopt, nb, nbmin : int = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    lquery = lwork == -1;
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if lda < max(1, n) {
        info.* = -4;
    } else if lwork < 1 && !lquery {
        info.* = -8;
    }

    if info.* == 0 {
        nb = ilaenv(1, ifx MT == C32 then "CHETRF_RK" else "ZHETRF_RK", cast(*u8) *uplo, n, -1, -1, -1);
        lwkopt = max(1, n*nb);
        fset(work, 1, _hetrf_rk_complex_from_real(_hetrf_rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
    }

    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    } else if lquery {
        return;
    }

    nbmin = 2;
    ldwork = n;
    if nb > 1 && nb < n {
        iws = ldwork*nb;
        if lwork < iws {
            nb = max(lwork / ldwork, 1);
            nbmin = max(2, ilaenv(2, ifx MT == C32 then "CHETRF_RK" else "ZHETRF_RK", cast(*u8) *uplo, n, -1, -1, -1));
        }
    } else {
        iws = 1;
    }
    if nb < nbmin {
        nb = n;
    }

    work_mat_raw := _hetrf_rk_matrix_storage(work.data, ldwork, ldwork, max(1, nb));
    work_mat := array2d_view(work_mat_raw, ldwork, row_major=false);

    if upper {
        k = n;
        while true {
            if k < 1 then break;

            if k > nb {
                a_block_raw := _hetrf_rk_matrix_view(a.*, lda, 1, 1, k, k);
                a_block := array2d_view(a_block_raw, lda, row_major=false);
                e_block := _hetrf_rk_vector_slice(e, 1, k);
                ipiv_block := _hetrf_rk_vector_slice(ipiv, 1, k);
                _hetrf_rk_lahef(uplo, k, nb, *kb, *a_block, lda, *e_block, *ipiv_block, *work_mat, ldwork, *iinfo, MT);
            } else {
                a_block_raw := _hetrf_rk_matrix_view(a.*, lda, 1, 1, k, k);
                a_block := array2d_view(a_block_raw, lda, row_major=false);
                e_block := _hetrf_rk_vector_slice(e, 1, k);
                ipiv_block := _hetrf_rk_vector_slice(ipiv, 1, k);
                _hetrf_rk_hetf2(uplo, k, *a_block, lda, *e_block, *ipiv_block, *iinfo, MT);
                kb = k;
            }

            if info.* == 0 && iinfo > 0 {
                info.* = iinfo;
            }

            if k < n {
                for <i:k - kb + 1..k {
                    ip = abs(fget(ipiv, i));
                    if ip != i {
                        row_i := _hetrf_rk_matrix_row(a.*, lda, i, k + 1, n - k);
                        row_ip := _hetrf_rk_matrix_row(a.*, lda, ip, k + 1, n - k);
                        _hetrf_rk_swap(n - k, *row_i, lda, *row_ip, lda, MT);
                    }
                }
            }

            k -= kb;
        }
    } else {
        k = 1;
        while true {
            if k > n then break;

            if k <= n - nb {
                a_block_raw := _hetrf_rk_matrix_view(a.*, lda, k, k, n - k + 1, n - k + 1);
                a_block := array2d_view(a_block_raw, lda, row_major=false);
                e_block := _hetrf_rk_vector_slice(e, k, n - k + 1);
                ipiv_block := _hetrf_rk_vector_slice(ipiv, k, n - k + 1);
                _hetrf_rk_lahef(uplo, n - k + 1, nb, *kb, *a_block, lda, *e_block, *ipiv_block, *work_mat, ldwork, *iinfo, MT);
            } else {
                a_block_raw := _hetrf_rk_matrix_view(a.*, lda, k, k, n - k + 1, n - k + 1);
                a_block := array2d_view(a_block_raw, lda, row_major=false);
                e_block := _hetrf_rk_vector_slice(e, k, n - k + 1);
                ipiv_block := _hetrf_rk_vector_slice(ipiv, k, n - k + 1);
                _hetrf_rk_hetf2(uplo, n - k + 1, *a_block, lda, *e_block, *ipiv_block, *iinfo, MT);
                kb = n - k + 1;
            }

            if info.* == 0 && iinfo > 0 {
                info.* = iinfo + k - 1;
            }

            for i:k..k + kb - 1 {
                if fget(ipiv, i) > 0 {
                    fset(ipiv, i, fget(ipiv.*, i) + k - 1);
                } else {
                    fset(ipiv, i, fget(ipiv.*, i) - k + 1);
                }
            }

            if k > 1 {
                for i:k..k + kb - 1 {
                    ip = abs(fget(ipiv, i));
                    if ip != i {
                        row_i := _hetrf_rk_matrix_row(a.*, lda, i, 1, k - 1);
                        row_ip := _hetrf_rk_matrix_row(a.*, lda, ip, 1, k - 1);
                        _hetrf_rk_swap(k - 1, *row_i, lda, *row_ip, lda, MT);
                    }
                }
            }

            k += kb;
        }
    }

    fset(work, 1, _hetrf_rk_complex_from_real(_hetrf_rk_roundup_lwork(ncast(RT, lwkopt), RT), RT, MT));
}
