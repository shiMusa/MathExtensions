// ? LAPACK *hetrs_aa_2stage()
// ? http://www.netlib.org/lapack/explore-html/

// ? original comment chetrs_aa_2stage/zhetrs_aa_2stage
// *HETRS_AA_2STAGE solves a system of linear equations A*X = B for Hermitian matrices.

chetrs_aa_2stage :: #bake_arguments _hetrs_aa_2stage(RT=float32, MT=C32);
zhetrs_aa_2stage :: #bake_arguments _hetrs_aa_2stage(RT=float64, MT=C64);

#scope_file
_hetrs_aa_2stage_raw_slice :: inline (data: *$T, $$count: int) -> [] T {
    view: [] T;
    view.data = data;
    view.count = ifx count < 0 then 0 else count;
    return view;
}

_hetrs_aa_2stage_matrix_storage :: inline (
    a: $A/ArrayType2D, $$lda: int, $$row: int, $$col: int, $$m: int, $$n: int
) -> [] A.data_type {
    count := ifx m <= 0 || n <= 0 then 0 else (n - 1)*lda + m;
    return _hetrs_aa_2stage_raw_slice(a.data + (col - 1)*lda + (row - 1), count);
}

_hetrs_aa_2stage_complex_from_real :: inline (x: $RT, $RT: Type, $MT: Type) -> MT #modify {
    return (RT == float32 || RT == float64) && MT == complex_with_inner_type(RT);
} {
    return cpx(x, zero(RT));
}

_hetrs_aa_2stage_trsm :: (
    $$side: u8, $$uplo: u8, $$transa: u8, $$diag: u8, $$m: int, $$n: int,
    alpha: $MT, a: $A/ArrayType2D, $$lda: int, b: *$B/ArrayType2D, $$ldb: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && B.data_type == MT;
} {
    #if MT == C32 {
        ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    } else {
        ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
    }
}

_hetrs_aa_2stage_laswp :: (
    $$n: int, a: *$A/ArrayType2D, $$lda: int, $$k1: int, $$k2: int,
    ipiv: $IPIV/ArrayType1D, $$incx: int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && A.data_type == MT && IPIV.data_type == int;
} {
    #if MT == C32 {
        claswp(n, a, lda, k1, k2, ipiv, incx);
    } else {
        zlaswp(n, a, lda, k1, k2, ipiv, incx);
    }
}

_hetrs_aa_2stage_gbtrs :: (
    $$trans: u8, $$n: int, $$kl: int, $$ku: int, $$nrhs: int,
    ab: $AB/ArrayType2D, $$ldab: int, ipiv: $IPIV/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, info: *int, $MT: Type
) #modify {
    return (MT == C32 || MT == C64) && AB.data_type == MT && IPIV.data_type == int && B.data_type == MT;
} {
    #if MT == C32 {
        cgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    } else {
        zgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    }
}

_hetrs_aa_2stage :: (
    $$uplo: u8, $$n: int, $$nrhs: int, a: $A/ArrayType2D, $$lda: int,
    tb: $TB/ArrayType1D, $$ltb: int, ipiv: $IPIV/ArrayType1D, ipiv2: $IPIV2/ArrayType1D,
    b: *$B/ArrayType2D, $$ldb: int, info: *int, $RT: Type, $MT: Type
) #modify {
    return (RT == float32 || RT == float64) &&
           MT == complex_with_inner_type(RT) &&
           A.data_type == MT && TB.data_type == MT &&
           IPIV.data_type == int && IPIV2.data_type == int && B.data_type == MT;
} {
    one_t :: #run _hetrs_aa_2stage_complex_from_real(one(RT), RT, MT);
    routine_name :: #run if MT == C32 { "CHETRS_AA_2STAGE" } else { "ZHETRS_AA_2STAGE" };

    ldttb, nb : int = ---;
    upper : bool = ---;

    info.* = 0;
    upper = lsame(uplo, #char "U");
    if !(upper || lsame(uplo, #char "L")) {
        info.* = -1;
    } else if n < 0 {
        info.* = -2;
    } else if nrhs < 0 {
        info.* = -3;
    } else if lda < max(1, n) {
        info.* = -5;
    } else if ltb < 4*n {
        info.* = -7;
    } else if ldb < max(1, n) {
        info.* = -11;
    }
    if info.* != 0 {
        xerbla(routine_name, -info.*);
        return;
    }

    if n == 0 || nrhs == 0 {
        return;
    }

    nb = cast(int) real(fget(tb, 1));
    ldttb = ltb / n;
    tb_raw := _hetrs_aa_2stage_raw_slice(tb.data, ldttb*n);
    tb_mat := array2d_view(tb_raw, ldttb, row_major=false);
    b_full_raw := _hetrs_aa_2stage_matrix_storage(b.*, ldb, 1, 1, n, nrhs);
    b_full := array2d_view(b_full_raw, ldb, row_major=false);

    if upper {
        if n > nb {
            _hetrs_aa_2stage_laswp(nrhs, *b_full, ldb, nb + 1, n, ipiv, 1, MT);

            a_u_raw := _hetrs_aa_2stage_matrix_storage(a, lda, 1, nb + 1, n - nb, n - nb);
            a_u := array2d_view(a_u_raw, lda, row_major=false);
            b_tail_raw := _hetrs_aa_2stage_matrix_storage(b.*, ldb, nb + 1, 1, n - nb, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _hetrs_aa_2stage_trsm(#char "L", #char "U", #char "C", #char "U", n - nb, nrhs, one_t, *a_u, lda, *b_tail, ldb, MT);
        }

        _hetrs_aa_2stage_gbtrs(#char "N", n, nb, nb, nrhs, *tb_mat, ldttb, ipiv2, *b_full, ldb, info, MT);

        if n > nb {
            a_u_raw := _hetrs_aa_2stage_matrix_storage(a, lda, 1, nb + 1, n - nb, n - nb);
            a_u := array2d_view(a_u_raw, lda, row_major=false);
            b_tail_raw := _hetrs_aa_2stage_matrix_storage(b.*, ldb, nb + 1, 1, n - nb, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _hetrs_aa_2stage_trsm(#char "L", #char "U", #char "N", #char "U", n - nb, nrhs, one_t, *a_u, lda, *b_tail, ldb, MT);
            _hetrs_aa_2stage_laswp(nrhs, *b_full, ldb, nb + 1, n, ipiv, -1, MT);
        }
    } else {
        if n > nb {
            _hetrs_aa_2stage_laswp(nrhs, *b_full, ldb, nb + 1, n, ipiv, 1, MT);

            a_l_raw := _hetrs_aa_2stage_matrix_storage(a, lda, nb + 1, 1, n - nb, n - nb);
            a_l := array2d_view(a_l_raw, lda, row_major=false);
            b_tail_raw := _hetrs_aa_2stage_matrix_storage(b.*, ldb, nb + 1, 1, n - nb, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _hetrs_aa_2stage_trsm(#char "L", #char "L", #char "N", #char "U", n - nb, nrhs, one_t, *a_l, lda, *b_tail, ldb, MT);
        }

        _hetrs_aa_2stage_gbtrs(#char "N", n, nb, nb, nrhs, *tb_mat, ldttb, ipiv2, *b_full, ldb, info, MT);

        if n > nb {
            a_l_raw := _hetrs_aa_2stage_matrix_storage(a, lda, nb + 1, 1, n - nb, n - nb);
            a_l := array2d_view(a_l_raw, lda, row_major=false);
            b_tail_raw := _hetrs_aa_2stage_matrix_storage(b.*, ldb, nb + 1, 1, n - nb, nrhs);
            b_tail := array2d_view(b_tail_raw, ldb, row_major=false);
            _hetrs_aa_2stage_trsm(#char "L", #char "L", #char "C", #char "U", n - nb, nrhs, one_t, *a_l, lda, *b_tail, ldb, MT);
            _hetrs_aa_2stage_laswp(nrhs, *b_full, ldb, nb + 1, n, ipiv, -1, MT);
        }
    }
}
